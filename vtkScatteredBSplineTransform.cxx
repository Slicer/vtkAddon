/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

=========================================================================auto=*/

#include <algorithm> // VTK 8.2.0 has bug for C++17 "clamp" function (algorithm must be included before vtMath.h)

#include "vtkScatteredBSplineTransform.h"

#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"

#include <cmath>
#include <vector>

vtkStandardNewMacro(vtkScatteredBSplineTransform);

vtkCxxSetObjectMacro(vtkScatteredBSplineTransform,InputPointSet,vtkPointSet);

namespace
{

// Refinement is stopped if the control point lattice would grow beyond this
// many nodes (each node stores 3 doubles, so this corresponds to about 400MB).
const vtkIdType MAXIMUM_NUMBER_OF_LATTICE_NODES = 16 * 1024 * 1024;

//----------------------------------------------------------------------------
// Values of the four uniform cubic b-spline basis functions that overlap a
// point at fractional position f within a lattice cell.
// This is the same polynomial, in the same order, that vtkBSplineTransform
// uses for evaluation (see vtkBSplineTransformWeights in
// VTK/Filters/Hybrid/vtkBSplineTransform.cxx), therefore the control point
// coefficients computed here can be handed to the parent class as-is.
inline void vtkScatteredBSplineWeights(double f, double w[4])
{
  const double sixth = 1.0/6.0;
  const double f2 = f*f;

  w[3] = f2*f*sixth;
  w[0] = (f2-f)*0.5 - w[3] + sixth;
  w[2] = f + w[0] - w[3]*2;
  w[1] = 1.0 - w[0] - w[2] - w[3];
}

//----------------------------------------------------------------------------
// A regular lattice of cubic b-spline control points. Each node stores a
// 3-component displacement coefficient. Node (i,j,k) is at
// Origin + (i,j,k)*Spacing and the coefficients are stored with the first axis
// varying fastest, which is the memory layout that vtkImageData expects.
class vtkScatteredBSplineLattice
{
public:
  double Origin[3];
  double Spacing[3];
  int Dims[3];
  std::vector<double> Phi;

  vtkScatteredBSplineLattice()
    {
    for (int axis = 0; axis < 3; axis++)
      {
      this->Origin[axis] = 0.0;
      this->Spacing[axis] = 1.0;
      this->Dims[axis] = 0;
      }
    }

  // Set up a lattice that spans [minCoordinates, maxCoordinates] using
  // nodeCounts nodes between the domain corners (inclusive) along each axis,
  // and zero all coefficients.
  // Three extra nodes are allocated along each axis: one below the domain and
  // two above it, so that every point of the domain - including a point
  // exactly at the upper corner - is covered by four in-range nodes.
  void SetGeometry(const double minCoordinates[3], const double maxCoordinates[3],
                   const int nodeCounts[3])
    {
    for (int axis = 0; axis < 3; axis++)
      {
      this->Spacing[axis] = (maxCoordinates[axis]-minCoordinates[axis])/(nodeCounts[axis]-1);
      this->Origin[axis] = minCoordinates[axis]-this->Spacing[axis];
      this->Dims[axis] = nodeCounts[axis]+3;
      }
    this->Phi.assign(3*this->GetNumberOfNodes(), 0.0);
    }

  vtkIdType GetNumberOfNodes() const
    {
    return static_cast<vtkIdType>(this->Dims[0])*this->Dims[1]*this->Dims[2];
    }

  vtkIdType GetNodeIndex(int i, int j, int k) const
    {
    return (static_cast<vtkIdType>(k)*this->Dims[1]+j)*this->Dims[0]+i;
    }

  // Index of the first of the four overlapping nodes along each axis, and the
  // corresponding basis function values.
  void GetWeights(const double point[3], int base[3], double weights[3][4]) const
    {
    for (int axis = 0; axis < 3; axis++)
      {
      const double u = (point[axis]-this->Origin[axis])/this->Spacing[axis];
      const double flooredU = std::floor(u);
      int index = static_cast<int>(flooredU)-1;
      // The lattice is padded so that this clamping only has an effect for
      // points on the domain boundary, where rounding may push the index out
      // of range by one.
      index = std::max(index, 0);
      index = std::min(index, this->Dims[axis]-4);
      base[axis] = index;
      vtkScatteredBSplineWeights(u-flooredU, weights[axis]);
      }
    }

  void Evaluate(const double point[3], double displacement[3]) const
    {
    int base[3] = { 0, 0, 0 };
    double weights[3][4] = { { 0.0 } };
    this->GetWeights(point, base, weights);

    displacement[0] = 0.0;
    displacement[1] = 0.0;
    displacement[2] = 0.0;
    for (int k = 0; k < 4; k++)
      {
      for (int j = 0; j < 4; j++)
        {
        const double weightJK = weights[1][j]*weights[2][k];
        const double* phi = &this->Phi[3*this->GetNodeIndex(base[0], base[1]+j, base[2]+k)];
        for (int i = 0; i < 4; i++)
          {
          const double weight = weights[0][i]*weightJK;
          displacement[0] += weight*phi[0];
          displacement[1] += weight*phi[1];
          displacement[2] += weight*phi[2];
          phi += 3;
          }
        }
      }
    }
};

//----------------------------------------------------------------------------
// Single-level b-spline approximation (the "BA" algorithm of Lee, Wolberg and
// Shin): each point pushes the value it needs onto the 64 control points that
// influence it, and each control point takes the weighted average of what was
// pushed onto it.
// The b-spline weights only depend on the position, so the three displacement
// components share the same denominator and can be fitted in a single pass.
void vtkScatteredBSplineFitLevel(const std::vector<double>& points,
                                 const std::vector<double>& values,
                                 vtkScatteredBSplineLattice& lattice)
{
  const vtkIdType numberOfNodes = lattice.GetNumberOfNodes();
  std::vector<double> delta(3*numberOfNodes, 0.0);
  std::vector<double> omega(numberOfNodes, 0.0);

  const vtkIdType numberOfPoints = static_cast<vtkIdType>(points.size()/3);
  for (vtkIdType pointIndex = 0; pointIndex < numberOfPoints; pointIndex++)
    {
    int base[3] = { 0, 0, 0 };
    double weights[3][4] = { { 0.0 } };
    lattice.GetWeights(&points[3*pointIndex], base, weights);

    double tensorWeights[64];
    double sumOfSquaredWeights = 0.0;
    int weightIndex = 0;
    for (int k = 0; k < 4; k++)
      {
      for (int j = 0; j < 4; j++)
        {
        const double weightJK = weights[1][j]*weights[2][k];
        for (int i = 0; i < 4; i++)
          {
          const double weight = weights[0][i]*weightJK;
          tensorWeights[weightIndex++] = weight;
          sumOfSquaredWeights += weight*weight;
          }
        }
      }
    if (sumOfSquaredWeights <= 0.0)
      {
      continue;
      }

    const double* value = &values[3*pointIndex];
    weightIndex = 0;
    for (int k = 0; k < 4; k++)
      {
      for (int j = 0; j < 4; j++)
        {
        const vtkIdType nodeIndex = lattice.GetNodeIndex(base[0], base[1]+j, base[2]+k);
        double* deltaNode = &delta[3*nodeIndex];
        double* omegaNode = &omega[nodeIndex];
        for (int i = 0; i < 4; i++)
          {
          const double weight = tensorWeights[weightIndex++];
          const double squaredWeight = weight*weight;
          // The value that this point alone would require at this node is
          // value*weight/sumOfSquaredWeights; nodes accumulate those requests
          // weighted by squaredWeight.
          const double scale = squaredWeight*weight/sumOfSquaredWeights;
          deltaNode[0] += scale*value[0];
          deltaNode[1] += scale*value[1];
          deltaNode[2] += scale*value[2];
          *omegaNode += squaredWeight;
          deltaNode += 3;
          omegaNode++;
          }
        }
      }
    }

  for (vtkIdType nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++)
    {
    const double weightSum = omega[nodeIndex];
    for (int component = 0; component < 3; component++)
      {
      lattice.Phi[3*nodeIndex+component] =
        (weightSum == 0.0 ? 0.0 : delta[3*nodeIndex+component]/weightSum);
      }
    }
}

//----------------------------------------------------------------------------
// Add the coarse lattice to the fine lattice, which must have exactly half the
// control point spacing of the coarse one over the same domain.
// Uniform cubic b-splines can be subdivided exactly, so the transform
// represented by the coarse lattice is preserved: the coefficient of each
// coarse node is spread onto the 5x5x5 surrounding fine nodes with the
// separable kernel {1/8, 1/2, 3/4, 1/2, 1/8}.
void vtkScatteredBSplineAddRefined(const vtkScatteredBSplineLattice& coarse,
                                   vtkScatteredBSplineLattice& fine)
{
  static const double kernel[5] = { 0.125, 0.5, 0.75, 0.5, 0.125 };

  for (int k = 0; k < coarse.Dims[2]; k++)
    {
    for (int j = 0; j < coarse.Dims[1]; j++)
      {
      for (int i = 0; i < coarse.Dims[0]; i++)
        {
        const double* coarsePhi = &coarse.Phi[3*coarse.GetNodeIndex(i, j, k)];
        if (coarsePhi[0] == 0.0 && coarsePhi[1] == 0.0 && coarsePhi[2] == 0.0)
          {
          continue;
          }
        for (int dk = 0; dk < 5; dk++)
          {
          const int fineK = 2*k+dk-3;
          if (fineK < 0 || fineK >= fine.Dims[2])
            {
            continue;
            }
          for (int dj = 0; dj < 5; dj++)
            {
            const int fineJ = 2*j+dj-3;
            if (fineJ < 0 || fineJ >= fine.Dims[1])
              {
              continue;
              }
            const double kernelJK = kernel[dj]*kernel[dk];
            for (int di = 0; di < 5; di++)
              {
              const int fineI = 2*i+di-3;
              if (fineI < 0 || fineI >= fine.Dims[0])
                {
                continue;
                }
              const double weight = kernel[di]*kernelJK;
              double* finePhi = &fine.Phi[3*fine.GetNodeIndex(fineI, fineJ, fineK)];
              finePhi[0] += weight*coarsePhi[0];
              finePhi[1] += weight*coarsePhi[1];
              finePhi[2] += weight*coarsePhi[2];
              }
            }
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
// Compute how far the lattice still is from the prescribed displacements and
// return the largest remaining displacement error.
double vtkScatteredBSplineComputeResiduals(const vtkScatteredBSplineLattice& lattice,
                                           const std::vector<double>& points,
                                           const std::vector<double>& displacements,
                                           std::vector<double>& residuals)
{
  double maximumResidual = 0.0;
  const vtkIdType numberOfPoints = static_cast<vtkIdType>(points.size()/3);
  for (vtkIdType pointIndex = 0; pointIndex < numberOfPoints; pointIndex++)
    {
    double fitted[3] = { 0.0, 0.0, 0.0 };
    lattice.Evaluate(&points[3*pointIndex], fitted);
    double* residual = &residuals[3*pointIndex];
    for (int component = 0; component < 3; component++)
      {
      residual[component] = displacements[3*pointIndex+component]-fitted[component];
      }
    maximumResidual = std::max(maximumResidual, vtkMath::Norm(residual));
    }
  return maximumResidual;
}

//----------------------------------------------------------------------------
// Least-squares fit of a linear function to the displacements:
// displacement[component] = coefficients[component][0..2]*position + coefficients[component][3]
// Falls back to the mean displacement if the normal equations are singular
// (for example, if all the points are on a plane or there are too few of them).
void vtkScatteredBSplineLinearApproximation(const std::vector<double>& points,
                                            const std::vector<double>& displacements,
                                            double coefficients[3][4])
{
  const vtkIdType numberOfPoints = static_cast<vtkIdType>(points.size()/3);

  // Normal equations of the overdetermined system: normalMatrix = X'*X, rightHandSide = X'*d
  double normalMatrix[4][4] = { { 0.0 } };
  double rightHandSide[3][4] = { { 0.0 } };
  for (vtkIdType pointIndex = 0; pointIndex < numberOfPoints; pointIndex++)
    {
    const double* point = &points[3*pointIndex];
    const double* displacement = &displacements[3*pointIndex];
    const double row[4] = { point[0], point[1], point[2], 1.0 };
    for (int i = 0; i < 4; i++)
      {
      for (int j = 0; j < 4; j++)
        {
        normalMatrix[i][j] += row[i]*row[j];
        }
      for (int component = 0; component < 3; component++)
        {
        rightHandSide[component][i] += row[i]*displacement[component];
        }
      }
    }

  // At least 4 points are needed to define a linear function in 3D
  bool solved = (numberOfPoints >= 4);
  for (int component = 0; component < 3 && solved; component++)
    {
    // vtkMath::SolveLinearSystem overwrites both the matrix and the vector
    double matrix[4][4];
    double* matrixRows[4] = { matrix[0], matrix[1], matrix[2], matrix[3] };
    std::copy(&normalMatrix[0][0], &normalMatrix[0][0]+16, &matrix[0][0]);
    std::copy(rightHandSide[component], rightHandSide[component]+4, coefficients[component]);
    if (!vtkMath::SolveLinearSystem(matrixRows, coefficients[component], 4))
      {
      solved = false;
      }
    }

  if (!solved)
    {
    // Not enough points or a singular system: use the mean displacement
    for (int component = 0; component < 3; component++)
      {
      coefficients[component][0] = 0.0;
      coefficients[component][1] = 0.0;
      coefficients[component][2] = 0.0;
      coefficients[component][3] = (numberOfPoints > 0 ?
        rightHandSide[component][3]/numberOfPoints : 0.0);
      }
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
vtkScatteredBSplineTransform::vtkScatteredBSplineTransform()
{
  this->InputPointSet = nullptr;
  this->DisplacementArrayName = nullptr;

  this->ControlPointSpacing[0] = 0.0;
  this->ControlPointSpacing[1] = 0.0;
  this->ControlPointSpacing[2] = 0.0;
  this->MinimumControlPointSpacing = 0.1;
  this->MaximumNumberOfLevels = 8;
  this->Tolerance = 0.01;
  this->UseLinearApproximation = false;

  this->GridBounds[0] = 0.0;
  this->GridBounds[1] = -1.0;
  this->GridBounds[2] = 0.0;
  this->GridBounds[3] = -1.0;
  this->GridBounds[4] = 0.0;
  this->GridBounds[5] = -1.0;
  this->AutomaticGridBounds = true;
  this->GridBoundsMargin = 0.0;

  this->Residual = -1.0;
  this->NumberOfLevels = 0;
}

//----------------------------------------------------------------------------
vtkScatteredBSplineTransform::~vtkScatteredBSplineTransform()
{
  this->SetInputPointSet(nullptr);
  this->SetDisplacementArrayName(nullptr);
}

//----------------------------------------------------------------------------
void vtkScatteredBSplineTransform::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "InputPointSet: " << this->InputPointSet << "\n";
  os << indent << "DisplacementArrayName: "
     << (this->DisplacementArrayName ? this->DisplacementArrayName : "(none)") << "\n";
  os << indent << "ControlPointSpacing: " << this->ControlPointSpacing[0] << " "
     << this->ControlPointSpacing[1] << " " << this->ControlPointSpacing[2] << "\n";
  os << indent << "MinimumControlPointSpacing: " << this->MinimumControlPointSpacing << "\n";
  os << indent << "MaximumNumberOfLevels: " << this->MaximumNumberOfLevels << "\n";
  os << indent << "Tolerance: " << this->Tolerance << "\n";
  os << indent << "UseLinearApproximation: " << (this->UseLinearApproximation ? "On" : "Off") << "\n";
  os << indent << "AutomaticGridBounds: " << (this->AutomaticGridBounds ? "On" : "Off") << "\n";
  os << indent << "GridBounds: " << this->GridBounds[0] << " " << this->GridBounds[1] << " "
     << this->GridBounds[2] << " " << this->GridBounds[3] << " "
     << this->GridBounds[4] << " " << this->GridBounds[5] << "\n";
  os << indent << "GridBoundsMargin: " << this->GridBoundsMargin << "\n";
  os << indent << "Residual: " << this->Residual << "\n";
  os << indent << "NumberOfLevels: " << this->NumberOfLevels << "\n";
}

//----------------------------------------------------------------------------
vtkAbstractTransform* vtkScatteredBSplineTransform::MakeTransform()
{
  return vtkScatteredBSplineTransform::New();
}

//----------------------------------------------------------------------------
vtkMTimeType vtkScatteredBSplineTransform::GetMTime()
{
  vtkMTimeType result = this->Superclass::GetMTime();
  if (this->InputPointSet)
    {
    result = std::max(result, this->InputPointSet->GetMTime());
    }
  return result;
}

//----------------------------------------------------------------------------
void vtkScatteredBSplineTransform::InternalUpdate()
{
  this->FitCoefficients();
  // The parent class picks up the fitted coefficients that FitCoefficients()
  // has installed as the coefficient data.
  this->Superclass::InternalUpdate();
}

//----------------------------------------------------------------------------
void vtkScatteredBSplineTransform::ClearCoefficients()
{
  this->SetCoefficientConnection(nullptr);
  // vtkBSplineTransform leaves its grid geometry untouched when there is no
  // coefficient data, so reset it here to keep the matrices that
  // vtkOrientedBSplineTransform::InternalUpdate() computes from them
  // well-defined
  for (int axis = 0; axis < 3; axis++)
    {
    this->Superclass::GridOrigin[axis] = 0.0;
    this->Superclass::GridSpacing[axis] = 1.0;
    this->Superclass::GridExtent[2*axis] = 0;
    this->Superclass::GridExtent[2*axis+1] = 0;
    this->Superclass::GridIncrements[axis] = 0;
    }
}

//----------------------------------------------------------------------------
bool vtkScatteredBSplineTransform::FitCoefficients()
{
  this->Residual = -1.0;
  this->NumberOfLevels = 0;

  if (!this->InputPointSet)
    {
    // Nothing to fit yet: the identity
    this->ClearCoefficients();
    return false;
    }

  vtkPoints* inputPoints = this->InputPointSet->GetPoints();
  const vtkIdType numberOfInputPoints = (inputPoints ? inputPoints->GetNumberOfPoints() : 0);
  if (numberOfInputPoints < 1)
    {
    vtkErrorMacro("vtkScatteredBSplineTransform: input point set is empty");
    this->ClearCoefficients();
    return false;
    }

  vtkDataArray* displacementArray = nullptr;
  vtkPointData* pointData = this->InputPointSet->GetPointData();
  if (this->DisplacementArrayName && this->DisplacementArrayName[0] != '\0')
    {
    displacementArray = pointData->GetArray(this->DisplacementArrayName);
    if (!displacementArray)
      {
      vtkErrorMacro("vtkScatteredBSplineTransform: input point set has no point data array named "
        << this->DisplacementArrayName);
      this->ClearCoefficients();
      return false;
      }
    }
  else
    {
    displacementArray = pointData->GetVectors();
    if (!displacementArray)
      {
      vtkErrorMacro("vtkScatteredBSplineTransform: DisplacementArrayName is not set"
        " and the input point set has no active vectors array");
      this->ClearCoefficients();
      return false;
      }
    }
  if (displacementArray->GetNumberOfComponents() != 3)
    {
    vtkErrorMacro("vtkScatteredBSplineTransform: displacement array "
      << displacementArray->GetName() << " must have 3 components, it has "
      << displacementArray->GetNumberOfComponents());
    this->ClearCoefficients();
    return false;
    }
  if (displacementArray->GetNumberOfTuples() != numberOfInputPoints)
    {
    vtkErrorMacro("vtkScatteredBSplineTransform: displacement array "
      << displacementArray->GetName() << " has " << displacementArray->GetNumberOfTuples()
      << " tuples, while the input point set has " << numberOfInputPoints << " points");
    this->ClearCoefficients();
    return false;
    }

  // Determine the fitting domain
  double minCoordinates[3] = { 0.0, 0.0, 0.0 };
  double maxCoordinates[3] = { 0.0, 0.0, 0.0 };
  if (this->AutomaticGridBounds)
    {
    double bounds[6] = { 0.0, -1.0, 0.0, -1.0, 0.0, -1.0 };
    inputPoints->GetBounds(bounds);
    for (int axis = 0; axis < 3; axis++)
      {
      minCoordinates[axis] = bounds[2*axis] - this->GridBoundsMargin;
      maxCoordinates[axis] = bounds[2*axis+1] + this->GridBoundsMargin;
      }
    }
  else
    {
    for (int axis = 0; axis < 3; axis++)
      {
      minCoordinates[axis] = this->GridBounds[2*axis];
      maxCoordinates[axis] = this->GridBounds[2*axis+1];
      if (minCoordinates[axis] > maxCoordinates[axis])
        {
        vtkErrorMacro("vtkScatteredBSplineTransform: invalid GridBounds, the minimum coordinate"
          " is larger than the maximum coordinate along axis " << axis);
        this->ClearCoefficients();
        return false;
        }
      }
    }

  double largestExtent = 0.0;
  for (int axis = 0; axis < 3; axis++)
    {
    largestExtent = std::max(largestExtent, maxCoordinates[axis]-minCoordinates[axis]);
    }
  if (largestExtent <= 0.0)
    {
    vtkErrorMacro("vtkScatteredBSplineTransform: the fitting domain is empty"
      " (all input points are at the same position)");
    this->ClearCoefficients();
    return false;
    }

  // Number of control points of the coarsest level
  int nodeCounts[3] = { 2, 2, 2 };
  for (int axis = 0; axis < 3; axis++)
    {
    double extent = maxCoordinates[axis]-minCoordinates[axis];
    double spacing = this->ControlPointSpacing[axis];
    if (spacing <= 0.0)
      {
      // Default to one tenth of the domain size. If the domain is flat along
      // this axis then use the size of the largest axis instead.
      spacing = (extent > 0.0 ? extent : largestExtent)/10.0;
      }
    if (extent < spacing)
      {
      // Flat or very thin domain: widen it so that the transform varies
      // smoothly along this axis, too
      const double center = 0.5*(minCoordinates[axis]+maxCoordinates[axis]);
      minCoordinates[axis] = center-0.5*spacing;
      maxCoordinates[axis] = center+0.5*spacing;
      extent = spacing;
      }
    // Computed in double to avoid overflow if a very small spacing is requested
    const double nodeCount = 1.0+std::floor(extent/spacing);
    nodeCounts[axis] = static_cast<int>(std::min(std::max(nodeCount, 2.0),
      static_cast<double>(MAXIMUM_NUMBER_OF_LATTICE_NODES)));
    }

  double numberOfNodes = 1.0;
  for (int axis = 0; axis < 3; axis++)
    {
    numberOfNodes *= nodeCounts[axis]+3;
    }
  if (numberOfNodes > MAXIMUM_NUMBER_OF_LATTICE_NODES)
    {
    vtkErrorMacro("vtkScatteredBSplineTransform: the requested control point spacing would need "
      << numberOfNodes << " control points, which is more than the " << MAXIMUM_NUMBER_OF_LATTICE_NODES
      << " limit. Increase ControlPointSpacing or reduce the size of the fitting domain.");
    this->ClearCoefficients();
    return false;
    }

  // Collect the points that are inside the fitting domain
  std::vector<double> points;
  std::vector<double> displacements;
  points.reserve(3*numberOfInputPoints);
  displacements.reserve(3*numberOfInputPoints);
  vtkIdType numberOfSkippedPoints = 0;
  for (vtkIdType pointIndex = 0; pointIndex < numberOfInputPoints; pointIndex++)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    inputPoints->GetPoint(pointIndex, point);
    bool inside = true;
    for (int axis = 0; axis < 3; axis++)
      {
      if (point[axis] < minCoordinates[axis] || point[axis] > maxCoordinates[axis])
        {
        inside = false;
        break;
        }
      }
    if (!inside)
      {
      numberOfSkippedPoints++;
      continue;
      }
    double displacement[3] = { 0.0, 0.0, 0.0 };
    displacementArray->GetTuple(pointIndex, displacement);
    points.insert(points.end(), point, point+3);
    displacements.insert(displacements.end(), displacement, displacement+3);
    }
  if (numberOfSkippedPoints > 0)
    {
    vtkWarningMacro("vtkScatteredBSplineTransform: " << numberOfSkippedPoints << " of "
      << numberOfInputPoints << " input points are outside the fitting domain,"
      " they are ignored");
    }
  if (points.empty())
    {
    vtkErrorMacro("vtkScatteredBSplineTransform: no input points are inside the fitting domain");
    this->ClearCoefficients();
    return false;
    }

  // Coarsest level
  vtkScatteredBSplineLattice lattice;
  lattice.SetGeometry(minCoordinates, maxCoordinates, nodeCounts);
  std::vector<double> residuals(displacements.size(), 0.0);

  if (this->UseLinearApproximation)
    {
    // Seed the lattice with the overall linear trend of the displacements.
    // Uniform cubic b-splines reproduce linear functions exactly, so sampling
    // the linear function at the control points is enough.
    double linearCoefficients[3][4] = { { 0.0 } };
    vtkScatteredBSplineLinearApproximation(points, displacements, linearCoefficients);
    for (int k = 0; k < lattice.Dims[2]; k++)
      {
      for (int j = 0; j < lattice.Dims[1]; j++)
        {
        for (int i = 0; i < lattice.Dims[0]; i++)
          {
          const double nodePosition[3] = {
            lattice.Origin[0]+i*lattice.Spacing[0],
            lattice.Origin[1]+j*lattice.Spacing[1],
            lattice.Origin[2]+k*lattice.Spacing[2] };
          double* phi = &lattice.Phi[3*lattice.GetNodeIndex(i, j, k)];
          for (int component = 0; component < 3; component++)
            {
            phi[component] = linearCoefficients[component][0]*nodePosition[0]
                           + linearCoefficients[component][1]*nodePosition[1]
                           + linearCoefficients[component][2]*nodePosition[2]
                           + linearCoefficients[component][3];
            }
          }
        }
      }
    vtkScatteredBSplineComputeResiduals(lattice, points, displacements, residuals);
    // Fit what is left over and add it to the seed
    vtkScatteredBSplineLattice correction;
    correction.SetGeometry(minCoordinates, maxCoordinates, nodeCounts);
    vtkScatteredBSplineFitLevel(points, residuals, correction);
    const vtkIdType numberOfCoefficients = static_cast<vtkIdType>(lattice.Phi.size());
    for (vtkIdType coefficientIndex = 0; coefficientIndex < numberOfCoefficients; coefficientIndex++)
      {
      lattice.Phi[coefficientIndex] += correction.Phi[coefficientIndex];
      }
    }
  else
    {
    vtkScatteredBSplineFitLevel(points, displacements, lattice);
    }
  double residual = vtkScatteredBSplineComputeResiduals(lattice, points, displacements, residuals);
  int numberOfLevels = 1;

  // Refinement levels
  while (numberOfLevels < this->MaximumNumberOfLevels && residual > this->Tolerance)
    {
    if (this->MinimumControlPointSpacing > 0.0)
      {
      bool spacingLimitReached = false;
      for (int axis = 0; axis < 3; axis++)
        {
        if (0.5*lattice.Spacing[axis] < this->MinimumControlPointSpacing)
          {
          spacingLimitReached = true;
          break;
          }
        }
      if (spacingLimitReached)
        {
        break;
        }
      }

    int refinedNodeCounts[3] = { 0, 0, 0 };
    double numberOfRefinedNodes = 1.0;
    for (int axis = 0; axis < 3; axis++)
      {
      refinedNodeCounts[axis] = 2*nodeCounts[axis]-1;
      numberOfRefinedNodes *= refinedNodeCounts[axis]+3;
      }
    if (numberOfRefinedNodes > MAXIMUM_NUMBER_OF_LATTICE_NODES)
      {
      vtkWarningMacro("vtkScatteredBSplineTransform: refinement stopped after " << numberOfLevels
        << " levels because the control point lattice would exceed "
        << MAXIMUM_NUMBER_OF_LATTICE_NODES << " nodes. Residual: " << residual);
      break;
      }

    vtkScatteredBSplineLattice refinedLattice;
    refinedLattice.SetGeometry(minCoordinates, maxCoordinates, refinedNodeCounts);
    // Fit the remaining error at the finer resolution, then fold in everything
    // that the coarser levels have already accounted for
    vtkScatteredBSplineFitLevel(points, residuals, refinedLattice);
    vtkScatteredBSplineAddRefined(lattice, refinedLattice);

    lattice.Phi.swap(refinedLattice.Phi);
    for (int axis = 0; axis < 3; axis++)
      {
      lattice.Origin[axis] = refinedLattice.Origin[axis];
      lattice.Spacing[axis] = refinedLattice.Spacing[axis];
      lattice.Dims[axis] = refinedLattice.Dims[axis];
      nodeCounts[axis] = refinedNodeCounts[axis];
      }

    residual = vtkScatteredBSplineComputeResiduals(lattice, points, displacements, residuals);
    numberOfLevels++;
    }

  vtkNew<vtkImageData> coefficients;
  coefficients->SetOrigin(lattice.Origin);
  coefficients->SetSpacing(lattice.Spacing);
  coefficients->SetDimensions(lattice.Dims);
  coefficients->AllocateScalars(VTK_DOUBLE, 3);
  std::copy(lattice.Phi.begin(), lattice.Phi.end(),
    static_cast<double*>(coefficients->GetScalarPointer()));
  this->SetCoefficientData(coefficients);

  this->Residual = residual;
  this->NumberOfLevels = numberOfLevels;
  return true;
}

//----------------------------------------------------------------------------
void vtkScatteredBSplineTransform::InternalDeepCopy(vtkAbstractTransform *transform)
{
  vtkScatteredBSplineTransform *scatteredBSplineTransform =
    static_cast<vtkScatteredBSplineTransform*>(transform);

  // The input point set is never modified by this class, so it is enough to
  // share it instead of copying it
  this->SetInputPointSet(scatteredBSplineTransform->GetInputPointSet());
  this->SetDisplacementArrayName(scatteredBSplineTransform->GetDisplacementArrayName());
  this->SetControlPointSpacing(scatteredBSplineTransform->GetControlPointSpacing());
  this->SetMinimumControlPointSpacing(scatteredBSplineTransform->GetMinimumControlPointSpacing());
  this->SetMaximumNumberOfLevels(scatteredBSplineTransform->GetMaximumNumberOfLevels());
  this->SetTolerance(scatteredBSplineTransform->GetTolerance());
  this->SetUseLinearApproximation(scatteredBSplineTransform->GetUseLinearApproximation());
  this->SetGridBounds(scatteredBSplineTransform->GetGridBounds());
  this->SetAutomaticGridBounds(scatteredBSplineTransform->GetAutomaticGridBounds());
  this->SetGridBoundsMargin(scatteredBSplineTransform->GetGridBoundsMargin());
  this->Residual = scatteredBSplineTransform->Residual;
  this->NumberOfLevels = scatteredBSplineTransform->NumberOfLevels;

  this->Superclass::InternalDeepCopy(transform);
}
