/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Laboratory for Percutaneous Surgery (PerkLab)
  Queen's University, Kingston, ON, Canada. All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was developed with funding from the Children's Hospital of
  Philadelphia Frontier Program.

==============================================================================*/

// vtkAddon includes
#include "vtkMeshDisplacementTransform.h"

// VTK includes
#include <vtkDataArray.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGenericCell.h>
#include <vtkIdTypeArray.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkStaticCellLocator.h>
#include <vtkStaticPointLocator.h>

// STD includes
#include <algorithm>
#include <cmath>
#include <vector>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkMeshDisplacementTransform);
vtkCxxSetObjectMacro(vtkMeshDisplacementTransform, InputPointSet, vtkPointSet);

//----------------------------------------------------------------------------
vtkMeshDisplacementTransform::vtkMeshDisplacementTransform() = default;

//----------------------------------------------------------------------------
vtkMeshDisplacementTransform::~vtkMeshDisplacementTransform()
{
  this->SetInputPointSet(nullptr);
  this->SetDisplacementArrayName(nullptr);
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "InputPointSet: " << this->InputPointSet << "\n";
  os << indent << "DisplacementArrayName: " << (this->DisplacementArrayName ? this->DisplacementArrayName : "(active vectors)") << "\n";
  os << indent << "DisplacementScale: " << this->DisplacementScale << "\n";
  os << indent << "BorderFalloffDistance: " << this->BorderFalloffDistance << "\n";
  os << indent << "AutomaticBorderFalloffDistance: " << (this->AutomaticBorderFalloffDistance ? "On" : "Off") << "\n";
  os << indent << "EffectiveBorderFalloffDistance: " << this->EffectiveBorderFalloffDistance << "\n";
  os << indent << "LogAllInverseConvergenceFailures: " << (this->LogAllInverseConvergenceFailures ? "On" : "Off") << "\n";
  os << indent << "NumberOfInverseConvergenceFailures: " << this->NumberOfInverseConvergenceFailures.load() << "\n";
}

//----------------------------------------------------------------------------
double vtkMeshDisplacementTransform::GetEffectiveBorderFalloffDistance()
{
  return this->EffectiveBorderFalloffDistance;
}

//----------------------------------------------------------------------------
vtkAbstractTransform* vtkMeshDisplacementTransform::MakeTransform()
{
  return vtkMeshDisplacementTransform::New();
}

//----------------------------------------------------------------------------
vtkMTimeType vtkMeshDisplacementTransform::GetMTime()
{
  vtkMTimeType mTime = this->Superclass::GetMTime();
  if (this->InputPointSet)
  {
    mTime = std::max(mTime, this->InputPointSet->GetMTime());
  }
  return mTime;
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::InternalDeepCopy(vtkAbstractTransform* transform)
{
  this->Superclass::InternalDeepCopy(transform);
  vtkMeshDisplacementTransform* meshTransform = static_cast<vtkMeshDisplacementTransform*>(transform);
  // The mesh is never modified by this class, so it is shared
  this->SetInputPointSet(meshTransform->GetInputPointSet());
  this->SetDisplacementArrayName(meshTransform->GetDisplacementArrayName());
  this->SetDisplacementScale(meshTransform->GetDisplacementScale());
  this->SetBorderFalloffDistance(meshTransform->GetBorderFalloffDistance());
  this->SetAutomaticBorderFalloffDistance(meshTransform->GetAutomaticBorderFalloffDistance());
  this->SetLogAllInverseConvergenceFailures(meshTransform->GetLogAllInverseConvergenceFailures());
  // vtkWarpTransform does not copy these. The inverse direction must be copied too, as an inverse transform is
  // updated by copying its forward transform and then inverting itself.
  this->SetInverseTolerance(meshTransform->GetInverseTolerance());
  this->SetInverseIterations(meshTransform->GetInverseIterations());
  if (this->InverseFlag != meshTransform->InverseFlag)
  {
    this->InverseFlag = meshTransform->InverseFlag;
    this->Modified();
  }
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::InternalUpdate()
{
  this->Displacements = nullptr;
  this->PointLocator = nullptr;
  this->CellLocator = nullptr;
  this->BoundarySurface = nullptr;
  this->BoundaryLocator = nullptr;
  this->BoundaryPointIds = nullptr;
  this->EffectiveBorderFalloffDistance = 0.0;
  this->NumberOfInverseConvergenceFailures = 0;
  vtkPointSet* mesh = this->InputPointSet;
  if (!mesh || mesh->GetNumberOfPoints() == 0)
  {
    return;
  }
  vtkDataArray* displacements = (this->DisplacementArrayName && this->DisplacementArrayName[0] != '\0' ? mesh->GetPointData()->GetArray(this->DisplacementArrayName) //
                                                                                                       : mesh->GetPointData()->GetVectors());
  if (!displacements || displacements->GetNumberOfComponents() != 3)
  {
    return;
  }
  this->Displacements = displacements;
  mesh->GetBounds(this->MeshBounds);

  // The falloff distance that keeps the transform invertible follows the size of the displacements
  if (this->AutomaticBorderFalloffDistance)
  {
    double largestDisplacement2 = 0.0;
    for (vtkIdType pointId = 0; pointId < displacements->GetNumberOfTuples(); ++pointId)
    {
      double vector[3] = { 0.0, 0.0, 0.0 };
      displacements->GetTuple(pointId, vector);
      largestDisplacement2 = std::max(largestDisplacement2, vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    }
    this->EffectiveBorderFalloffDistance = 2.0 * std::abs(this->DisplacementScale) * std::sqrt(largestDisplacement2);
  }
  else
  {
    this->EffectiveBorderFalloffDistance = std::max(0.0, this->BorderFalloffDistance);
  }

  // The diagonal of the bounding box of the mesh
  double diagonal = mesh->GetLength();
  if (diagonal <= 0.0)
  {
    diagonal = 1.0;
  }
  // A mesh point that is transformed has exactly the coordinates of the point; the tolerance
  // only absorbs rounding
  this->PointTolerance2 = (1e-9 * diagonal) * (1e-9 * diagonal);
  this->CellTolerance2 = (1e-6 * diagonal) * (1e-6 * diagonal);
  this->DerivativeStep = 1e-4 * diagonal;

  this->PointLocator = vtkSmartPointer<vtkStaticPointLocator>::New();
  this->PointLocator->SetDataSet(mesh);
  this->PointLocator->BuildLocator();
  if (mesh->GetNumberOfCells() > 0)
  {
    this->MaximumCellSize = mesh->GetMaxCellSize();
    this->CellLocator = vtkSmartPointer<vtkStaticCellLocator>::New();
    this->CellLocator->SetDataSet(mesh);
    this->CellLocator->BuildLocator();

    // The boundary surface, for the closest point in the falloff zone
    if (this->EffectiveBorderFalloffDistance > 0.0)
    {
      vtkPolyData* polyData = vtkPolyData::SafeDownCast(mesh);
      if (polyData)
      {
        this->BoundarySurface = polyData;
      }
      else
      {
        const char originalPointIdsName[] = "vtkMeshDisplacementTransformPointIds";
        vtkNew<vtkDataSetSurfaceFilter> surfaceFilter;
        surfaceFilter->SetInputData(mesh);
        // The corner points of the faces only, so that every surface point is a mesh point
        surfaceFilter->SetNonlinearSubdivisionLevel(0);
        surfaceFilter->PassThroughPointIdsOn();
        surfaceFilter->SetOriginalPointIdsName(originalPointIdsName);
        surfaceFilter->Update();
        this->BoundarySurface = surfaceFilter->GetOutput();
        this->BoundaryPointIds = vtkIdTypeArray::SafeDownCast(this->BoundarySurface->GetPointData()->GetArray(originalPointIdsName));
        if (!this->BoundaryPointIds)
        {
          this->BoundarySurface = nullptr;
        }
      }
      if (this->BoundarySurface && this->BoundarySurface->GetNumberOfCells() > 0)
      {
        this->MaximumBoundaryCellSize = this->BoundarySurface->GetMaxCellSize();
        this->BoundaryLocator = vtkSmartPointer<vtkStaticCellLocator>::New();
        this->BoundaryLocator->SetDataSet(this->BoundarySurface);
        this->BoundaryLocator->BuildLocator();
      }
    }
  }
}

//----------------------------------------------------------------------------
bool vtkMeshDisplacementTransform::GetBorderFalloffDisplacement(const double position[3], double borderFalloffDistance, double displacement[3])
{
  if (!this->BoundaryLocator)
  {
    return false;
  }
  vtkPolyData* surface = this->BoundarySurface;
  vtkNew<vtkGenericCell> face;
  double x[3] = { position[0], position[1], position[2] };
  double closestPoint[3] = { 0.0, 0.0, 0.0 };
  vtkIdType faceId = -1;
  int subId = 0;
  double distance2 = 0.0;
  int inside = 0;
  if (!this->BoundaryLocator->FindClosestPointWithinRadius(x, borderFalloffDistance, closestPoint, face, faceId, subId, distance2, inside))
  {
    return false;
  }
  const double t = std::sqrt(distance2) / borderFalloffDistance;
  if (t >= 1.0)
  {
    return false;
  }
  // Falls from 1 to 0 with zero slope at both ends. This was chosen with vtkMeshDisplacementTransformTest3, which
  // measures the inverse consistency (round trip through the inverse) and the smoothness (gradient, and its change
  // along rays leaving the mesh) of each border treatment on a convex, a non-convex and a tetrahedral mesh with a
  // non-linear field. Alternatives that were tried and dropped, so that the class offers one treatment:
  // - A jump to zero at the border (no falloff): the inverse fails at most positions near the border. Kept as the
  //   behavior at distance 0 only because it is exact at the mesh and costs nothing.
  // - Other falloff profiles. Linear (1 - t) inverts from the smallest distance (the largest displacement, against
  //   1.5 times for the cubic), but its slope jumps at both ends, so the gradient changes 3 times faster than with
  //   the cubic. The cosine (1 + cos(pi t)) / 2 is within 5 % of the cubic in every measure. The quintic (zero
  //   curvature at both ends as well) needs 1.875 times the largest displacement to invert and is no smoother in
  //   practice.
  // - Keeping the displacement of the closest point constant everywhere outside (VTK_DOUBLE_MAX distance, as
  //   vtkBSplineTransform does in Edge border mode). The smoothest on a convex mesh, but at a concave corner the
  //   closest point switches between the two inner faces and the displacement jumps there, which folds the
  //   transform: the inverse failed at 20 of 2664 positions around an L-shaped prism. It remains available by
  //   setting the distance to VTK_DOUBLE_MAX, for convex meshes.
  // - Continuing the field of the closest cell to first order (its gradient times the offset from the closest
  //   point). It removes the kink at the border between the field inside and the constant displacement outside,
  //   but tetrahedra have a constant gradient each, which jumps between cells, so the continued displacement jumps
  //   across the extension of every cell boundary: with it, the inverse failed at 113 positions around a
  //   tetrahedral cube.
  const double falloffWeight = (1.0 - t) * (1.0 - t) * (1.0 + 2.0 * t);

  // Interpolation weights at the closest point, in the closest face. The locator leaves the last face that it
  // tested in the generic cell, not the closest one, so the closest face is taken again.
  surface->GetCell(faceId, face);
  std::vector<double> weights(static_cast<size_t>(std::max(1, this->MaximumBoundaryCellSize)));
  double pcoords[3] = { 0.0, 0.0, 0.0 };
  double closestPointInFace[3] = { 0.0, 0.0, 0.0 };
  double faceDistance2 = 0.0;
  face->EvaluatePosition(closestPoint, closestPointInFace, subId, pcoords, faceDistance2, weights.data());
  // The closest point is on the face, so the weights are a convex combination up to rounding. They are made one, so
  // that a position outside the mesh never moves farther than the longest vector of the mesh.
  const vtkIdType numberOfFacePoints = face->GetNumberOfPoints();
  double sumOfWeights = 0.0;
  for (vtkIdType facePointIndex = 0; facePointIndex < numberOfFacePoints; ++facePointIndex)
  {
    weights[facePointIndex] = std::min(1.0, std::max(0.0, weights[facePointIndex]));
    sumOfWeights += weights[facePointIndex];
  }
  if (sumOfWeights <= 0.0)
  {
    return false;
  }
  vtkDataArray* displacements = this->Displacements;
  for (vtkIdType facePointIndex = 0; facePointIndex < numberOfFacePoints; ++facePointIndex)
  {
    vtkIdType pointId = face->GetPointId(facePointIndex);
    if (this->BoundaryPointIds)
    {
      pointId = this->BoundaryPointIds->GetValue(pointId);
    }
    double vector[3] = { 0.0, 0.0, 0.0 };
    displacements->GetTuple(pointId, vector);
    for (int axis = 0; axis < 3; ++axis)
    {
      displacement[axis] += falloffWeight * weights[facePointIndex] / sumOfWeights * vector[axis];
    }
  }
  return true;
}

//----------------------------------------------------------------------------
bool vtkMeshDisplacementTransform::GetDisplacement(const double position[3], double displacement[3])
{
  displacement[0] = displacement[1] = displacement[2] = 0.0;
  vtkDataArray* displacements = this->Displacements;
  vtkPointSet* mesh = this->InputPointSet;
  if (!displacements || !mesh || !this->PointLocator)
  {
    return false;
  }

  // Positions farther from the bounding box of the mesh than the falloff distance do not move
  const double borderFalloffDistance = (this->CellLocator ? this->EffectiveBorderFalloffDistance : 0.0);
  const double margin = borderFalloffDistance + std::sqrt(this->CellTolerance2);
  for (int axis = 0; axis < 3; ++axis)
  {
    if (position[axis] < this->MeshBounds[2 * axis] - margin || position[axis] > this->MeshBounds[2 * axis + 1] + margin)
    {
      return false;
    }
  }

  // A point of the mesh: its own vector, as a warp would use it
  double distance2 = 0.0;
  vtkIdType pointId = this->PointLocator->FindClosestPointWithinRadius(std::sqrt(this->PointTolerance2), position, distance2);
  if (pointId >= 0)
  {
    displacements->GetTuple(pointId, displacement);
  }
  else
  {
    // Anywhere else: interpolated in the cell that contains the position. Transforms may be
    // evaluated from several threads at once, so the scratch space belongs to this call.
    if (!this->CellLocator)
    {
      return false;
    }
    vtkNew<vtkGenericCell> cell;
    std::vector<double> weights(static_cast<size_t>(std::max(1, this->MaximumCellSize)));
    int subId = 0;
    double pcoords[3] = { 0.0, 0.0, 0.0 };
    double x[3] = { position[0], position[1], position[2] };
    vtkIdType cellId = this->CellLocator->FindCell(x, this->CellTolerance2, cell, subId, pcoords, weights.data());
    if (cellId < 0)
    {
      // Outside the cells: the displacement at the closest point of the mesh, weighted by the distance
      if (borderFalloffDistance <= 0.0 || !this->GetBorderFalloffDisplacement(position, borderFalloffDistance, displacement))
      {
        return false;
      }
    }
    else
    {
      vtkIdList* cellPointIds = cell->GetPointIds();
      for (vtkIdType cellPointIndex = 0; cellPointIndex < cellPointIds->GetNumberOfIds(); ++cellPointIndex)
      {
        double vector[3] = { 0.0, 0.0, 0.0 };
        displacements->GetTuple(cellPointIds->GetId(cellPointIndex), vector);
        for (int axis = 0; axis < 3; ++axis)
        {
          displacement[axis] += weights[cellPointIndex] * vector[axis];
        }
      }
    }
  }
  for (int axis = 0; axis < 3; ++axis)
  {
    displacement[axis] *= this->DisplacementScale;
  }
  return true;
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::ForwardTransformPoint(const double in[3], double out[3])
{
  double displacement[3] = { 0.0, 0.0, 0.0 };
  this->GetDisplacement(in, displacement);
  out[0] = in[0] + displacement[0];
  out[1] = in[1] + displacement[1];
  out[2] = in[2] + displacement[2];
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::ForwardTransformPoint(const float in[3], float out[3])
{
  double inDouble[3] = { in[0], in[1], in[2] };
  double outDouble[3] = { 0.0, 0.0, 0.0 };
  this->ForwardTransformPoint(inDouble, outDouble);
  out[0] = static_cast<float>(outDouble[0]);
  out[1] = static_cast<float>(outDouble[1]);
  out[2] = static_cast<float>(outDouble[2]);
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::ForwardTransformDerivative(const double in[3], double out[3], double derivative[3][3])
{
  this->ForwardTransformPoint(in, out);
  // Central differences of the displacement; the field is only piecewise smooth, so this is
  // what the iterative inverse can use everywhere, also across cell boundaries
  double step = this->DerivativeStep;
  for (int column = 0; column < 3; ++column)
  {
    double before[3] = { in[0], in[1], in[2] };
    double after[3] = { in[0], in[1], in[2] };
    before[column] -= step;
    after[column] += step;
    double displacementBefore[3] = { 0.0, 0.0, 0.0 };
    double displacementAfter[3] = { 0.0, 0.0, 0.0 };
    this->GetDisplacement(before, displacementBefore);
    this->GetDisplacement(after, displacementAfter);
    for (int row = 0; row < 3; ++row)
    {
      derivative[row][column] = (row == column ? 1.0 : 0.0) + (displacementAfter[row] - displacementBefore[row]) / (2.0 * step);
    }
  }
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::ForwardTransformDerivative(const float in[3], float out[3], float derivative[3][3])
{
  double inDouble[3] = { in[0], in[1], in[2] };
  double outDouble[3] = { 0.0, 0.0, 0.0 };
  double derivativeDouble[3][3];
  this->ForwardTransformDerivative(inDouble, outDouble, derivativeDouble);
  for (int row = 0; row < 3; ++row)
  {
    out[row] = static_cast<float>(outDouble[row]);
    for (int column = 0; column < 3; ++column)
    {
      derivative[row][column] = static_cast<float>(derivativeDouble[row][column]);
    }
  }
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::ReportInverseConvergenceFailure(const double position[3], double error, int numberOfIterations)
{
  ++this->NumberOfInverseConvergenceFailures;
  // As vtkOrientedGridTransform: the first failure after the transform is modified is logged, unless all are requested
  bool log = this->LogAllInverseConvergenceFailures;
  if (!log)
  {
    const vtkMTimeType mTime = this->GetMTime();
    vtkMTimeType lastWarningTime = this->LastInverseConvergenceWarningTime.load();
    log = (mTime > lastWarningTime && this->LastInverseConvergenceWarningTime.compare_exchange_strong(lastWarningTime, mTime));
  }
  if (!log)
  {
    return;
  }
  vtkWarningMacro("InverseTransformPoint: no convergence (" << position[0] << ", " << position[1] << ", " << position[2] << ") error = " << error << " after " << numberOfIterations
                                                            << " iterations."
                                                            << (this->LogAllInverseConvergenceFailures ? ""
                                                                                                       : "  Further convergence warnings suppressed until the transform is "
                                                                                                         "modified (see GetNumberOfInverseConvergenceFailures())."));
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::InverseTransformPoint(const double in[3], double out[3])
{
  const double tolerance2 = this->InverseTolerance * this->InverseTolerance;
  // Squared distance between where a position is moved and the target, and the difference
  auto getResidual = [this, in](const double position[3], double residual[3])
  {
    this->ForwardTransformPoint(position, residual);
    for (int axis = 0; axis < 3; ++axis)
    {
      residual[axis] -= in[axis];
    }
    return vtkMath::Dot(residual, residual);
  };

  // First estimate: moved back by the displacement at the target
  double displacement[3] = { 0.0, 0.0, 0.0 };
  this->GetDisplacement(in, displacement);
  double position[3] = { in[0] - displacement[0], in[1] - displacement[1], in[2] - displacement[2] };
  double residual[3] = { 0.0, 0.0, 0.0 };
  double error2 = getResidual(position, residual);
  double bestPosition[3] = { position[0], position[1], position[2] };
  double bestError2 = error2;

  // Newton's method. The Jacobian costs 6 evaluations, so it is reused while the error decreases fast; a step that would
  // increase the error is shortened. A fresh Jacobian that cannot decrease the error means that the iteration is stuck
  // (where the transform folds, or at a discontinuity), so it stops there instead of running all the iterations.
  // Converged when the error is below the tolerance both where the position is moved and at the position itself (the
  // Newton step), as in vtkWarpTransform.
  double jacobian[3][3];
  bool hasJacobian = false;
  bool isJacobianFresh = false;
  bool converged = false;
  int iteration = 0;
  for (; iteration < this->InverseIterations; ++iteration)
  {
    if (!hasJacobian)
    {
      double moved[3] = { 0.0, 0.0, 0.0 };
      this->ForwardTransformDerivative(position, moved, jacobian);
      hasJacobian = true;
      isJacobianFresh = true;
    }
    double step[3] = { 0.0, 0.0, 0.0 };
    vtkMath::LinearSolve3x3(jacobian, residual, step);
    const bool isStepValid = std::isfinite(step[0]) && std::isfinite(step[1]) && std::isfinite(step[2]);
    if (isStepValid && error2 < tolerance2 && vtkMath::Dot(step, step) < tolerance2)
    {
      converged = true;
      break;
    }
    bool decreased = false;
    double nextPosition[3] = { 0.0, 0.0, 0.0 };
    double nextResidual[3] = { 0.0, 0.0, 0.0 };
    double nextError2 = error2;
    double stepFraction = 1.0;
    for (int halving = 0; isStepValid && halving < 6; ++halving, stepFraction *= 0.5)
    {
      for (int axis = 0; axis < 3; ++axis)
      {
        nextPosition[axis] = position[axis] - stepFraction * step[axis];
      }
      nextError2 = getResidual(nextPosition, nextResidual);
      if (nextError2 < error2)
      {
        decreased = true;
        break;
      }
    }
    if (!decreased)
    {
      if (isJacobianFresh)
      {
        break;
      }
      hasJacobian = false;
      continue;
    }
    // With a Jacobian that is no longer accurate enough for fast convergence, it is computed again for the next step
    if (nextError2 > 0.25 * error2)
    {
      hasJacobian = false;
    }
    isJacobianFresh = false;
    for (int axis = 0; axis < 3; ++axis)
    {
      position[axis] = nextPosition[axis];
      residual[axis] = nextResidual[axis];
    }
    error2 = nextError2;
    if (error2 < bestError2)
    {
      bestError2 = error2;
      std::copy(position, position + 3, bestPosition);
    }
  }

  if (!converged)
  {
    this->ReportInverseConvergenceFailure(in, std::sqrt(bestError2), iteration);
  }
  std::copy(bestPosition, bestPosition + 3, out);
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::InverseTransformPoint(const float in[3], float out[3])
{
  double inDouble[3] = { in[0], in[1], in[2] };
  double outDouble[3] = { 0.0, 0.0, 0.0 };
  this->InverseTransformPoint(inDouble, outDouble);
  out[0] = static_cast<float>(outDouble[0]);
  out[1] = static_cast<float>(outDouble[1]);
  out[2] = static_cast<float>(outDouble[2]);
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::InverseTransformDerivative(const double in[3], double out[3], double derivative[3][3])
{
  // The derivative of the forward transform at the inverse position; vtkWarpTransform inverts it
  this->InverseTransformPoint(in, out);
  double moved[3] = { 0.0, 0.0, 0.0 };
  this->ForwardTransformDerivative(out, moved, derivative);
}

//----------------------------------------------------------------------------
void vtkMeshDisplacementTransform::InverseTransformDerivative(const float in[3], float out[3], float derivative[3][3])
{
  double inDouble[3] = { in[0], in[1], in[2] };
  double outDouble[3] = { 0.0, 0.0, 0.0 };
  double derivativeDouble[3][3];
  this->InverseTransformDerivative(inDouble, outDouble, derivativeDouble);
  for (int row = 0; row < 3; ++row)
  {
    out[row] = static_cast<float>(outDouble[row]);
    for (int column = 0; column < 3; ++column)
    {
      derivative[row][column] = static_cast<float>(derivativeDouble[row][column]);
    }
  }
}
