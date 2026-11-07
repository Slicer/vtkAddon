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

// Measures how well vtkMeshDisplacementTransform can be inverted, and how
// smooth it is, around the border of convex and non-convex meshes, for each
// border treatment, and checks that the default of the class is the best one.
// The results that decided the default are described in
// vtkMeshDisplacementTransform::GetDisplacement().

// vtkAddon includes
#include <vtkAddonTestingMacros.h>
#include <vtkMeshDisplacementTransform.h>
#include <vtkTestingOutputWindow.h>

// VTK includes
#include <vtkBSplineTransform.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkImageData.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>

// STD includes
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{

// Half-size of the cube that the meshes fill
const double MeshHalfSize = 50.0;
// Round trip errors above this count as a failed inversion
const double FailureTolerance = 0.01;
// The inverse stops within its tolerance (0.001 mm, vtkWarpTransform's default) of the position and of
// the target, so a converged round trip is within twice that
const double MaximumConvergedRoundTripError = 2e-3;

//----------------------------------------------------------------------------
// Smooth, non-linear reference transform, as in vtkMeshDisplacementTransformTest2
void CreateReferenceTransform(vtkBSplineTransform* transform)
{
  vtkNew<vtkImageData> coefficients;
  coefficients->SetOrigin(-80.0, -80.0, -80.0);
  coefficients->SetSpacing(20.0, 20.0, 20.0);
  coefficients->SetDimensions(9, 9, 9);
  coefficients->AllocateScalars(VTK_DOUBLE, 3);
  for (int k = 0; k < 9; k++)
    {
    for (int j = 0; j < 9; j++)
      {
      for (int i = 0; i < 9; i++)
        {
        const double x = -80.0+i*20.0;
        const double y = -80.0+j*20.0;
        const double z = -80.0+k*20.0;
        double* coefficient = static_cast<double*>(coefficients->GetScalarPointer(i, j, k));
        coefficient[0] = 6.0*sin(x/30.0)*cos(z/45.0);
        coefficient[1] = 5.0*cos(y/25.0)*sin(z/35.0+0.3);
        coefficient[2] = 4.0*sin((x+y)/40.0)+0.002*z*z;
        }
      }
    }
  transform->SetCoefficientData(coefficients);
  transform->Update();
}

//----------------------------------------------------------------------------
// Grid of hexahedra, or of tetrahedra (6 per hexahedron, sharing its diagonal),
// that fills the cube, or the cube without the x > 0, y > 0 quadrant (an
// L-shaped prism). Only the points of the cells are created.
void CreateGrid(vtkUnstructuredGrid* grid, int resolution, bool lShaped, bool tetrahedra = false)
{
  const double step = 2.0*MeshHalfSize/(resolution-1);
  vtkNew<vtkPoints> points;
  std::map<int, vtkIdType> pointIds;
  auto getPointId = [&](int i, int j, int k) -> vtkIdType
    {
    const int key = (k*resolution+j)*resolution+i;
    auto it = pointIds.find(key);
    if (it != pointIds.end())
      {
      return it->second;
      }
    vtkIdType id = points->InsertNextPoint(-MeshHalfSize+i*step, -MeshHalfSize+j*step, -MeshHalfSize+k*step);
    pointIds[key] = id;
    return id;
    };

  grid->Allocate((resolution-1)*(resolution-1)*(resolution-1)*(tetrahedra ? 6 : 1));
  vtkNew<vtkHexahedron> hexahedron;
  vtkNew<vtkTetra> tetrahedron;
  // Corners as bit masks (x, y, z), in the VTK hexahedron order, and the 6 tetrahedra
  const int hexahedronCorners[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
  const int tetrahedronCorners[6][4] = { { 0, 1, 3, 7 }, { 0, 1, 5, 7 }, { 0, 2, 3, 7 }, { 0, 2, 6, 7 }, { 0, 4, 5, 7 }, { 0, 4, 6, 7 } };
  const int half = (resolution-1)/2;
  for (int k = 0; k < resolution-1; k++)
    {
    for (int j = 0; j < resolution-1; j++)
      {
      for (int i = 0; i < resolution-1; i++)
        {
        if (lShaped && i >= half && j >= half)
          {
          continue;
          }
        vtkIdType cornerIds[8];
        for (int corner = 0; corner < 8; corner++)
          {
          cornerIds[corner] = getPointId(i+(corner&1), j+((corner>>1)&1), k+((corner>>2)&1));
          }
        if (tetrahedra)
          {
          for (int tetrahedronIndex = 0; tetrahedronIndex < 6; tetrahedronIndex++)
            {
            for (int cornerIndex = 0; cornerIndex < 4; cornerIndex++)
              {
              tetrahedron->GetPointIds()->SetId(cornerIndex, cornerIds[tetrahedronCorners[tetrahedronIndex][cornerIndex]]);
              }
            grid->InsertNextCell(tetrahedron->GetCellType(), tetrahedron->GetPointIds());
            }
          }
        else
          {
          for (int cornerIndex = 0; cornerIndex < 8; cornerIndex++)
            {
            hexahedron->GetPointIds()->SetId(cornerIndex, cornerIds[hexahedronCorners[cornerIndex]]);
            }
          grid->InsertNextCell(hexahedron->GetCellType(), hexahedron->GetPointIds());
          }
        }
      }
    }
  grid->SetPoints(points);
}

//----------------------------------------------------------------------------
// Add the displacement of the reference transform at each mesh point; returns
// the largest displacement
double AddDisplacementArray(vtkPointSet* mesh, vtkAbstractTransform* referenceTransform)
{
  vtkPoints* points = mesh->GetPoints();
  vtkNew<vtkDoubleArray> displacementArray;
  displacementArray->SetName("Displacement");
  displacementArray->SetNumberOfComponents(3);
  displacementArray->SetNumberOfTuples(points->GetNumberOfPoints());
  double largestDisplacement = 0.0;
  for (vtkIdType pointIndex = 0; pointIndex < points->GetNumberOfPoints(); pointIndex++)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    points->GetPoint(pointIndex, point);
    double transformedPoint[3] = { 0.0, 0.0, 0.0 };
    referenceTransform->TransformPoint(point, transformedPoint);
    double displacement[3] = { transformedPoint[0]-point[0], transformedPoint[1]-point[1], transformedPoint[2]-point[2] };
    displacementArray->SetTuple(pointIndex, displacement);
    largestDisplacement = std::max(largestDisplacement, vtkMath::Norm(displacement));
    }
  mesh->GetPointData()->AddArray(displacementArray);
  return largestDisplacement;
}

//----------------------------------------------------------------------------
// Signed distance of a position from the surface of the cube (negative inside)
double DistanceFromCube(const double position[3])
{
  double outside[3] = { 0.0, 0.0, 0.0 };
  double inside = -VTK_DOUBLE_MAX;
  for (int axis = 0; axis < 3; axis++)
    {
    const double d = std::abs(position[axis])-MeshHalfSize;
    outside[axis] = std::max(0.0, d);
    inside = std::max(inside, d);
    }
  const double outsideDistance = vtkMath::Norm(outside);
  return (outsideDistance > 0.0 ? outsideDistance : inside);
}

//----------------------------------------------------------------------------
// Positions along rays that leave the cube through its faces, edges and
// corners, at every depth from a few mm inside to the given distance outside,
// so that each falloff distance is sampled at every depth. None of them is on
// a face or at a point of the mesh.
std::vector<std::vector<double>> CreateBandPositions(double bandOutside)
{
  // Start point and direction of each ray. Two of the offsets are just past a cell boundary of the mesh
  // (12.5 mm cells), so that the derivative, by central differences over a shorter step, straddles it:
  // this catches a displacement that jumps where the closest cell changes.
  std::vector<std::pair<std::vector<double>, std::vector<double>>> rays;
  const double offsets[6] = { -40.0+0.37, -25.0+0.005, 0.37, 12.5+0.005, 20.0+0.37, 40.0+0.37 };
  for (int axis = 0; axis < 3; axis++)
    {
    const int axisU = (axis+1)%3;
    const int axisV = (axis+2)%3;
    for (int sign = -1; sign <= 1; sign += 2)
      {
      // Faces
      for (double u : offsets)
        {
        for (double v : offsets)
          {
          std::vector<double> start(3, 0.0), direction(3, 0.0);
          start[axis] = sign*MeshHalfSize;
          start[axisU] = u;
          start[axisV] = v+0.34;
          direction[axis] = sign;
          rays.push_back({ start, direction });
          }
        }
      // Edges along this axis
      for (int signV = -1; signV <= 1; signV += 2)
        {
        for (double w : offsets)
          {
          std::vector<double> start(3, 0.0), direction(3, 0.0);
          start[axis] = w;
          start[axisU] = sign*MeshHalfSize;
          start[axisV] = signV*MeshHalfSize;
          direction[axisU] = sign/std::sqrt(2.0);
          direction[axisV] = signV/std::sqrt(2.0);
          rays.push_back({ start, direction });
          }
        }
      }
    }
  // Corners
  for (int signX = -1; signX <= 1; signX += 2)
    {
    for (int signY = -1; signY <= 1; signY += 2)
      {
      for (int signZ = -1; signZ <= 1; signZ += 2)
        {
        rays.push_back({ { signX*MeshHalfSize, signY*MeshHalfSize, signZ*MeshHalfSize },
                         { signX/std::sqrt(3.0), signY/std::sqrt(3.0), signZ/std::sqrt(3.0) } });
        }
      }
    }

  std::vector<std::vector<double>> positions;
  const double depthStep = 5.0;
  for (const auto& ray : rays)
    {
    for (double depth = -6.0; depth <= bandOutside; depth += depthStep)
      {
      positions.push_back({ ray.first[0]+depth*ray.second[0], ray.first[1]+depth*ray.second[1], ray.first[2]+depth*ray.second[2] });
      }
    }
  return positions;
}

//----------------------------------------------------------------------------
struct Metrics
{
  int NumberOfFailures = 0;
  double MaximumRoundTripError = 0.0;
  double RMSRoundTripError = 0.0;
  double MaximumGradient = 0.0;
  double RMSGradient = 0.0;
  double MaximumGradientChange = 0.0;
  std::string GradientChangeLocation; // where the largest change of the gradient is
  double Seconds = 0.0;
};

//----------------------------------------------------------------------------
double FrobeniusNorm(const double matrix[3][3])
{
  double sum = 0.0;
  for (int row = 0; row < 3; row++)
    {
    for (int column = 0; column < 3; column++)
      {
      sum += matrix[row][column]*matrix[row][column];
      }
    }
  return std::sqrt(sum);
}

//----------------------------------------------------------------------------
// Size of the deformation at a position: ||J - I||, the gradient of the displacement
double GetGradient(vtkMeshDisplacementTransform* transform, const double position[3])
{
  double transformed[3] = { 0.0, 0.0, 0.0 };
  double jacobian[3][3];
  transform->InternalTransformDerivative(position, transformed, jacobian);
  for (int axis = 0; axis < 3; axis++)
    {
    jacobian[axis][axis] -= 1.0;
    }
  return FrobeniusNorm(jacobian);
}

//----------------------------------------------------------------------------
// Largest change of the gradient per mm along rays that leave the cube through
// its faces, from a few mm inside to the given distance outside. The slope of
// a falloff profile jumps where it starts and where it ends, on surfaces that
// such a ray crosses at a known place, which positions on a grid would miss.
double GetMaximumGradientChangeAlongRays(vtkMeshDisplacementTransform* transform, double rayLength, std::string& location)
{
  const double step = 0.5;
  double maximumChange = 0.0;
  location.clear();
  for (int axis = 0; axis < 3; axis++)
    {
    for (int sign = -1; sign <= 1; sign += 2)
      {
      const double offsets[6] = { -40.0+0.37, -25.0+0.005, 0.37, 12.5+0.005, 20.0+0.37, 40.0+0.37 };
      for (double u : offsets)
        {
        for (double v : offsets)
          {
          double previousGradient = 0.0;
          for (double s = -3.0; s <= rayLength; s += step)
            {
            double position[3] = { 0.0, 0.0, 0.0 };
            position[axis] = sign*(MeshHalfSize+s);
            position[(axis+1)%3] = u;
            position[(axis+2)%3] = v;
            const double gradient = GetGradient(transform, position);
            if (s > -3.0 && std::abs(gradient-previousGradient)/step > maximumChange)
              {
              maximumChange = std::abs(gradient-previousGradient)/step;
              std::ostringstream where;
              where << std::fixed << std::setprecision(1) << "(" << position[0] << ", " << position[1] << ", " << position[2]
                    << "), " << s << " mm along " << (sign > 0 ? "+" : "-") << "xyz"[axis] << ", gradient "
                    << std::setprecision(2) << previousGradient << " -> " << gradient;
              location = where.str();
              }
            previousGradient = gradient;
            }
          }
        }
      }
    }
  return maximumChange;
}

//----------------------------------------------------------------------------
// Inverse consistency and smoothness of the transform at the positions
Metrics Evaluate(vtkMeshDisplacementTransform* transform, const std::vector<std::vector<double>>& positions, double rayLength)
{
  Metrics metrics;
  auto startTime = std::chrono::steady_clock::now();
  vtkAbstractTransform* inverse = transform->GetInverse();
  double sumRoundTripError2 = 0.0;
  double sumGradient2 = 0.0;
  for (size_t positionIndex = 0; positionIndex < positions.size(); positionIndex++)
    {
    const double* position = positions[positionIndex].data();

    // Round trip
    double transformed[3] = { 0.0, 0.0, 0.0 };
    double roundTripped[3] = { 0.0, 0.0, 0.0 };
    transform->TransformPoint(position, transformed);
    inverse->TransformPoint(transformed, roundTripped);
    const double roundTripError = std::sqrt(vtkMath::Distance2BetweenPoints(position, roundTripped));
    metrics.MaximumRoundTripError = std::max(metrics.MaximumRoundTripError, roundTripError);
    sumRoundTripError2 += roundTripError*roundTripError;
    if (roundTripError > FailureTolerance)
      {
      metrics.NumberOfFailures++;
      }

    // Gradient of the displacement
    const double gradient = GetGradient(transform, position);
    metrics.MaximumGradient = std::max(metrics.MaximumGradient, gradient);
    sumGradient2 += gradient*gradient;
    }
  metrics.RMSRoundTripError = std::sqrt(sumRoundTripError2/positions.size());
  metrics.RMSGradient = std::sqrt(sumGradient2/positions.size());
  metrics.MaximumGradientChange = GetMaximumGradientChangeAlongRays(transform, rayLength, metrics.GradientChangeLocation);
  metrics.Seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-startTime).count();
  return metrics;
}

//----------------------------------------------------------------------------
struct Strategy
{
  std::string Name;
  double DistanceFactor; // multiple of the largest displacement, 0 for no falloff, negative for constant extension
};

//----------------------------------------------------------------------------
void PrintHeader()
{
  std::cout << std::left << std::setw(16) << "strategy" << std::right
            << std::setw(9) << "failures" << std::setw(11) << "maxRT[mm]" << std::setw(11) << "rmsRT[mm]"
            << std::setw(9) << "maxGrad" << std::setw(9) << "rmsGrad" << std::setw(11) << "maxGradChg"
            << std::setw(8) << "sec" << std::endl;
}

//----------------------------------------------------------------------------
void PrintMetrics(const Strategy& strategy, const Metrics& metrics)
{
  std::cout << std::left << std::setw(16) << strategy.Name << std::right << std::fixed
            << std::setw(9) << metrics.NumberOfFailures
            << std::setw(11) << std::setprecision(4) << metrics.MaximumRoundTripError
            << std::setw(11) << std::setprecision(5) << metrics.RMSRoundTripError
            << std::setw(9) << std::setprecision(3) << metrics.MaximumGradient
            << std::setw(9) << std::setprecision(3) << metrics.RMSGradient
            << std::setw(11) << std::setprecision(3) << metrics.MaximumGradientChange
            << std::setw(8) << std::setprecision(1) << metrics.Seconds << std::endl;
  if (metrics.MaximumGradientChange > 1.0)
    {
    std::cout << "    largest gradient change at " << metrics.GradientChangeLocation << std::endl;
    }
}

//----------------------------------------------------------------------------
std::vector<Strategy> CreateStrategies()
{
  std::vector<Strategy> strategies;
  strategies.push_back({ "none", 0.0 });
  // The falloff at multiples of the largest displacement around its theoretical limit of
  // invertibility (1.5, the largest slope of the cubic weight), and larger
  for (double factor : { 1.0, 1.5, 2.0, 3.0 })
    {
    strategies.push_back({ "falloff x" + std::to_string(factor).substr(0, 3), factor });
    }
  strategies.push_back({ "constant", -1.0 });
  return strategies;
}

//----------------------------------------------------------------------------
void ApplyStrategy(vtkMeshDisplacementTransform* transform, const Strategy& strategy, double largestDisplacement)
{
  transform->AutomaticBorderFalloffDistanceOff();
  if (strategy.DistanceFactor < 0.0)
    {
    transform->SetBorderFalloffDistance(VTK_DOUBLE_MAX);
    }
  else
    {
    transform->SetBorderFalloffDistance(strategy.DistanceFactor*largestDisplacement);
    }
}

//----------------------------------------------------------------------------
// Run every strategy on a mesh; returns the metrics in the order of the strategies
std::vector<Metrics> EvaluateMesh(const char* meshName, bool lShaped, bool tetrahedra, vtkBSplineTransform* referenceTransform,
                                  const std::vector<Strategy>& strategies)
{
  vtkNew<vtkUnstructuredGrid> mesh;
  CreateGrid(mesh, 9, lShaped, tetrahedra);
  const double largestDisplacement = AddDisplacementArray(mesh, referenceTransform);
  const double bandOutside = 3.0*largestDisplacement+5.0;
  std::vector<std::vector<double>> positions = CreateBandPositions(bandOutside);
  std::cout << std::endl << meshName << ": " << mesh->GetNumberOfCells() << " cells, largest displacement "
            << std::setprecision(3) << largestDisplacement << " mm, " << positions.size() << " positions" << std::endl;
  PrintHeader();

  vtkNew<vtkMeshDisplacementTransform> transform;
  transform->SetInputPointSet(mesh);
  transform->SetDisplacementArrayName("Displacement");
  std::vector<Metrics> results;
  for (const Strategy& strategy : strategies)
    {
    ApplyStrategy(transform, strategy, largestDisplacement);
    results.push_back(Evaluate(transform, positions, bandOutside));
    PrintMetrics(strategy, results.back());
    }
  return results;
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkMeshDisplacementTransformTest3(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  TESTING_OUTPUT_INIT();

  vtkNew<vtkBSplineTransform> referenceTransform;
  CreateReferenceTransform(referenceTransform);
  std::vector<Strategy> strategies = CreateStrategies();

  // Strategies that cannot be inverted make vtkWarpTransform warn about the iterations
  TESTING_OUTPUT_IGNORE_WARNINGS_ERRORS_BEGIN();
  std::vector<Metrics> cubeResults = EvaluateMesh("Cube of hexahedra", false, false, referenceTransform, strategies);
  std::vector<Metrics> lResults = EvaluateMesh("L-shaped prism of hexahedra", true, false, referenceTransform, strategies);
  // The default strategy and the constant extension on tetrahedra, whose gradient jumps between cells
  std::vector<Strategy> tetrahedraStrategies;
  for (const Strategy& strategy : strategies)
    {
    if (strategy.DistanceFactor < 0.0 || strategy.DistanceFactor >= 2.0)
      {
      tetrahedraStrategies.push_back(strategy);
      }
    }
  std::vector<Metrics> tetrahedraResults = EvaluateMesh("Cube of tetrahedra", false, true, referenceTransform, tetrahedraStrategies);
  TESTING_OUTPUT_IGNORE_WARNINGS_ERRORS_END();

  // The default of the class: the falloff at twice the largest displacement (the automatic distance)
  int defaultIndex = -1;
  for (size_t index = 0; index < strategies.size(); index++)
    {
    if (strategies[index].DistanceFactor == 2.0)
      {
      defaultIndex = static_cast<int>(index);
      }
    }
  CHECK_BOOL(defaultIndex >= 0, true);
  std::cout << std::endl << "Default strategy: " << strategies[defaultIndex].Name << std::endl;

  // Without a falloff the inverse fails near the border; with the default it does not, on either mesh
  CHECK_BOOL(cubeResults[0].NumberOfFailures > 0, true);
  CHECK_INT(cubeResults[defaultIndex].NumberOfFailures, 0);
  CHECK_INT(lResults[defaultIndex].NumberOfFailures, 0);
  CHECK_BOOL(cubeResults[defaultIndex].MaximumRoundTripError < MaximumConvergedRoundTripError, true);
  CHECK_BOOL(lResults[defaultIndex].MaximumRoundTripError < MaximumConvergedRoundTripError, true);

  // The constant extension is the smoothest on a convex mesh, but the closest mesh point switches
  // abruptly at the concave corner of the L, which folds the transform there: a finite falloff
  // is the robust choice
  CHECK_INT(cubeResults.back().NumberOfFailures, 0);
  CHECK_BOOL(lResults.back().NumberOfFailures > 0, true);

  // The default strategy inverts on tetrahedra as well, as smoothly as on hexahedra
  CHECK_INT(tetrahedraResults[0].NumberOfFailures, 0);
  CHECK_BOOL(tetrahedraResults[0].MaximumRoundTripError < MaximumConvergedRoundTripError, true);
  CHECK_BOOL(tetrahedraResults[0].MaximumGradientChange < 1.2*cubeResults[defaultIndex].MaximumGradientChange, true);

  // Below the theoretical limit of the cubic weight (1.5 times the largest displacement) the inverse
  // fails; a larger distance makes the falloff smoother (its gradient scales with 1 / distance)
  CHECK_BOOL(cubeResults[1].NumberOfFailures > 0, true);
  CHECK_INT(cubeResults[2].NumberOfFailures, 0);
  CHECK_BOOL(cubeResults[4].MaximumGradient < cubeResults[defaultIndex].MaximumGradient, true);
  CHECK_BOOL(cubeResults[defaultIndex].MaximumGradient < cubeResults[2].MaximumGradient, true);

  std::cout << "vtkMeshDisplacementTransformTest3 completed successfully" << std::endl;
  return EXIT_SUCCESS;
}
