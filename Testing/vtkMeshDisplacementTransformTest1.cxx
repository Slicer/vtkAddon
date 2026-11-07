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

// Checks vtkMeshDisplacementTransform on simple meshes with a linear
// displacement field, which the transform must reproduce exactly. Covers the
// transform, its derivative and inverse inside the mesh, each border treatment
// outside it, surface meshes and point clouds, and invalid inputs.

// vtkAddon includes
#include <vtkAddonTestingMacros.h>
#include <vtkMeshDisplacementTransform.h>
#include <vtkTestingOutputWindow.h>

// VTK includes
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>

// STD includes
#include <chrono>
#include <iostream>

namespace
{

//----------------------------------------------------------------------------
// Linear displacement field: interpolation in a hexahedron reproduces it
// exactly, so it is the expected displacement everywhere inside the mesh.
void GetReferenceDisplacement(const double point[3], double displacement[3])
{
  displacement[0] = 2.0 + 0.05*point[1];
  displacement[1] = -1.0 + 0.03*point[2];
  displacement[2] = 0.5 + 0.04*point[0];
}

//----------------------------------------------------------------------------
// Jacobian of the transform that the reference displacement defines
// (identity + gradient of the displacement), row by row.
const double ReferenceDerivative[3][3] = {
  { 1.0, 0.05, 0.0 },
  { 0.0, 1.0, 0.03 },
  { 0.04, 0.0, 1.0 } };

//----------------------------------------------------------------------------
// Build a volumetric mesh: a regular grid of hexahedral cells that fills the
// cube of the given half-size.
void CreateHexahedralGrid(vtkUnstructuredGrid* grid, int resolution, double halfSize)
{
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(resolution*resolution*resolution);
  const double step = 2.0*halfSize/(resolution-1);
  for (int k = 0; k < resolution; k++)
    {
    for (int j = 0; j < resolution; j++)
      {
      for (int i = 0; i < resolution; i++)
        {
        points->SetPoint((k*resolution+j)*resolution+i,
          -halfSize+i*step, -halfSize+j*step, -halfSize+k*step);
        }
      }
    }
  grid->SetPoints(points);

  grid->Allocate((resolution-1)*(resolution-1)*(resolution-1));
  vtkNew<vtkHexahedron> hexahedron;
  for (int k = 0; k < resolution-1; k++)
    {
    for (int j = 0; j < resolution-1; j++)
      {
      for (int i = 0; i < resolution-1; i++)
        {
        const vtkIdType corner = (k*resolution+j)*resolution+i;
        const vtkIdType cornerIds[8] = {
          corner, corner+1, corner+resolution+1, corner+resolution,
          corner+resolution*resolution, corner+resolution*resolution+1,
          corner+resolution*resolution+resolution+1, corner+resolution*resolution+resolution };
        for (int cornerIndex = 0; cornerIndex < 8; cornerIndex++)
          {
          hexahedron->GetPointIds()->SetId(cornerIndex, cornerIds[cornerIndex]);
          }
        grid->InsertNextCell(hexahedron->GetCellType(), hexahedron->GetPointIds());
        }
      }
    }
}

//----------------------------------------------------------------------------
// Add a 3-component point data array that contains the reference displacement
// at each point of the mesh.
void AddDisplacementArray(vtkPointSet* mesh, const char* arrayName, bool setAsActiveVectors = false)
{
  vtkPoints* points = mesh->GetPoints();
  vtkNew<vtkDoubleArray> displacementArray;
  displacementArray->SetName(arrayName);
  displacementArray->SetNumberOfComponents(3);
  displacementArray->SetNumberOfTuples(points->GetNumberOfPoints());
  for (vtkIdType pointIndex = 0; pointIndex < points->GetNumberOfPoints(); pointIndex++)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    points->GetPoint(pointIndex, point);
    double displacement[3] = { 0.0, 0.0, 0.0 };
    GetReferenceDisplacement(point, displacement);
    displacementArray->SetTuple(pointIndex, displacement);
    }
  mesh->GetPointData()->AddArray(displacementArray);
  if (setAsActiveVectors)
    {
    mesh->GetPointData()->SetActiveVectors(arrayName);
    }
}

//----------------------------------------------------------------------------
// Largest distance between the transformed mesh points and the points that the
// displacement array prescribes, with the displacements multiplied by scale.
double GetMaximumErrorAtMeshPoints(vtkAbstractTransform* transform, vtkPointSet* mesh,
                                   const char* arrayName, double scale = 1.0)
{
  vtkPoints* points = mesh->GetPoints();
  vtkDataArray* displacementArray = mesh->GetPointData()->GetArray(arrayName);
  double maximumError = 0.0;
  for (vtkIdType pointIndex = 0; pointIndex < points->GetNumberOfPoints(); pointIndex++)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    points->GetPoint(pointIndex, point);
    double displacement[3] = { 0.0, 0.0, 0.0 };
    displacementArray->GetTuple(pointIndex, displacement);
    double expected[3] = { point[0]+scale*displacement[0], point[1]+scale*displacement[1], point[2]+scale*displacement[2] };
    double actual[3] = { 0.0, 0.0, 0.0 };
    transform->TransformPoint(point, actual);
    maximumError = std::max(maximumError, sqrt(vtkMath::Distance2BetweenPoints(expected, actual)));
    }
  return maximumError;
}

//----------------------------------------------------------------------------
// Distance between where the transform moves a point and where it is expected.
double GetErrorAtPoint(vtkAbstractTransform* transform, const double point[3], const double expected[3])
{
  double actual[3] = { 0.0, 0.0, 0.0 };
  transform->TransformPoint(point, actual);
  return sqrt(vtkMath::Distance2BetweenPoints(expected, actual));
}

//----------------------------------------------------------------------------
// Thin shell of hexahedra on a saddle surface, one cell thick, as the wall of a
// vessel or a leaflet is meshed: its cell faces are not planar, and finding the
// parametric coordinates of a position near such a cell can diverge.
void CreateThinSaddleShell(vtkUnstructuredGrid* grid)
{
  const int resolution = 25;
  const double size = 10.0;
  const double thickness = 0.05;
  vtkNew<vtkPoints> points;
  for (int layer = 0; layer < 2; layer++)
    {
    for (int j = 0; j < resolution; j++)
      {
      for (int i = 0; i < resolution; i++)
        {
        const double x = -size/2.0+size*i/(resolution-1);
        const double y = -size/2.0+size*j/(resolution-1);
        points->InsertNextPoint(x, y, 0.15*(x*x-y*y)+layer*thickness);
        }
      }
    }
  grid->SetPoints(points);
  vtkNew<vtkHexahedron> hexahedron;
  const vtkIdType layerOffset = resolution*resolution;
  for (int j = 0; j < resolution-1; j++)
    {
    for (int i = 0; i < resolution-1; i++)
      {
      const vtkIdType corner = j*resolution+i;
      const vtkIdType cornerIds[8] = { corner, corner+1, corner+resolution+1, corner+resolution,
        corner+layerOffset, corner+layerOffset+1, corner+layerOffset+resolution+1, corner+layerOffset+resolution };
      for (int cornerIndex = 0; cornerIndex < 8; cornerIndex++)
        {
        hexahedron->GetPointIds()->SetId(cornerIndex, cornerIds[cornerIndex]);
        }
      grid->InsertNextCell(hexahedron->GetCellType(), hexahedron->GetPointIds());
      }
    }
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkMeshDisplacementTransformTest1(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  TESTING_OUTPUT_INIT();

  // A 100 mm cube of 4x4x4 cells, each 25 mm wide
  vtkNew<vtkUnstructuredGrid> volumetricMesh;
  CreateHexahedralGrid(volumetricMesh, 5, 50.0);
  AddDisplacementArray(volumetricMesh, "Displacement");

  vtkNew<vtkMeshDisplacementTransform> transform;
  transform->SetInputPointSet(volumetricMesh);
  transform->SetDisplacementArrayName("Displacement");
  transform->Update();
  CHECK_POINTER(transform->GetInputPointSet(), volumetricMesh.GetPointer());
  CHECK_STRING(transform->GetDisplacementArrayName(), "Displacement");
  CHECK_DOUBLE(transform->GetDisplacementScale(), 1.0);

  //---------------------------------------------------------------------------
  // The mesh points move exactly by their own vectors, as vtkWarpVector moves them
  std::cout << "Checking the mesh points" << std::endl;
  CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(transform, volumetricMesh, "Displacement"), 0.0, 1e-9);

  //---------------------------------------------------------------------------
  // Between the mesh points the displacement is interpolated in the cell that
  // contains the position
  std::cout << "Checking interpolation inside cells" << std::endl;
  double maximumInterpolationError = 0.0;
  for (int i = 0; i < 10; i++)
    {
    // Positions along a diagonal of the cube, none of them a mesh point
    double point[3] = { -46.0+i*9.7, -41.0+i*8.3, 47.0-i*9.9 };
    double displacement[3] = { 0.0, 0.0, 0.0 };
    GetReferenceDisplacement(point, displacement);
    double expected[3] = { point[0]+displacement[0], point[1]+displacement[1], point[2]+displacement[2] };
    maximumInterpolationError = std::max(maximumInterpolationError, GetErrorAtPoint(transform, point, expected));

    double transformDisplacement[3] = { 0.0, 0.0, 0.0 };
    CHECK_BOOL(transform->GetDisplacement(point, transformDisplacement), true);
    }
  std::cout << "  maximum interpolation error: " << maximumInterpolationError << std::endl;
  CHECK_DOUBLE_TOLERANCE(maximumInterpolationError, 0.0, 1e-9);

  //---------------------------------------------------------------------------
  // Outside the cells nothing moves
  std::cout << "Checking positions outside the mesh" << std::endl;
  {
    double outsidePoint[3] = { 80.0, 0.0, 0.0 };
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(transform, outsidePoint, outsidePoint), 0.0, 1e-12);
    double displacement[3] = { 1.0, 1.0, 1.0 };
    CHECK_BOOL(transform->GetDisplacement(outsidePoint, displacement), false);
    CHECK_DOUBLE(displacement[0], 0.0);
    CHECK_DOUBLE(displacement[1], 0.0);
    CHECK_DOUBLE(displacement[2], 0.0);
  }

  //---------------------------------------------------------------------------
  // With a falloff distance the displacement outside the cells is that of the
  // closest mesh point, falling smoothly to zero at the falloff distance
  std::cout << "Checking the falloff outside the mesh" << std::endl;
  CHECK_DOUBLE(transform->GetBorderFalloffDistance(), 0.0);
  transform->SetBorderFalloffDistance(30.0);
  {
    // Position outside the mesh, closest point on the mesh, expected weight
    const double falloffPoints[4][7] = {
      { 65.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.5 },
      { 57.5, 10.0, -20.0, 50.0, 10.0, -20.0, 0.84375 },
      { -30.0, -74.0, 5.0, -30.0, -50.0, 5.0, 0.104 },
      { 60.0, 60.0, 60.0, 50.0, 50.0, 50.0, -1.0 } };
    for (const double* falloffPoint : falloffPoints)
      {
      const double* point = falloffPoint;
      const double* closestPoint = falloffPoint + 3;
      double weight = falloffPoint[6];
      if (weight < 0.0)
        {
        double t = sqrt(vtkMath::Distance2BetweenPoints(point, closestPoint)) / 30.0;
        weight = (1.0-t)*(1.0-t)*(1.0+2.0*t);
        }
      double closestDisplacement[3] = { 0.0, 0.0, 0.0 };
      GetReferenceDisplacement(closestPoint, closestDisplacement);
      double expected[3] = { point[0]+weight*closestDisplacement[0], point[1]+weight*closestDisplacement[1],
        point[2]+weight*closestDisplacement[2] };
      CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(transform, point, expected), 0.0, 1e-9);
      double displacement[3] = { 0.0, 0.0, 0.0 };
      CHECK_BOOL(transform->GetDisplacement(point, displacement), true);
      }

    // Continuous at the boundary of the mesh
    double insidePoint[3] = { 50.0-1e-6, 12.0, 7.0 };
    double outsidePoint[3] = { 50.0+1e-6, 12.0, 7.0 };
    double insideDisplacement[3] = { 0.0, 0.0, 0.0 };
    double outsideDisplacement[3] = { 0.0, 0.0, 0.0 };
    transform->GetDisplacement(insidePoint, insideDisplacement);
    transform->GetDisplacement(outsidePoint, outsideDisplacement);
    CHECK_DOUBLE_TOLERANCE(sqrt(vtkMath::Distance2BetweenPoints(insideDisplacement, outsideDisplacement)), 0.0, 1e-5);

    // Zero at and beyond the falloff distance
    double farPoints[2][3] = { { 80.0, 0.0, 0.0 }, { 300.0, -200.0, 0.0 } };
    for (double* farPoint : farPoints)
      {
      CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(transform, farPoint, farPoint), 0.0, 1e-12);
      double displacement[3] = { 1.0, 1.0, 1.0 };
      CHECK_BOOL(transform->GetDisplacement(farPoint, displacement), false);
      }
  }

  //---------------------------------------------------------------------------
  // With the largest falloff distance the displacement of the closest mesh
  // point is kept everywhere outside the mesh
  std::cout << "Checking constant displacement outside the mesh" << std::endl;
  transform->SetBorderFalloffDistance(VTK_DOUBLE_MAX);
  {
    // Position outside the mesh, closest point on the mesh
    const double constantPoints[4][6] = {
      { 80.0, 0.0, 0.0, 50.0, 0.0, 0.0 },
      { 300.0, -200.0, 10.0, 50.0, -50.0, 10.0 },
      { -1e5, 20.0, -30.0, -50.0, 20.0, -30.0 },
      { 1e6, 1e6, -1e6, 50.0, 50.0, -50.0 } };
    for (const double* constantPoint : constantPoints)
      {
      const double* point = constantPoint;
      double closestDisplacement[3] = { 0.0, 0.0, 0.0 };
      GetReferenceDisplacement(constantPoint + 3, closestDisplacement);
      double expected[3] = { point[0]+closestDisplacement[0], point[1]+closestDisplacement[1], point[2]+closestDisplacement[2] };
      CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(transform, point, expected), 0.0, 1e-6);

      // The inverse works far from the mesh too
      double roundTrippedPoint[3] = { 0.0, 0.0, 0.0 };
      transform->GetInverse()->TransformPoint(expected, roundTrippedPoint);
      CHECK_DOUBLE_TOLERANCE(sqrt(vtkMath::Distance2BetweenPoints(point, roundTrippedPoint)), 0.0, 1e-2);
      }

    // Finding the closest mesh point is fast also far from the mesh
    auto startTime = std::chrono::steady_clock::now();
    const int numberOfFarPoints = 10000;
    for (int i = 0; i < numberOfFarPoints; i++)
      {
      double point[3] = { 100.0+i*37.0, -2000.0+i*0.5, 1000.0-i*11.0 };
      double transformedPoint[3] = { 0.0, 0.0, 0.0 };
      transform->TransformPoint(point, transformedPoint);
      }
    double elapsedSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
    std::cout << "  " << numberOfFarPoints << " far positions transformed in " << elapsedSeconds << " s" << std::endl;
    CHECK_BOOL(elapsedSeconds < 5.0, true);
  }

  //---------------------------------------------------------------------------
  // The automatic falloff distance is twice the largest displacement, and it
  // follows the displacement scale
  std::cout << "Checking the automatic falloff distance" << std::endl;
  CHECK_BOOL(transform->GetAutomaticBorderFalloffDistance(), false);
  transform->SetBorderFalloffDistance(30.0);
  transform->AutomaticBorderFalloffDistanceOn();
  {
    double largestDisplacement = 0.0;
    vtkDataArray* displacementArray = volumetricMesh->GetPointData()->GetArray("Displacement");
    for (vtkIdType pointIndex = 0; pointIndex < displacementArray->GetNumberOfTuples(); pointIndex++)
      {
      double vector[3] = { 0.0, 0.0, 0.0 };
      displacementArray->GetTuple(pointIndex, vector);
      largestDisplacement = std::max(largestDisplacement, vtkMath::Norm(vector));
      }
    transform->Update();
    CHECK_DOUBLE_TOLERANCE(transform->GetEffectiveBorderFalloffDistance(), 2.0*largestDisplacement, 1e-9);
    // The explicit distance is kept, but not used
    CHECK_DOUBLE(transform->GetBorderFalloffDistance(), 30.0);
    // Halfway to the automatic distance the weight is 0.5
    double closestPoint[3] = { 50.0, 0.0, 0.0 };
    double point[3] = { 50.0+transform->GetEffectiveBorderFalloffDistance()/2.0, 0.0, 0.0 };
    double closestDisplacement[3] = { 0.0, 0.0, 0.0 };
    GetReferenceDisplacement(closestPoint, closestDisplacement);
    double expected[3] = { point[0]+0.5*closestDisplacement[0], point[1]+0.5*closestDisplacement[1], point[2]+0.5*closestDisplacement[2] };
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(transform, point, expected), 0.0, 1e-9);
    // Scaled displacements, in either direction, scale the distance
    transform->SetDisplacementScale(-0.5);
    transform->Update();
    CHECK_DOUBLE_TOLERANCE(transform->GetEffectiveBorderFalloffDistance(), largestDisplacement, 1e-9);
    transform->SetDisplacementScale(1.0);
    // Turned off, the explicit distance is used again
    transform->AutomaticBorderFalloffDistanceOff();
    transform->Update();
    CHECK_DOUBLE(transform->GetEffectiveBorderFalloffDistance(), 30.0);
  }

  //---------------------------------------------------------------------------
  // Near the boundary of the mesh the inverse only works with a falloff: without
  // it, points moved out of the mesh are not moved back, as nothing moves there
  std::cout << "Checking the inverse near the boundary of the mesh" << std::endl;
  for (double borderFalloffDistance : { 0.0, 30.0 })
    {
    transform->SetBorderFalloffDistance(borderFalloffDistance);
    double maximumBoundaryInverseError = 0.0;
    for (int i = 0; i < 12; i++)
      {
      for (double y : { -30.0, 0.0, 45.0 })
        {
        double point[3] = { 38.0+i*3.0, y, 0.5*y-20.0 };
        double transformedPoint[3] = { 0.0, 0.0, 0.0 };
        transform->TransformPoint(point, transformedPoint);
        double roundTrippedPoint[3] = { 0.0, 0.0, 0.0 };
        transform->GetInverse()->TransformPoint(transformedPoint, roundTrippedPoint);
        maximumBoundaryInverseError = std::max(maximumBoundaryInverseError,
          sqrt(vtkMath::Distance2BetweenPoints(point, roundTrippedPoint)));
        }
      }
    std::cout << "  falloff distance " << borderFalloffDistance << ": maximum inverse error " << maximumBoundaryInverseError << std::endl;
    CHECK_BOOL(maximumBoundaryInverseError > 1.0, borderFalloffDistance == 0.0);
    CHECK_BOOL(maximumBoundaryInverseError < 1e-2, borderFalloffDistance > 0.0);
    }
  transform->SetBorderFalloffDistance(0.0);

  //---------------------------------------------------------------------------
  // The derivative is that of the displacement field (computed by finite
  // differences, which are exact for a linear field)
  std::cout << "Checking the derivative" << std::endl;
  {
    double point[3] = { 3.0, -7.0, 11.0 };
    double transformedPoint[3] = { 0.0, 0.0, 0.0 };
    double derivative[3][3];
    transform->InternalTransformDerivative(point, transformedPoint, derivative);
    for (int row = 0; row < 3; row++)
      {
      for (int column = 0; column < 3; column++)
        {
        CHECK_DOUBLE_TOLERANCE(derivative[row][column], ReferenceDerivative[row][column], 1e-6);
        }
      }
    // The single-precision variant gives the same
    float pointFloat[3] = { 3.0f, -7.0f, 11.0f };
    float transformedPointFloat[3] = { 0.0f, 0.0f, 0.0f };
    float derivativeFloat[3][3];
    transform->InternalTransformDerivative(pointFloat, transformedPointFloat, derivativeFloat);
    for (int axis = 0; axis < 3; axis++)
      {
      CHECK_DOUBLE_TOLERANCE(transformedPointFloat[axis], transformedPoint[axis], 1e-4);
      CHECK_DOUBLE_TOLERANCE(derivativeFloat[axis][axis], ReferenceDerivative[axis][axis], 1e-4);
      }
  }

  //---------------------------------------------------------------------------
  // The inverse undoes the forward transform
  std::cout << "Checking the inverse transform" << std::endl;
  double maximumInverseError = 0.0;
  for (int i = 0; i < 10; i++)
    {
    double point[3] = { -30.0+i*6.1, 25.0-i*5.3, -20.0+i*4.7 };
    double transformedPoint[3] = { 0.0, 0.0, 0.0 };
    transform->TransformPoint(point, transformedPoint);
    double roundTrippedPoint[3] = { 0.0, 0.0, 0.0 };
    transform->GetInverse()->TransformPoint(transformedPoint, roundTrippedPoint);
    maximumInverseError = std::max(maximumInverseError, sqrt(vtkMath::Distance2BetweenPoints(point, roundTrippedPoint)));
    }
  std::cout << "  maximum inverse error: " << maximumInverseError << std::endl;
  CHECK_BOOL(maximumInverseError < 1e-3, true);

  // The inverse follows changes of the transform (the inverse is updated by
  // copying the transform, which must keep its direction)
  for (double scale : { 0.5, 1.0 })
    {
    transform->SetDisplacementScale(scale);
    double point[3] = { 12.0, -8.0, 21.0 };
    double transformedPoint[3] = { 0.0, 0.0, 0.0 };
    transform->TransformPoint(point, transformedPoint);
    double roundTrippedPoint[3] = { 0.0, 0.0, 0.0 };
    transform->GetInverse()->TransformPoint(transformedPoint, roundTrippedPoint);
    CHECK_DOUBLE_TOLERANCE(sqrt(vtkMath::Distance2BetweenPoints(point, roundTrippedPoint)), 0.0, 1e-3);
    }

  //---------------------------------------------------------------------------
  // The displacements can be scaled
  std::cout << "Checking the displacement scale" << std::endl;
  transform->SetDisplacementScale(0.5);
  CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(transform, volumetricMesh, "Displacement", 0.5), 0.0, 1e-9);
  transform->SetDisplacementScale(1.0);

  //---------------------------------------------------------------------------
  // On a surface mesh the points move by their vectors; a position off the
  // surface is not in any cell, so it does not move
  std::cout << "Checking a surface mesh" << std::endl;
  vtkNew<vtkSphereSource> sphereSource;
  sphereSource->SetRadius(50.0);
  sphereSource->SetThetaResolution(20);
  sphereSource->SetPhiResolution(20);
  sphereSource->Update();
  vtkNew<vtkPolyData> surfaceMesh;
  surfaceMesh->DeepCopy(sphereSource->GetOutput());
  AddDisplacementArray(surfaceMesh, "Displacement");

  vtkNew<vtkMeshDisplacementTransform> surfaceTransform;
  surfaceTransform->SetInputPointSet(surfaceMesh);
  surfaceTransform->SetDisplacementArrayName("Displacement");
  CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(surfaceTransform, surfaceMesh, "Displacement"), 0.0, 1e-9);
  {
    double centerPoint[3] = { 0.0, 0.0, 0.0 };
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(surfaceTransform, centerPoint, centerPoint), 0.0, 1e-12);
  }

  //---------------------------------------------------------------------------
  // A point cloud without cells: only its points move
  std::cout << "Checking a point cloud" << std::endl;
  vtkNew<vtkPolyData> pointCloud;
  pointCloud->SetPoints(volumetricMesh->GetPoints());
  AddDisplacementArray(pointCloud, "Displacement");
  vtkNew<vtkMeshDisplacementTransform> pointCloudTransform;
  pointCloudTransform->SetInputPointSet(pointCloud);
  pointCloudTransform->SetDisplacementArrayName("Displacement");
  CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(pointCloudTransform, pointCloud, "Displacement"), 0.0, 1e-9);
  {
    double betweenPoints[3] = { 10.0, 10.0, 10.0 };
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(pointCloudTransform, betweenPoints, betweenPoints), 0.0, 1e-12);
    // The falloff is not used without cells
    pointCloudTransform->SetBorderFalloffDistance(30.0);
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(pointCloudTransform, betweenPoints, betweenPoints), 0.0, 1e-12);
    CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(pointCloudTransform, pointCloud, "Displacement"), 0.0, 1e-9);
  }

  //---------------------------------------------------------------------------
  // If no displacement array name is specified then the active vectors are used
  std::cout << "Checking the active vectors array fallback" << std::endl;
  vtkNew<vtkUnstructuredGrid> activeVectorsMesh;
  CreateHexahedralGrid(activeVectorsMesh, 5, 50.0);
  AddDisplacementArray(activeVectorsMesh, "Displacement", true);
  vtkNew<vtkMeshDisplacementTransform> activeVectorsTransform;
  activeVectorsTransform->SetInputPointSet(activeVectorsMesh);
  CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(activeVectorsTransform, activeVectorsMesh, "Displacement"), 0.0, 1e-9);

  //---------------------------------------------------------------------------
  // Modifying the input mesh is noticed, without setting it again
  std::cout << "Checking that the transform follows changes of the input mesh" << std::endl;
  {
    double firstPoint[3] = { 0.0, 0.0, 0.0 };
    volumetricMesh->GetPoints()->GetPoint(0, firstPoint);
    vtkMTimeType transformMTime = transform->GetMTime();
    vtkDataArray* displacementArray = volumetricMesh->GetPointData()->GetArray("Displacement");
    double originalDisplacement[3] = { 0.0, 0.0, 0.0 };
    displacementArray->GetTuple(0, originalDisplacement);
    double modifiedDisplacement[3] = { 20.0, -20.0, 20.0 };
    displacementArray->SetTuple(0, modifiedDisplacement);
    displacementArray->Modified();
    volumetricMesh->Modified();
    CHECK_BOOL(transform->GetMTime() > transformMTime, true);
    double expected[3] = { firstPoint[0]+20.0, firstPoint[1]-20.0, firstPoint[2]+20.0 };
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(transform, firstPoint, expected), 0.0, 1e-9);
    // Restore the original displacement
    displacementArray->SetTuple(0, originalDisplacement);
    displacementArray->Modified();
    volumetricMesh->Modified();
  }

  //---------------------------------------------------------------------------
  // DeepCopy and MakeTransform: the copy shares the mesh, which is never modified
  std::cout << "Checking DeepCopy" << std::endl;
  {
    transform->SetDisplacementScale(2.0);
    transform->SetBorderFalloffDistance(15.0);
    transform->AutomaticBorderFalloffDistanceOn();
    vtkAbstractTransform* madeTransform = transform->MakeTransform();
    vtkMeshDisplacementTransform* madeMeshTransform = vtkMeshDisplacementTransform::SafeDownCast(madeTransform);
    CHECK_NOT_NULL(madeMeshTransform);
    madeTransform->DeepCopy(transform);
    CHECK_POINTER(madeMeshTransform->GetInputPointSet(), volumetricMesh.GetPointer());
    CHECK_STRING(madeMeshTransform->GetDisplacementArrayName(), "Displacement");
    CHECK_DOUBLE(madeMeshTransform->GetDisplacementScale(), 2.0);
    CHECK_DOUBLE(madeMeshTransform->GetBorderFalloffDistance(), 15.0);
    CHECK_BOOL(madeMeshTransform->GetAutomaticBorderFalloffDistance(), true);
    CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(madeTransform, volumetricMesh, "Displacement", 2.0), 0.0, 1e-9);
    madeTransform->Delete();
    transform->SetDisplacementScale(1.0);
    transform->SetBorderFalloffDistance(0.0);
    transform->AutomaticBorderFalloffDistanceOff();
  }

  //---------------------------------------------------------------------------
  // Around a thin shell with curved cells the falloff never gives a displacement
  // longer than the longest vector of the mesh: the displacement outside is taken
  // from the closest point of the mesh, whose interpolation weights must not be
  // extrapolated
  std::cout << "Checking the falloff around a thin curved shell" << std::endl;
  {
    vtkNew<vtkUnstructuredGrid> shell;
    CreateThinSaddleShell(shell);
    // A non-linear field: extrapolating it from a wrong cell gives a wrong displacement (a linear field
    // would be reproduced whichever cell it is interpolated in)
    vtkNew<vtkDoubleArray> displacementArray;
    displacementArray->SetName("Displacement");
    displacementArray->SetNumberOfComponents(3);
    displacementArray->SetNumberOfTuples(shell->GetNumberOfPoints());
    for (vtkIdType pointIndex = 0; pointIndex < shell->GetNumberOfPoints(); pointIndex++)
      {
      double point[3] = { 0.0, 0.0, 0.0 };
      shell->GetPoint(pointIndex, point);
      displacementArray->SetTuple3(pointIndex, 2.0*sin(1.3*point[0]), 2.0*cos(1.1*point[1]), 2.0*sin(0.7*point[0]*point[1]));
      }
    shell->GetPointData()->AddArray(displacementArray);
    double largestDisplacement = 0.0;
    for (vtkIdType pointIndex = 0; pointIndex < displacementArray->GetNumberOfTuples(); pointIndex++)
      {
      double vector[3] = { 0.0, 0.0, 0.0 };
      displacementArray->GetTuple(pointIndex, vector);
      largestDisplacement = std::max(largestDisplacement, vtkMath::Norm(vector));
      }
    vtkNew<vtkMeshDisplacementTransform> shellTransform;
    shellTransform->SetInputPointSet(shell);
    shellTransform->SetDisplacementArrayName("Displacement");
    shellTransform->AutomaticBorderFalloffDistanceOn();
    shellTransform->Update();
    const double reach = shellTransform->GetEffectiveBorderFalloffDistance();
    double bounds[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    shell->GetBounds(bounds);
    double longestDisplacement = 0.0;
    int numberOfMovedPositions = 0;
    const int samplesPerAxis = 40;
    for (int k = 0; k < samplesPerAxis; k++)
      {
      for (int j = 0; j < samplesPerAxis; j++)
        {
        for (int i = 0; i < samplesPerAxis; i++)
          {
          double position[3] = { 0.0, 0.0, 0.0 };
          const int indices[3] = { i, j, k };
          for (int axis = 0; axis < 3; axis++)
            {
            const double low = bounds[2*axis]-0.5*reach;
            const double high = bounds[2*axis+1]+0.5*reach;
            position[axis] = low+(high-low)*(indices[axis]+0.37)/samplesPerAxis;
            }
          double displacement[3] = { 0.0, 0.0, 0.0 };
          if (shellTransform->GetDisplacement(position, displacement))
            {
            numberOfMovedPositions++;
            }
          longestDisplacement = std::max(longestDisplacement, vtkMath::Norm(displacement));
          }
        }
      }
    std::cout << "  largest displacement of the mesh " << largestDisplacement << ", longest displacement around it "
              << longestDisplacement << ", " << numberOfMovedPositions << " positions moved" << std::endl;
    CHECK_BOOL(numberOfMovedPositions > 0, true);
    CHECK_BOOL(longestDisplacement <= largestDisplacement*(1.0+1e-9), true);
  }

  //---------------------------------------------------------------------------
  // Where the inverse does not converge (here because it may only take one
  // iteration), the failures are counted, but only the first one after a
  // modification is logged, unless all are requested
  std::cout << "Checking inverse convergence failures" << std::endl;
  {
    vtkNew<vtkMeshDisplacementTransform> limitedTransform;
    limitedTransform->SetInputPointSet(volumetricMesh);
    limitedTransform->SetDisplacementArrayName("Displacement");
    limitedTransform->SetBorderFalloffDistance(30.0);
    limitedTransform->SetInverseIterations(1);
    vtkAbstractTransform* limitedInverse = limitedTransform->GetInverse();
    // Positions around the border of the cube, where the displacement varies
    auto invertAround = [&]()
      {
      for (int i = 0; i < 20; i++)
        {
        for (int j = 0; j < 20; j++)
          {
          double target[3] = { 40.0+i*1.5, -40.0+j*4.0, 10.0 };
          double inverted[3] = { 0.0, 0.0, 0.0 };
          limitedInverse->TransformPoint(target, inverted);
          }
        }
      };
    TESTING_OUTPUT_RESET();
    invertAround();
    vtkMeshDisplacementTransform* limitedInverseMesh = vtkMeshDisplacementTransform::SafeDownCast(limitedInverse);
    CHECK_NOT_NULL(limitedInverseMesh);
    const vtkIdType numberOfFailures = limitedInverseMesh->GetNumberOfInverseConvergenceFailures();
    std::cout << "  " << numberOfFailures << " failures" << std::endl;
    CHECK_BOOL(numberOfFailures > 1, true);
    TESTING_OUTPUT_ASSERT_WARNINGS(1);
    TESTING_OUTPUT_RESET();
    // Again: already reported
    invertAround();
    TESTING_OUTPUT_ASSERT_WARNINGS(0);
    // Modified: reported again, and counted from zero
    limitedTransform->SetDisplacementScale(1.01);
    invertAround();
    TESTING_OUTPUT_ASSERT_WARNINGS(1);
    CHECK_BOOL(limitedInverseMesh->GetNumberOfInverseConvergenceFailures() <= numberOfFailures + 10, true);
    TESTING_OUTPUT_RESET();
    // All of them
    limitedTransform->LogAllInverseConvergenceFailuresOn();
    invertAround();
    CHECK_BOOL(limitedInverseMesh->GetLogAllInverseConvergenceFailures(), true);
    const vtkIdType numberOfLoggedFailures = limitedInverseMesh->GetNumberOfInverseConvergenceFailures();
    TESTING_OUTPUT_ASSERT_WARNINGS(static_cast<int>(numberOfLoggedFailures));
    TESTING_OUTPUT_RESET();
  }

  //---------------------------------------------------------------------------
  // Without usable displacements the transform is the identity
  std::cout << "Checking invalid inputs" << std::endl;
  {
    double point[3] = { 10.0, 20.0, 30.0 };

    // No input yet
    vtkNew<vtkMeshDisplacementTransform> invalidTransform;
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(invalidTransform, point, point), 0.0, 1e-12);

    // Displacement array that does not exist
    invalidTransform->SetInputPointSet(volumetricMesh);
    invalidTransform->SetDisplacementArrayName("NoSuchArray");
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(invalidTransform, point, point), 0.0, 1e-12);

    // Displacement array with the wrong number of components
    vtkNew<vtkUnstructuredGrid> scalarMesh;
    CreateHexahedralGrid(scalarMesh, 3, 50.0);
    vtkNew<vtkDoubleArray> scalarArray;
    scalarArray->SetName("Displacement");
    scalarArray->SetNumberOfComponents(1);
    scalarArray->SetNumberOfTuples(scalarMesh->GetNumberOfPoints());
    scalarArray->FillComponent(0, 5.0);
    scalarMesh->GetPointData()->AddArray(scalarArray);
    invalidTransform->SetInputPointSet(scalarMesh);
    invalidTransform->SetDisplacementArrayName("Displacement");
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(invalidTransform, point, point), 0.0, 1e-12);

    // Empty mesh
    vtkNew<vtkPolyData> emptyMesh;
    vtkNew<vtkPoints> emptyPoints;
    emptyMesh->SetPoints(emptyPoints);
    invalidTransform->SetInputPointSet(emptyMesh);
    CHECK_DOUBLE_TOLERANCE(GetErrorAtPoint(invalidTransform, point, point), 0.0, 1e-12);
  }

  // None of the above is reported as a warning or error
  TESTING_OUTPUT_ASSERT_WARNINGS_ERRORS(0);

  std::cout << "vtkMeshDisplacementTransformTest1 completed successfully" << std::endl;
  return EXIT_SUCCESS;
}
