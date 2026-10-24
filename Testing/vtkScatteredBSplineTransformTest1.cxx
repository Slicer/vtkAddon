/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// vtkAddon includes
#include <vtkAddonTestingMacros.h>
#include <vtkScatteredBSplineTransform.h>
#include <vtkTestingOutputWindow.h>

// VTK includes
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkImageData.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>

// STD includes
#include <iostream>

namespace
{

//----------------------------------------------------------------------------
// Smooth analytic warp that is used as ground truth throughout the test.
void GetReferenceDisplacement(const double point[3], double displacement[3])
{
  const double amplitude = 3.0;
  const double frequency = 0.05;
  displacement[0] = amplitude*sin(frequency*point[1]);
  displacement[1] = amplitude*sin(frequency*point[2]);
  displacement[2] = amplitude*sin(frequency*point[0]);
}

//----------------------------------------------------------------------------
// Build a volumetric mesh: a regular grid of hexahedral cells that fills the
// cube of the given half-size.
void CreateHexahedralGrid(vtkUnstructuredGrid* grid, int resolution, double halfSize)
{
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(resolution*resolution*resolution);
  for (int k = 0; k < resolution; k++)
    {
    for (int j = 0; j < resolution; j++)
      {
      for (int i = 0; i < resolution; i++)
        {
        const double step = 2.0*halfSize/(resolution-1);
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
// displacement array prescribes.
double GetMaximumErrorAtMeshPoints(vtkAbstractTransform* transform, vtkPointSet* mesh,
                                   const char* arrayName)
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
    double expected[3] = { point[0]+displacement[0], point[1]+displacement[1], point[2]+displacement[2] };
    double actual[3] = { 0.0, 0.0, 0.0 };
    transform->TransformPoint(point, actual);
    maximumError = std::max(maximumError, sqrt(vtkMath::Distance2BetweenPoints(expected, actual)));
    }
  return maximumError;
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkScatteredBSplineTransformTest1(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  TESTING_OUTPUT_INIT();

  const double tolerance = 0.01;

  //---------------------------------------------------------------------------
  // Fit a surface mesh
  std::cout << "Fitting a surface mesh" << std::endl;

  vtkNew<vtkSphereSource> sphereSource;
  sphereSource->SetRadius(50.0);
  sphereSource->SetThetaResolution(30);
  sphereSource->SetPhiResolution(30);
  sphereSource->Update();
  vtkNew<vtkPolyData> surfaceMesh;
  surfaceMesh->DeepCopy(sphereSource->GetOutput());
  AddDisplacementArray(surfaceMesh, "Displacement");

  vtkNew<vtkScatteredBSplineTransform> transform;
  transform->SetInputPointSet(surfaceMesh);
  transform->SetDisplacementArrayName("Displacement");
  transform->SetTolerance(tolerance);
  transform->Update();

  std::cout << "  levels: " << transform->GetNumberOfLevels()
            << ", residual: " << transform->GetResidual() << std::endl;
  CHECK_BOOL(transform->GetNumberOfLevels() > 0, true);
  CHECK_BOOL(transform->GetResidual() <= tolerance, true);
  // The residual that the transform reports must match what the transform
  // actually does at the mesh points
  CHECK_DOUBLE_TOLERANCE(GetMaximumErrorAtMeshPoints(transform, surfaceMesh, "Displacement"), 0.0, tolerance);

  //---------------------------------------------------------------------------
  // The fitted transform must approximate the reference warp between the mesh
  // points, too, not only at them
  std::cout << "Checking interpolation between mesh points" << std::endl;

  double maximumInterpolationError = 0.0;
  for (int i = 0; i < 10; i++)
    {
    for (int j = 0; j < 10; j++)
      {
      // A point on the sphere that is not one of the mesh points
      const double theta = (i+0.5)*vtkMath::Pi()/10.0;
      const double phi = (j+0.5)*2.0*vtkMath::Pi()/10.0;
      double point[3] = { 50.0*sin(theta)*cos(phi), 50.0*sin(theta)*sin(phi), 50.0*cos(theta) };
      double displacement[3] = { 0.0, 0.0, 0.0 };
      GetReferenceDisplacement(point, displacement);
      double expected[3] = { point[0]+displacement[0], point[1]+displacement[1], point[2]+displacement[2] };
      double actual[3] = { 0.0, 0.0, 0.0 };
      transform->TransformPoint(point, actual);
      maximumInterpolationError = std::max(maximumInterpolationError,
        sqrt(vtkMath::Distance2BetweenPoints(expected, actual)));
      }
    }
  std::cout << "  maximum interpolation error: " << maximumInterpolationError << std::endl;
  CHECK_BOOL(maximumInterpolationError < 0.5, true);

  //---------------------------------------------------------------------------
  // The inverse must undo the forward transform
  std::cout << "Checking the inverse transform" << std::endl;

  double maximumInverseError = 0.0;
  vtkPoints* surfacePoints = surfaceMesh->GetPoints();
  for (vtkIdType pointIndex = 0; pointIndex < surfacePoints->GetNumberOfPoints(); pointIndex += 17)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    surfacePoints->GetPoint(pointIndex, point);
    double transformedPoint[3] = { 0.0, 0.0, 0.0 };
    transform->TransformPoint(point, transformedPoint);
    double roundTrippedPoint[3] = { 0.0, 0.0, 0.0 };
    transform->GetInverse()->TransformPoint(transformedPoint, roundTrippedPoint);
    maximumInverseError = std::max(maximumInverseError,
      sqrt(vtkMath::Distance2BetweenPoints(point, roundTrippedPoint)));
    }
  std::cout << "  maximum inverse error: " << maximumInverseError << std::endl;
  CHECK_BOOL(maximumInverseError < 1e-4, true);

  //---------------------------------------------------------------------------
  // Fit a volumetric mesh
  std::cout << "Fitting a volumetric mesh" << std::endl;

  vtkNew<vtkUnstructuredGrid> volumetricMesh;
  CreateHexahedralGrid(volumetricMesh, 10, 50.0);
  AddDisplacementArray(volumetricMesh, "Displacement");

  vtkNew<vtkScatteredBSplineTransform> volumetricTransform;
  volumetricTransform->SetInputPointSet(volumetricMesh);
  volumetricTransform->SetDisplacementArrayName("Displacement");
  volumetricTransform->SetTolerance(tolerance);
  volumetricTransform->Update();

  std::cout << "  levels: " << volumetricTransform->GetNumberOfLevels()
            << ", residual: " << volumetricTransform->GetResidual() << std::endl;
  CHECK_BOOL(volumetricTransform->GetResidual() <= tolerance, true);
  CHECK_DOUBLE_TOLERANCE(
    GetMaximumErrorAtMeshPoints(volumetricTransform, volumetricMesh, "Displacement"), 0.0, tolerance);

  //---------------------------------------------------------------------------
  // If no displacement array name is specified then the active vectors array
  // is used
  std::cout << "Checking the active vectors array fallback" << std::endl;

  vtkNew<vtkPolyData> activeVectorsMesh;
  activeVectorsMesh->DeepCopy(sphereSource->GetOutput());
  AddDisplacementArray(activeVectorsMesh, "Displacement", true);

  vtkNew<vtkScatteredBSplineTransform> activeVectorsTransform;
  activeVectorsTransform->SetInputPointSet(activeVectorsMesh);
  activeVectorsTransform->SetTolerance(tolerance);
  activeVectorsTransform->Update();
  CHECK_DOUBLE_TOLERANCE(activeVectorsTransform->GetResidual(), transform->GetResidual(), 1e-9);

  //---------------------------------------------------------------------------
  // Linear approximation must not degrade the fit at the mesh points
  std::cout << "Checking the linear initial approximation" << std::endl;

  vtkNew<vtkScatteredBSplineTransform> linearApproximationTransform;
  linearApproximationTransform->SetInputPointSet(surfaceMesh);
  linearApproximationTransform->SetDisplacementArrayName("Displacement");
  linearApproximationTransform->SetTolerance(tolerance);
  linearApproximationTransform->UseLinearApproximationOn();
  linearApproximationTransform->Update();
  std::cout << "  levels: " << linearApproximationTransform->GetNumberOfLevels()
            << ", residual: " << linearApproximationTransform->GetResidual() << std::endl;
  CHECK_BOOL(linearApproximationTransform->GetResidual() <= tolerance, true);
  CHECK_DOUBLE_TOLERANCE(
    GetMaximumErrorAtMeshPoints(linearApproximationTransform, surfaceMesh, "Displacement"), 0.0, tolerance);

  //---------------------------------------------------------------------------
  // Explicit control point spacing must be reflected in the fitted lattice and
  // must not be altered by the fit
  std::cout << "Checking explicit control point spacing" << std::endl;

  // The volumetric mesh spans exactly 100mm along each axis, so a 25mm request
  // results in 5 control points and therefore exactly 25mm spacing
  vtkNew<vtkScatteredBSplineTransform> coarseTransform;
  coarseTransform->SetInputPointSet(volumetricMesh);
  coarseTransform->SetDisplacementArrayName("Displacement");
  coarseTransform->SetControlPointSpacing(25.0, 25.0, 25.0);
  coarseTransform->SetMaximumNumberOfLevels(1);
  coarseTransform->Update();
  CHECK_INT(coarseTransform->GetNumberOfLevels(), 1);
  CHECK_NOT_NULL(coarseTransform->GetCoefficientData());
  double fittedSpacing[3] = { 0.0, 0.0, 0.0 };
  coarseTransform->GetCoefficientData()->GetSpacing(fittedSpacing);
  std::cout << "  fitted lattice spacing: " << fittedSpacing[0] << " "
            << fittedSpacing[1] << " " << fittedSpacing[2] << std::endl;
  for (int axis = 0; axis < 3; axis++)
    {
    CHECK_DOUBLE_TOLERANCE(fittedSpacing[axis], 25.0, 0.01);
    // The requested spacing must not have been overwritten by the fit
    CHECK_DOUBLE_TOLERANCE(coarseTransform->GetControlPointSpacing()[axis], 25.0, 1e-12);
    }

  // The minimum control point spacing must stop the refinement
  vtkNew<vtkScatteredBSplineTransform> spacingLimitedTransform;
  spacingLimitedTransform->SetInputPointSet(volumetricMesh);
  spacingLimitedTransform->SetDisplacementArrayName("Displacement");
  spacingLimitedTransform->SetControlPointSpacing(25.0, 25.0, 25.0);
  spacingLimitedTransform->SetMinimumControlPointSpacing(13.0);
  spacingLimitedTransform->SetTolerance(1e-12);
  spacingLimitedTransform->Update();
  std::cout << "  levels with a 13mm minimum spacing: "
            << spacingLimitedTransform->GetNumberOfLevels() << std::endl;
  CHECK_INT(spacingLimitedTransform->GetNumberOfLevels(), 1);

  //---------------------------------------------------------------------------
  // Points outside an explicitly specified fitting domain are ignored
  std::cout << "Checking explicit grid bounds" << std::endl;

  vtkNew<vtkScatteredBSplineTransform> boundedTransform;
  boundedTransform->SetInputPointSet(volumetricMesh);
  boundedTransform->SetDisplacementArrayName("Displacement");
  boundedTransform->AutomaticGridBoundsOff();
  boundedTransform->SetGridBounds(-30.0, 30.0, -30.0, 30.0, -60.0, 60.0);
  TESTING_OUTPUT_ASSERT_WARNINGS_BEGIN();
  boundedTransform->Update();
  TESTING_OUTPUT_ASSERT_WARNINGS_END();
  double coefficientOrigin[3] = { 0.0, 0.0, 0.0 };
  boundedTransform->GetCoefficientData()->GetOrigin(coefficientOrigin);
  // The lattice starts one control point below the requested domain
  CHECK_BOOL(coefficientOrigin[0] < -30.0, true);
  CHECK_BOOL(boundedTransform->GetResidual() >= 0.0, true);

  //---------------------------------------------------------------------------
  // The automatic fitting domain can extend beyond the points by a margin
  std::cout << "Checking the grid bounds margin" << std::endl;

  vtkNew<vtkScatteredBSplineTransform> marginTransform;
  marginTransform->SetInputPointSet(volumetricMesh);
  marginTransform->SetDisplacementArrayName("Displacement");
  CHECK_DOUBLE(marginTransform->GetGridBoundsMargin(), 0.0);
  marginTransform->Update();
  double originWithoutMargin[3] = { 0.0, 0.0, 0.0 };
  marginTransform->GetCoefficientData()->GetOrigin(originWithoutMargin);
  double spacingWithoutMargin[3] = { 0.0, 0.0, 0.0 };
  marginTransform->GetCoefficientData()->GetSpacing(spacingWithoutMargin);
  marginTransform->SetGridBoundsMargin(20.0);
  marginTransform->Update();
  double originWithMargin[3] = { 0.0, 0.0, 0.0 };
  marginTransform->GetCoefficientData()->GetOrigin(originWithMargin);
  double meshBounds[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  volumetricMesh->GetBounds(meshBounds);
  // The lattice starts one control point below the domain, which starts 20 mm below the points
  CHECK_BOOL(originWithMargin[0] < meshBounds[0] - 20.0, true);
  CHECK_BOOL(originWithMargin[0] < originWithoutMargin[0], true);
  CHECK_BOOL(marginTransform->GetResidual() >= 0.0, true);
  // Copied with the transform
  vtkNew<vtkScatteredBSplineTransform> marginCopy;
  marginCopy->DeepCopy(marginTransform);
  CHECK_DOUBLE(marginCopy->GetGridBoundsMargin(), 20.0);
  // Negative margins are not allowed
  marginTransform->SetGridBoundsMargin(-5.0);
  CHECK_DOUBLE(marginTransform->GetGridBoundsMargin(), 0.0);

  //---------------------------------------------------------------------------
  // Modifying the input mesh must trigger a new fit
  std::cout << "Checking that the transform follows changes of the input mesh" << std::endl;

  double firstPoint[3] = { 0.0, 0.0, 0.0 };
  surfaceMesh->GetPoints()->GetPoint(0, firstPoint);
  double transformedFirstPoint[3] = { 0.0, 0.0, 0.0 };
  transform->TransformPoint(firstPoint, transformedFirstPoint);

  vtkDataArray* displacementArray = surfaceMesh->GetPointData()->GetArray("Displacement");
  double modifiedDisplacement[3] = { 20.0, -20.0, 20.0 };
  displacementArray->SetTuple(0, modifiedDisplacement);
  displacementArray->Modified();
  surfaceMesh->Modified();

  transform->Update();
  double transformedFirstPointAfterChange[3] = { 0.0, 0.0, 0.0 };
  transform->TransformPoint(firstPoint, transformedFirstPointAfterChange);
  std::cout << "  transformed point moved by "
            << sqrt(vtkMath::Distance2BetweenPoints(transformedFirstPoint, transformedFirstPointAfterChange))
            << std::endl;
  CHECK_BOOL(vtkMath::Distance2BetweenPoints(
    transformedFirstPoint, transformedFirstPointAfterChange) > 1.0, true);

  // Restore the original displacement
  double originalDisplacement[3] = { 0.0, 0.0, 0.0 };
  GetReferenceDisplacement(firstPoint, originalDisplacement);
  displacementArray->SetTuple(0, originalDisplacement);
  displacementArray->Modified();
  surfaceMesh->Modified();
  transform->Update();

  //---------------------------------------------------------------------------
  // DeepCopy and MakeTransform
  std::cout << "Checking DeepCopy" << std::endl;

  vtkAbstractTransform* madeTransform = transform->MakeTransform();
  CHECK_NOT_NULL(vtkScatteredBSplineTransform::SafeDownCast(madeTransform));
  madeTransform->DeepCopy(transform);
  madeTransform->Update();
  CHECK_DOUBLE_TOLERANCE(
    GetMaximumErrorAtMeshPoints(madeTransform, surfaceMesh, "Displacement"), 0.0, tolerance);
  madeTransform->Delete();

  //---------------------------------------------------------------------------
  // Invalid inputs must be reported and must leave the transform as identity
  std::cout << "Checking error handling" << std::endl;

  {
    // No input yet: the identity, quietly, as a transform that is being set up is not an error
    vtkNew<vtkScatteredBSplineTransform> invalidTransform;
    invalidTransform->Update();
    CHECK_DOUBLE(invalidTransform->GetResidual(), -1.0);
    double point[3] = { 10.0, 20.0, 30.0 };
    double transformedPoint[3] = { 0.0, 0.0, 0.0 };
    invalidTransform->TransformPoint(point, transformedPoint);
    CHECK_DOUBLE_TOLERANCE(sqrt(vtkMath::Distance2BetweenPoints(point, transformedPoint)), 0.0, 1e-9);
  }

  {
    // Displacement array with the wrong number of components
    vtkNew<vtkPolyData> invalidMesh;
    invalidMesh->DeepCopy(sphereSource->GetOutput());
    vtkNew<vtkDoubleArray> scalarArray;
    scalarArray->SetName("Displacement");
    scalarArray->SetNumberOfComponents(1);
    scalarArray->SetNumberOfTuples(invalidMesh->GetNumberOfPoints());
    scalarArray->FillComponent(0, 0.0);
    invalidMesh->GetPointData()->AddArray(scalarArray);

    vtkNew<vtkScatteredBSplineTransform> invalidTransform;
    invalidTransform->SetInputPointSet(invalidMesh);
    invalidTransform->SetDisplacementArrayName("Displacement");
    TESTING_OUTPUT_ASSERT_ERRORS_BEGIN();
    invalidTransform->Update();
    TESTING_OUTPUT_ASSERT_ERRORS_END();
  }

  {
    // Displacement array that does not exist
    vtkNew<vtkScatteredBSplineTransform> invalidTransform;
    invalidTransform->SetInputPointSet(surfaceMesh);
    invalidTransform->SetDisplacementArrayName("NoSuchArray");
    TESTING_OUTPUT_ASSERT_ERRORS_BEGIN();
    invalidTransform->Update();
    TESTING_OUTPUT_ASSERT_ERRORS_END();
  }

  {
    // Control point spacing so small that the lattice would not fit in memory
    vtkNew<vtkScatteredBSplineTransform> invalidTransform;
    invalidTransform->SetInputPointSet(surfaceMesh);
    invalidTransform->SetDisplacementArrayName("Displacement");
    invalidTransform->SetControlPointSpacing(1e-6, 1e-6, 1e-6);
    TESTING_OUTPUT_ASSERT_ERRORS_BEGIN();
    invalidTransform->Update();
    TESTING_OUTPUT_ASSERT_ERRORS_END();
  }

  {
    // Empty mesh
    vtkNew<vtkPolyData> emptyMesh;
    vtkNew<vtkPoints> emptyPoints;
    emptyMesh->SetPoints(emptyPoints);
    vtkNew<vtkScatteredBSplineTransform> invalidTransform;
    invalidTransform->SetInputPointSet(emptyMesh);
    TESTING_OUTPUT_ASSERT_ERRORS_BEGIN();
    invalidTransform->Update();
    TESTING_OUTPUT_ASSERT_ERRORS_END();
  }

  std::cout << "vtkScatteredBSplineTransformTest1 completed successfully" << std::endl;
  return EXIT_SUCCESS;
}
