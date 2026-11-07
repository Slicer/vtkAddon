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

// Checks that vtkMeshDisplacementTransform approximates a non-linear
// displacement field (a B-spline transform) on hexahedral and tetrahedral
// meshes, with an error that decreases at the expected rate as the mesh is
// refined.

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
#include <iostream>
#include <vector>

namespace
{

// Half-size of the cube that the meshes fill
const double MeshHalfSize = 50.0;

//----------------------------------------------------------------------------
// Smooth, non-linear reference transform: a cubic B-spline transform whose
// coefficient grid extends beyond the meshes, with displacements of a few mm
// that vary on a length scale comparable to the mesh size.
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
// Volumetric mesh that fills the cube: a regular grid of hexahedra, or of
// tetrahedra (each hexahedron split into 6 that share the main diagonal, so
// that neighboring cells share faces).
void CreateGrid(vtkUnstructuredGrid* grid, int resolution, bool tetrahedra)
{
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(resolution*resolution*resolution);
  const double step = 2.0*MeshHalfSize/(resolution-1);
  for (int k = 0; k < resolution; k++)
    {
    for (int j = 0; j < resolution; j++)
      {
      for (int i = 0; i < resolution; i++)
        {
        points->SetPoint((k*resolution+j)*resolution+i,
          -MeshHalfSize+i*step, -MeshHalfSize+j*step, -MeshHalfSize+k*step);
        }
      }
    }
  grid->SetPoints(points);

  // Point offsets of the corners of a hexahedron, indexed by (x bit, y bit, z bit)
  const vtkIdType xOffset = 1;
  const vtkIdType yOffset = resolution;
  const vtkIdType zOffset = resolution*resolution;
  // Corners of the 6 tetrahedra, as bit masks of the corner: each goes from
  // corner 0 to corner 7 by changing one coordinate at a time
  const int tetrahedronCorners[6][4] = {
    { 0, 1, 3, 7 }, { 0, 1, 5, 7 }, { 0, 2, 3, 7 },
    { 0, 2, 6, 7 }, { 0, 4, 5, 7 }, { 0, 4, 6, 7 } };

  grid->Allocate((resolution-1)*(resolution-1)*(resolution-1)*(tetrahedra ? 6 : 1));
  vtkNew<vtkHexahedron> hexahedron;
  vtkNew<vtkTetra> tetrahedron;
  for (int k = 0; k < resolution-1; k++)
    {
    for (int j = 0; j < resolution-1; j++)
      {
      for (int i = 0; i < resolution-1; i++)
        {
        const vtkIdType origin = (k*resolution+j)*resolution+i;
        vtkIdType cornerIds[8];
        for (int corner = 0; corner < 8; corner++)
          {
          cornerIds[corner] = origin
            + ((corner & 1) ? xOffset : 0) + ((corner & 2) ? yOffset : 0) + ((corner & 4) ? zOffset : 0);
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
          // VTK hexahedron point order
          const int hexahedronCorners[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
          for (int cornerIndex = 0; cornerIndex < 8; cornerIndex++)
            {
            hexahedron->GetPointIds()->SetId(cornerIndex, cornerIds[hexahedronCorners[cornerIndex]]);
            }
          grid->InsertNextCell(hexahedron->GetCellType(), hexahedron->GetPointIds());
          }
        }
      }
    }
}

//----------------------------------------------------------------------------
// Add the displacement of the reference transform at each mesh point
void AddDisplacementArray(vtkPointSet* mesh, vtkAbstractTransform* referenceTransform)
{
  vtkPoints* points = mesh->GetPoints();
  vtkNew<vtkDoubleArray> displacementArray;
  displacementArray->SetName("Displacement");
  displacementArray->SetNumberOfComponents(3);
  displacementArray->SetNumberOfTuples(points->GetNumberOfPoints());
  for (vtkIdType pointIndex = 0; pointIndex < points->GetNumberOfPoints(); pointIndex++)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    points->GetPoint(pointIndex, point);
    double transformedPoint[3] = { 0.0, 0.0, 0.0 };
    referenceTransform->TransformPoint(point, transformedPoint);
    double displacement[3] = { transformedPoint[0]-point[0], transformedPoint[1]-point[1], transformedPoint[2]-point[2] };
    displacementArray->SetTuple(pointIndex, displacement);
    }
  mesh->GetPointData()->AddArray(displacementArray);
}

//----------------------------------------------------------------------------
// Positions where the transforms are compared: spread over the inside of the
// meshes, none of them a mesh point of any of the meshes used
std::vector<std::vector<double>> GetTestPoints()
{
  std::vector<std::vector<double>> testPoints;
  for (int k = 0; k < 7; k++)
    {
    for (int j = 0; j < 7; j++)
      {
      for (int i = 0; i < 7; i++)
        {
        testPoints.push_back({ -47.0+i*13.9, -45.3+j*14.3, -46.1+k*13.7 });
        }
      }
    }
  return testPoints;
}

//----------------------------------------------------------------------------
// Maximum error at the mesh points, root mean square errors at the test points
struct Errors
{
  double MeshPoints = 0.0;
  double Interpolation = 0.0;
  double Derivative = 0.0;
  double Inverse = 0.0;
};

//----------------------------------------------------------------------------
// Compare the mesh transform with the reference transform
Errors GetErrors(vtkMeshDisplacementTransform* transform, vtkBSplineTransform* referenceTransform)
{
  Errors errors;

  // At the mesh points
  vtkPoints* points = transform->GetInputPointSet()->GetPoints();
  for (vtkIdType pointIndex = 0; pointIndex < points->GetNumberOfPoints(); pointIndex++)
    {
    double point[3] = { 0.0, 0.0, 0.0 };
    points->GetPoint(pointIndex, point);
    double expected[3] = { 0.0, 0.0, 0.0 };
    referenceTransform->TransformPoint(point, expected);
    double actual[3] = { 0.0, 0.0, 0.0 };
    transform->TransformPoint(point, actual);
    errors.MeshPoints = std::max(errors.MeshPoints, sqrt(vtkMath::Distance2BetweenPoints(expected, actual)));
    }

  // Between the mesh points (root mean square errors, which depend less on
  // where single test points are in their cells than the maximum)
  int numberOfTestPoints = 0;
  int numberOfInverseTestPoints = 0;
  for (const std::vector<double>& testPoint : GetTestPoints())
    {
    const double* point = testPoint.data();
    double expected[3] = { 0.0, 0.0, 0.0 };
    double expectedDerivative[3][3];
    referenceTransform->InternalTransformDerivative(point, expected, expectedDerivative);
    double actual[3] = { 0.0, 0.0, 0.0 };
    double actualDerivative[3][3];
    transform->InternalTransformDerivative(point, actual, actualDerivative);
    errors.Interpolation += vtkMath::Distance2BetweenPoints(expected, actual);
    for (int row = 0; row < 3; row++)
      {
      for (int column = 0; column < 3; column++)
        {
        errors.Derivative += (expectedDerivative[row][column]-actualDerivative[row][column])
          * (expectedDerivative[row][column]-actualDerivative[row][column]) / 9.0;
        }
      }
    numberOfTestPoints++;

    // The inverse of the mesh transform maps the reference transformed point
    // back. Only checked where the transformed point is inside the mesh, as
    // elsewhere the mesh transform does not move anything.
    if (std::max({ fabs(expected[0]), fabs(expected[1]), fabs(expected[2]) }) > MeshHalfSize - 5.0)
      {
      continue;
      }
    double roundTrippedPoint[3] = { 0.0, 0.0, 0.0 };
    transform->GetInverse()->TransformPoint(expected, roundTrippedPoint);
    errors.Inverse += vtkMath::Distance2BetweenPoints(point, roundTrippedPoint);
    numberOfInverseTestPoints++;
    }
  errors.Interpolation = sqrt(errors.Interpolation / numberOfTestPoints);
  errors.Derivative = sqrt(errors.Derivative / numberOfTestPoints);
  errors.Inverse = sqrt(errors.Inverse / numberOfInverseTestPoints);

  return errors;
}

//----------------------------------------------------------------------------
// Check the errors of meshes of increasing resolution with one cell type
int TestRefinement(vtkBSplineTransform* referenceTransform, bool tetrahedra)
{
  std::cout << "Checking " << (tetrahedra ? "tetrahedral" : "hexahedral") << " meshes" << std::endl;

  // Cell sizes 25, 12.5 and 6.25 mm
  const int resolutions[3] = { 5, 9, 17 };
  Errors errors[3];
  for (int resolutionIndex = 0; resolutionIndex < 3; resolutionIndex++)
    {
    vtkNew<vtkUnstructuredGrid> mesh;
    CreateGrid(mesh, resolutions[resolutionIndex], tetrahedra);
    AddDisplacementArray(mesh, referenceTransform);

    vtkNew<vtkMeshDisplacementTransform> transform;
    transform->SetInputPointSet(mesh);
    transform->SetDisplacementArrayName("Displacement");

    errors[resolutionIndex] = GetErrors(transform, referenceTransform);
    std::cout << "  cell size " << 2.0*MeshHalfSize/(resolutions[resolutionIndex]-1) << " mm:"
              << " mesh point error " << errors[resolutionIndex].MeshPoints
              << ", interpolation error " << errors[resolutionIndex].Interpolation
              << ", derivative error " << errors[resolutionIndex].Derivative
              << ", inverse error " << errors[resolutionIndex].Inverse << std::endl;

    // The mesh points move exactly as the reference transform moves them
    CHECK_DOUBLE_TOLERANCE(errors[resolutionIndex].MeshPoints, 0.0, 1e-9);

    if (resolutionIndex == 2)
      {
      // With a falloff outside the mesh the inverse also works for the test
      // points that the reference transform moves out of the mesh
      transform->SetBorderFalloffDistance(30.0);
      double maximumRoundTripError = 0.0;
      for (const std::vector<double>& testPoint : GetTestPoints())
        {
        double transformedPoint[3] = { 0.0, 0.0, 0.0 };
        transform->TransformPoint(testPoint.data(), transformedPoint);
        double roundTrippedPoint[3] = { 0.0, 0.0, 0.0 };
        transform->GetInverse()->TransformPoint(transformedPoint, roundTrippedPoint);
        maximumRoundTripError = std::max(maximumRoundTripError,
          sqrt(vtkMath::Distance2BetweenPoints(testPoint.data(), roundTrippedPoint)));
        }
      std::cout << "  with falloff: maximum round trip error " << maximumRoundTripError << std::endl;
      CHECK_BOOL(maximumRoundTripError < 1e-2, true);
      }
    }

  // The field is really non-linear: the coarsest mesh cannot reproduce it
  CHECK_BOOL(errors[0].Interpolation > 0.1, true);

  for (int resolutionIndex = 1; resolutionIndex < 3; resolutionIndex++)
    {
    // Linear interpolation is second order accurate: halving the cell size
    // divides the error by about 4
    CHECK_BOOL(errors[resolutionIndex].Interpolation < errors[resolutionIndex-1].Interpolation/3.0, true);
    // The derivative is first order accurate: halving the cell size divides
    // the error by about 2
    CHECK_BOOL(errors[resolutionIndex].Derivative < errors[resolutionIndex-1].Derivative/1.5, true);
    // The inverse is as accurate as the forward transform
    CHECK_BOOL(errors[resolutionIndex].Inverse < errors[resolutionIndex-1].Inverse/3.0, true);
    }

  // On the finest mesh the transforms are close
  CHECK_BOOL(errors[2].Interpolation < 0.05, true);
  CHECK_BOOL(errors[2].Derivative < 0.01, true);
  CHECK_BOOL(errors[2].Inverse < 0.05, true);

  return EXIT_SUCCESS;
}

} // end anonymous namespace

//----------------------------------------------------------------------------
int vtkMeshDisplacementTransformTest2(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  TESTING_OUTPUT_INIT();

  vtkNew<vtkBSplineTransform> referenceTransform;
  CreateReferenceTransform(referenceTransform);

  CHECK_EXIT_SUCCESS(TestRefinement(referenceTransform, false));
  CHECK_EXIT_SUCCESS(TestRefinement(referenceTransform, true));

  // None of the above is reported as a warning or error
  TESTING_OUTPUT_ASSERT_WARNINGS_ERRORS(0);

  std::cout << "vtkMeshDisplacementTransformTest2 completed successfully" << std::endl;
  return EXIT_SUCCESS;
}
