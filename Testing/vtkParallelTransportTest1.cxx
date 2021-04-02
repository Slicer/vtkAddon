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
#include <vtkParallelTransportFrame.h>

// VTK includes
#include <vtkArcSource.h>
#include <vtkDoubleArray.h>
#include <vtkLineSource.h>
#include <vtkNew.h>
#include <vtkPointData.h>

//----------------------------------------------------------------------------
int vtkParallelTransportTest1(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  vtkNew<vtkParallelTransportFrame> parallelTransportFrame;
  double tolerance = 1e-3;
  bool verbose = false;

  // Test computation with a straight line input (zero curvature)

  double lineDirection[3] = { 0.0, 0.0, 1.0 };
  vtkNew<vtkLineSource> line;
  line->SetPoint1(10,20,30);
  line->SetPoint1(30,70,-80);
  line->SetResolution(10);
  parallelTransportFrame->SetInputConnection(line->GetOutputPort());
  parallelTransportFrame->Update();

  {
    vtkPointData* pointData = parallelTransportFrame->GetOutput()->GetPointData();
    vtkDoubleArray* normalsArray = vtkDoubleArray::SafeDownCast(pointData->GetArray("Normals"));
    vtkDoubleArray* binormalsArray = vtkDoubleArray::SafeDownCast(pointData->GetArray("Binormals"));
    vtkDoubleArray* tangentsArray = vtkDoubleArray::SafeDownCast(pointData->GetArray("Tangents"));
    CHECK_NOT_NULL(normalsArray);
    CHECK_NOT_NULL(binormalsArray);
    CHECK_NOT_NULL(tangentsArray);
    for (vtkIdType tupleIndex = 0; tupleIndex < normalsArray->GetNumberOfTuples(); tupleIndex++)
      {
      double* v = normalsArray->GetTuple3(tupleIndex);
      CHECK_DOUBLE_TOLERANCE(v[0], 0.963584, tolerance);
      CHECK_DOUBLE_TOLERANCE(v[1], -0.176089, tolerance);
      CHECK_DOUBLE_TOLERANCE(v[2], 0.201244, tolerance);
      v = binormalsArray->GetTuple3(tupleIndex);
      CHECK_DOUBLE_TOLERANCE(v[0], 0.0, tolerance);
      CHECK_DOUBLE_TOLERANCE(v[1], 0.752577, tolerance);
      CHECK_DOUBLE_TOLERANCE(v[2], 0.658505, tolerance);
      v = tangentsArray->GetTuple3(tupleIndex);
      CHECK_DOUBLE_TOLERANCE(v[0], -0.267407, tolerance);
      CHECK_DOUBLE_TOLERANCE(v[1], -0.634524, tolerance);
      CHECK_DOUBLE_TOLERANCE(v[2], 0.725171, tolerance);
      }
  }

  // Test with an arc section

  vtkNew<vtkArcSource> arc;
  double arcCenter[3] = { 40, 50, 60 };
  double arcNormal[3] = { 0.3, -0.4, 0.9 };
  vtkMath::Normalize(arcNormal);
  arc->SetAngle(90);
  arc->SetCenter(arcCenter);
  arc->SetPolarVector(5, 15, -3);
  arc->UseNormalAndAngleOn();
  arc->SetResolution(200);

  // Make the initial radial direction as preferred initial normal direction.
  // This will make the binormal direction parallel with the arc's plane normal,
  // which will remain constant over the entire arc.
  arc->Update();
  double arcFirstPoint[3] = { 0,0,0 };
  arc->GetOutput()->GetPoint(0, arcFirstPoint);
  double preferredNormalDirection[3] = { arcCenter[0] - arcFirstPoint[0], arcCenter[1] - arcFirstPoint[1], arcCenter[2] - arcFirstPoint[2] };
  parallelTransportFrame->SetPreferredInitialNormalVector(preferredNormalDirection);
  
  parallelTransportFrame->SetInputConnection(arc->GetOutputPort());
  parallelTransportFrame->Update();

  {
    vtkPoints* points = parallelTransportFrame->GetOutput()->GetPoints();
    vtkPointData* pointData = parallelTransportFrame->GetOutput()->GetPointData();
    vtkDoubleArray* normalsArray = vtkDoubleArray::SafeDownCast(pointData->GetArray("Normals"));
    vtkDoubleArray* binormalsArray = vtkDoubleArray::SafeDownCast(pointData->GetArray("Binormals"));
    vtkDoubleArray* tangentsArray = vtkDoubleArray::SafeDownCast(pointData->GetArray("Tangents"));
    double constantBinormal[3] = { 0.0, 0.0, 0.0 };
    binormalsArray->GetTypedTuple(0, constantBinormal);
    for (vtkIdType tupleIndex = 1; tupleIndex < normalsArray->GetNumberOfTuples(); tupleIndex++)
      {
      double normal[3] = { 0.0, 0.0, 0.0 };
      normalsArray->GetTypedTuple(tupleIndex, normal);
      double binormal[3] = { 0.0, 0.0, 0.0 };
      binormalsArray->GetTypedTuple(tupleIndex, binormal);

      double pointPosition[3];
      points->GetPoint(tupleIndex, pointPosition);
      double radialDirection[3] = { arcCenter[0] - pointPosition[0], arcCenter[1] - pointPosition[1], arcCenter[2] - pointPosition[2] };
      vtkMath::Normalize(radialDirection);

      if (verbose)
        {
        std::cout << "Circle point [" << tupleIndex << "]" << std::endl;
        std::cout << "  normal: " << normal[0] << ", " << normal[1] << ", " << normal[2] << std::endl;
        std::cout << "  binormal: " << binormal[0] << ", " << binormal[1] << ", " << binormal[2] << std::endl;
        std::cout << "  radial: " << radialDirection[0] << ", " << radialDirection[1] << ", " << radialDirection[2] << std::endl;
        }

      // Use a bit less tight tolerance for the normal (radial) comparison,
      // because the arc is estimated from straight line segments
      CHECK_DOUBLE_TOLERANCE(normal[0], radialDirection[0], tolerance * 10.0);
      CHECK_DOUBLE_TOLERANCE(normal[1], radialDirection[1], tolerance * 10.0);
      CHECK_DOUBLE_TOLERANCE(normal[2], radialDirection[2], tolerance * 10.0);

      CHECK_DOUBLE_TOLERANCE(binormal[0], constantBinormal[0], tolerance);
      CHECK_DOUBLE_TOLERANCE(binormal[1], constantBinormal[1], tolerance);
      CHECK_DOUBLE_TOLERANCE(binormal[2], constantBinormal[2], tolerance);
      }
  }

  std::cout << "Test succeeded." << std::endl;
  return EXIT_SUCCESS;
}
