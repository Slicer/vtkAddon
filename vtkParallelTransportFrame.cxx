/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

=========================================================================auto=*/

/*
  Portions of vtkParallelTransportFrame::ComputeAxisDirections method
  are covered under the VMTK copyright and BSD license:

    Copyright (c) Luca Antiga, David Steinman. All rights reserved.
    See https://github.com/vmtk/vmtk/blob/master/LICENSE file for details.
*/

#include <algorithm> // VTK 8.2.0 has bug for C++17 "clamp" function (algorithm must be included before vtMath.h)

#include "vtkParallelTransportFrame.h"

#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyLine.h"

vtkStandardNewMacro(vtkParallelTransportFrame);

//----------------------------------------------------------------------------
vtkParallelTransportFrame::vtkParallelTransportFrame()
{
  this->SetTangentsArrayName("Tangents");
  this->SetNormalsArrayName("Normals");
  this->SetBinormalsArrayName("Binormals");
}

//----------------------------------------------------------------------------
vtkParallelTransportFrame::~vtkParallelTransportFrame()
{
  this->SetTangentsArrayName(nullptr);
  this->SetNormalsArrayName(nullptr);
  this->SetBinormalsArrayName(nullptr);
}

//----------------------------------------------------------------------------
void vtkParallelTransportFrame::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "TangentsArrayName: " << (this->TangentsArrayName ? this->TangentsArrayName : "(none)") << "\n";
  os << indent << "NormalsArrayName: " << (this->NormalsArrayName ? this->NormalsArrayName : "(none)") << "\n";
  os << indent << "BinormalsArrayName: " << (this->BinormalsArrayName ? this->BinormalsArrayName : "(none)") << "\n";

  os << indent << "Tolerance: " << this->Tolerance << "\n";
  os << indent << "MinimumDistance: " << this->MinimumDistance << "\n";

  os << indent << "PreferredInitialNormalVector: ("
    << this->PreferredInitialNormalVector[0] << ", " << this->PreferredInitialNormalVector[1] << ", " << this->PreferredInitialNormalVector[2] << ")\n";
  os << indent << "PreferredInitialBinormalVector: ("
    << this->PreferredInitialBinormalVector[0] << ", " << this->PreferredInitialBinormalVector[1] << ", " << this->PreferredInitialBinormalVector[2] << ")\n";
}

//----------------------------------------------------------------------------
int vtkParallelTransportFrame::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the input and ouptut
  vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!this->TangentsArrayName)
    {
    vtkErrorMacro(<< "TangentsArrayName is not specified");
    return 0;
    }
  if (!this->NormalsArrayName)
    {
    vtkErrorMacro(<< "NormalsArrayName is not specified");
    return 0;
    }
  if (!this->BinormalsArrayName)
    {
    vtkErrorMacro(<< "BinormalsArrayName is not specified");
    return 0;
    }

  output->DeepCopy(input);

  vtkNew<vtkDoubleArray> tangentsArray;
  tangentsArray->SetNumberOfComponents(3);
  tangentsArray->SetName(this->TangentsArrayName);

  vtkNew<vtkDoubleArray> normalsArray;
  normalsArray->SetNumberOfComponents(3);
  normalsArray->SetName(this->NormalsArrayName);

  vtkNew<vtkDoubleArray> binormalsArray;
  binormalsArray->SetNumberOfComponents(3);
  binormalsArray->SetName(this->BinormalsArrayName);

  vtkIdType numberOfPoints = input->GetNumberOfPoints();
  tangentsArray->SetNumberOfTuples(numberOfPoints);
  tangentsArray->Fill(0.0);
  normalsArray->SetNumberOfTuples(numberOfPoints);
  normalsArray->Fill(0.0);
  binormalsArray->SetNumberOfTuples(numberOfPoints);
  binormalsArray->Fill(0.0);

  int numberOfCells = input->GetNumberOfCells();
  for (int cellIndex = 0; cellIndex < numberOfCells; cellIndex++)
    {
    this->ComputeAxisDirections2(input, cellIndex, tangentsArray, normalsArray, binormalsArray);
    }

  output->GetPointData()->AddArray(tangentsArray);
  output->GetPointData()->AddArray(normalsArray);
  output->GetPointData()->AddArray(binormalsArray);

  return 1;
}

//----------------------------------------------------------------------------
void vtkParallelTransportFrame::ComputeAxisDirections(vtkPolyData* input, vtkIdType cellIndex, vtkDoubleArray* tangentsArray, vtkDoubleArray* normalsArray, vtkDoubleArray* binormalsArray)
{
  vtkPolyLine* polyLine = vtkPolyLine::SafeDownCast(input->GetCell(cellIndex));
  if (!polyLine)
    {
    return;
    }
  vtkIdType numberOfPointsInCell = polyLine->GetNumberOfPoints();
  if (numberOfPointsInCell < 2)
    {
    return;
    }

  double tangent0[3] = { 0.0, 0.0, 0.0 };
  vtkIdType pointId0 = polyLine->GetPointId(0);
  double pointPosition0[3];
  input->GetPoint(pointId0, pointPosition0);

  // Find tangent by direction vector by moving a minimal distance from the initial point
  for (int pointIndex = 1; pointIndex < numberOfPointsInCell; pointIndex++)
    {
    vtkIdType pointId1 = polyLine->GetPointId(pointIndex);
    double pointPosition1[3];
    input->GetPoint(pointId1, pointPosition1);
    tangent0[0] = pointPosition1[0] - pointPosition0[0];
    tangent0[1] = pointPosition1[1] - pointPosition0[1];
    tangent0[2] = pointPosition1[2] - pointPosition0[2];
    if (vtkMath::Norm(tangent0) >= this->MinimumDistance)
      {
      break;
      }
    }
  vtkMath::Normalize(tangent0);

  // Compute initial normal and binormal directions from the initial tangent and preferred
  // normal/binormal directions.
  double normal0[3] = {0.0, 0.0, 0.0};
  double binormal0[3] = {0.0, 0.0, 0.0};
  vtkMath::Cross(tangent0, this->PreferredInitialNormalVector, binormal0);
  if (vtkMath::Norm(binormal0) > this->Tolerance)
    {
    vtkMath::Normalize(binormal0);
    vtkMath::Cross(binormal0, tangent0, normal0);
    }
  else
    {
    vtkMath::Cross(this->PreferredInitialBinormalVector, tangent0, normal0);
    vtkMath::Normalize(normal0);
    vtkMath::Cross(tangent0, normal0, binormal0);
    }

  tangentsArray->SetTuple(pointId0, tangent0);
  normalsArray->SetTuple(pointId0, normal0);
  binormalsArray->SetTuple(pointId0, binormal0);

  vtkIdType pointId2 = -1;
  double tangent1[3] = { tangent0[0], tangent0[1], tangent0[2] };
  double normal1[3] = { normal0[0], normal0[1], normal0[2] };
  double binormal1[3] = { binormal0[0], binormal0[1], binormal0[2] };
  for (int i = 1; i < numberOfPointsInCell - 1; i++)
    {
    vtkIdType pointId1 = polyLine->GetPointId(i);
    pointId2 = polyLine->GetPointId(i+1);
    double pointPosition1[3];
    double pointPosition2[3];
    input->GetPoint(pointId1, pointPosition1);
    input->GetPoint(pointId2, pointPosition2);

    tangent1[0] = pointPosition2[0] - pointPosition1[0];
    tangent1[1] = pointPosition2[1] - pointPosition1[1];
    tangent1[2] = pointPosition2[2] - pointPosition1[2];

    vtkMath::Normalize(tangent0);
    vtkMath::Normalize(tangent1);

    double dot = vtkMath::Dot(tangent0, tangent1);
    double theta = 0.0;
   if ((1 - dot) < this->Tolerance)
      {
      theta = 0.0;
      }
    else
      {
      theta = acos(dot);
      }

    double rotationAxis[3];
    vtkMath::Cross(tangent0, tangent1, rotationAxis);

    vtkParallelTransportFrame::RotateVector(normal0, normal1, rotationAxis, theta);

    dot = vtkMath::Dot(tangent1, normal1);
    normal1[0] -= dot * tangent1[0];
    normal1[1] -= dot * tangent1[1];
    normal1[2] -= dot * tangent1[2];

    vtkMath::Normalize(normal1);
    vtkMath::Cross(tangent1, normal1, binormal1);

    tangentsArray->SetTuple(pointId1, tangent1);
    normalsArray->SetTuple(pointId1, normal1);
    binormalsArray->SetTuple(pointId1, binormal1);

    // Save current data for next iteration
    tangent0[0] = tangent1[0];
    tangent0[1] = tangent1[1];
    tangent0[2] = tangent1[2];
    normal0[0] = normal1[0];
    normal0[1] = normal1[1];
    normal0[2] = normal1[2];
    }

  if (pointId2 >= 0)
    {
    tangentsArray->SetTuple(pointId2, tangent1);
    normalsArray->SetTuple(pointId2, normal1);
    binormalsArray->SetTuple(pointId2, binormal1);
    }
}

//----------------------------------------------------------------------------
void vtkParallelTransportFrame::ComputeAxisDirections2(vtkPolyData* input, vtkIdType cellIndex, vtkDoubleArray* tangentsArray, vtkDoubleArray* normalsArray, vtkDoubleArray* binormalsArray)
{
  vtkPolyLine* polyLine = vtkPolyLine::SafeDownCast(input->GetCell(cellIndex));
  if (!polyLine)
    {
    return;
    }
  vtkIdType numberOfPointsInCell = polyLine->GetNumberOfPoints();
  if (numberOfPointsInCell < 2)
    {
    return;
    }

  double tangent0[3] = { 0.0, 0.0, 0.0 };
  vtkIdType pointId0 = polyLine->GetPointId(0);
  double pointPosition0[3];
  input->GetPoint(pointId0, pointPosition0);

  // Find tangent by direction vector by moving a minimal distance from the initial point
  for (int pointIndex = 1; pointIndex < numberOfPointsInCell; pointIndex++)
    {
    vtkIdType pointId1 = polyLine->GetPointId(pointIndex);
    double pointPosition1[3];
    input->GetPoint(pointId1, pointPosition1);
    tangent0[0] = pointPosition1[0] - pointPosition0[0];
    tangent0[1] = pointPosition1[1] - pointPosition0[1];
    tangent0[2] = pointPosition1[2] - pointPosition0[2];
    if (vtkMath::Norm(tangent0) >= this->MinimumDistance)
      {
      break;
      }
    }
  vtkMath::Normalize(tangent0);

  // Compute initial normal and binormal directions from the initial tangent and preferred
  // normal/binormal directions.
  double normal0[3] = {0.0, 0.0, 0.0};
  double binormal0[3] = {0.0, 0.0, 0.0};
  vtkMath::Cross(tangent0, this->PreferredInitialNormalVector, binormal0);
  if (vtkMath::Norm(binormal0) > this->Tolerance)
    {
    vtkMath::Normalize(binormal0);
    vtkMath::Cross(binormal0, tangent0, normal0);
    }
  else
    {
    vtkMath::Cross(this->PreferredInitialBinormalVector, tangent0, normal0);
    vtkMath::Normalize(normal0);
    vtkMath::Cross(tangent0, normal0, binormal0);
    }

  tangentsArray->SetTuple(pointId0, tangent0);
  normalsArray->SetTuple(pointId0, normal0);
  binormalsArray->SetTuple(pointId0, binormal0);





  vtkIdType pointId2 = -1;
  double tangent1[3] = { tangent0[0], tangent0[1], tangent0[2] };
  for (int i = 1; i < numberOfPointsInCell - 1; i++)
    {
    vtkIdType pointId1 = polyLine->GetPointId(i);
    pointId2 = polyLine->GetPointId(i+1);
    double pointPosition1[3];
    double pointPosition2[3];
    input->GetPoint(pointId1, pointPosition1);
    input->GetPoint(pointId2, pointPosition2);

    tangent1[0] = pointPosition2[0] - pointPosition1[0];
    tangent1[1] = pointPosition2[1] - pointPosition1[1];
    tangent1[2] = pointPosition2[2] - pointPosition1[2];

    vtkMath::Normalize(tangent1);
    
    tangentsArray->SetTuple(pointId1, tangent1);

    // Save current data for next iteration
    tangent0[0] = tangent1[0];
    tangent0[1] = tangent1[1];
    tangent0[2] = tangent1[2];
    }

  if (pointId2 >= 0)
    tangentsArray->SetTuple(pointId2, tangent1);



  vtkIdType pointIdi = -1;
  vtkIdType pointIdip1 = -1;
  double ri[3] = {normal0[0],normal0[1],normal0[2]};
  for (int i = 0; i < numberOfPointsInCell - 1; i++)
    {
    pointIdi = polyLine->GetPointId(i);
    pointIdip1 = polyLine->GetPointId(i+1);
    // Get positions
    double xi[3];
    double xip1[3];
    input->GetPoint(pointIdi, xi);
    input->GetPoint(pointIdip1, xip1);
    // Get tangent vectors
    double ti[3];
    double tip1[3];
    tangentsArray->GetTuple(pointIdi, ti);
    tangentsArray->GetTuple(pointIdip1, tip1);

    // Compute reflection vector of R1
    double v1[3], v1_[3] = {0,0,0};
    v1[0] = v1_[0] = xip1[0] - xi[0];
    v1[1] = v1_[1] = xip1[1] - xi[1];
    v1[2] = v1_[2] = xip1[2] - xi[2];
    double c1 = 0;
    c1 = vtkMath::Dot(v1,v1);

    // Compute rLi = R1*ri
    double v1ri = 0;
    v1ri = vtkMath::Dot(v1,ri);
    vtkMath::MultiplyScalar(v1,v1ri);
    vtkMath::MultiplyScalar(v1,-2/c1);
    double rLi[3] = {0,0,0};
    vtkMath::Add(ri,v1,rLi);

    // Compute tLi = R1*ti
    double v1ti = 0;
    v1ti = vtkMath::Dot(v1_,ti);
    vtkMath::MultiplyScalar(v1_,v1ti);
    vtkMath::MultiplyScalar(v1_,-2/c1);
    double tLi[3] = {0,0,0};
    vtkMath::Add(ti,v1_,tLi);

    // Compute reflection vector of R2
    double v2[3] = {0,0,0};
    vtkMath::Subtract(tip1,tLi,v2);

    double c2 = 0;
    c2 = vtkMath::Dot(v2,v2);

    // Compute rip1 = R2*rLi (normal vector)
    double v2rLi = 0;
    v2rLi = vtkMath::Dot(v2,rLi);
    vtkMath::MultiplyScalar(v2,v2rLi);
    vtkMath::MultiplyScalar(v2,-2/c2);
    double rip1[3] = {0,0,0};
    vtkMath::Add(rLi,v2,rip1);
    vtkMath::Normalize(rip1);

    // Compute vector sip1 of the frame i+1 (binormal vector)
    double sip1[3] = {0,0,0};
    vtkMath::Cross(tip1,rip1,sip1);
    vtkMath::Normalize(sip1);

    //tangentsArray->SetTuple(pointIdip1, tip1);
    normalsArray->SetTuple(pointIdip1, rip1);
    binormalsArray->SetTuple(pointIdip1, sip1);
    }
}

//----------------------------------------------------------------------------
void vtkParallelTransportFrame::RotateVector(double* inVector, double* outVector, const double* axis, double angle)
{
  double UdotN = vtkMath::Dot(inVector, axis);
  double NcrossU[3];
  vtkMath::Cross(axis, inVector, NcrossU);
  for (int comp = 0; comp < 3; comp++)
    {
    outVector[comp] = cos(angle) * inVector[comp]
      + (1 - cos(angle)) * UdotN * axis[comp]
      + sin(angle) * NcrossU[comp];
    }
}
