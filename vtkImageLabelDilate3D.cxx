// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause

#include "vtkImageLabelDilate3D.h"

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <map>
#include <vector>

vtkStandardNewMacro(vtkImageLabelDilate3D);

//------------------------------------------------------------------------------
// Construct an instance of vtkImageLabelDilate3D filter.
vtkImageLabelDilate3D::vtkImageLabelDilate3D()
{
  this->BackgroundValue = 0;
  this->SetKernelSize(1, 1, 1);
  this->HandleBoundaries = 1;
}

//------------------------------------------------------------------------------
vtkImageLabelDilate3D::~vtkImageLabelDilate3D() = default;

//------------------------------------------------------------------------------
void vtkImageLabelDilate3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "BackgroundValue: " << this->BackgroundValue << endl;
}

//------------------------------------------------------------------------------
// This method sets the size of the neighborhood.  It also sets the
// default middle of the neighborhood
void vtkImageLabelDilate3D::SetKernelSize(int size0, int size1, int size2)
{
  if (this->KernelSize[0] == size0 && this->KernelSize[1] == size1 && this->KernelSize[2] == size2)
    {
    // no change
    return;
    }

  // Set the kernel size and middle
  this->KernelSize[0] = size0;
  this->KernelMiddle[0] = size0 / 2;
  this->KernelSize[1] = size1;
  this->KernelMiddle[1] = size1 / 2;
  this->KernelSize[2] = size2;
  this->KernelMiddle[2] = size2 / 2;

  this->Modified();
}

namespace
{

//------------------------------------------------------------------------------
// Compute the most frequently occurring number in the vector
template <typename T>
T vtkComputeModeOfArray(const std::vector<T>& values) {
  // A map to store the count of each value
  std::map<T, int> valueCount;
  T maxValue = 0; // the most frequent value found so far
  std::size_t maxCount = 0; // how many of this value are found so far
  typename std::vector<T>::size_type countForMajority = typename std::vector<T>::size_type(values.size() / 2); // more than this count means majority
  for (T value : values)
    {
    std::size_t count = ++valueCount[value];
    if (count > maxCount)
      {
      maxValue = value;
      maxCount = count;
      if (count > countForMajority)
        {
        // already more than half of the values are this, so we don't need to continue the search
        break;
        }
      }
    }
  return maxValue;
}

} // end anonymous namespace

//------------------------------------------------------------------------------
// This method contains the second switch statement that calls the correct
// templated function for the mask types.
template <typename T>
void vtkImageLabelDilate3DExecute(vtkImageLabelDilate3D* self, vtkImageData* inData, T* inPtr,
  vtkImageData* outData, T* outPtr, int outExt[6], int id, vtkDataArray* inArray)
{
  if (!inArray)
    {
    return;
    }

  // The logic of iterating through a box neighborhood was adopted from vtkImageMedian3D filter.

  // Get information to march through data
  vtkIdType inInc0, inInc1, inInc2;
  inData->GetIncrements(inInc0, inInc1, inInc2);
  vtkIdType outIncX, outIncY, outIncZ;
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  int* kernelMiddle = self->GetKernelMiddle();
  int* kernelSize = self->GetKernelSize();

  int numComp = inArray->GetNumberOfComponents();

  // For looping through hood pixels
  int hoodMin0 = outExt[0] - kernelMiddle[0];
  int hoodMin1 = outExt[2] - kernelMiddle[1];
  int hoodMin2 = outExt[4] - kernelMiddle[2];
  int hoodMax0 = kernelSize[0] + hoodMin0 - 1;
  int hoodMax1 = kernelSize[1] + hoodMin1 - 1;
  int hoodMax2 = kernelSize[2] + hoodMin2 - 1;

  // Clip by the input image extent
  int* inExt = inData->GetExtent();
  hoodMin0 = (hoodMin0 > inExt[0]) ? hoodMin0 : inExt[0];
  hoodMin1 = (hoodMin1 > inExt[2]) ? hoodMin1 : inExt[2];
  hoodMin2 = (hoodMin2 > inExt[4]) ? hoodMin2 : inExt[4];
  hoodMax0 = (hoodMax0 < inExt[1]) ? hoodMax0 : inExt[1];
  hoodMax1 = (hoodMax1 < inExt[3]) ? hoodMax1 : inExt[3];
  hoodMax2 = (hoodMax2 < inExt[5]) ? hoodMax2 : inExt[5];

  // Save the starting neighborhood dimensions (2 loops only once)
  int hoodStartMin0 = hoodMin0;
  int hoodStartMax0 = hoodMax0;
  int hoodStartMin1 = hoodMin1;
  int hoodStartMax1 = hoodMax1;

  // The portion of the output that needs no boundary computation.
  int middleMin0 = inExt[0] + kernelMiddle[0];
  int middleMax0 = inExt[1] - (kernelSize[0] - 1) + kernelMiddle[0];
  int middleMin1 = inExt[2] + kernelMiddle[1];
  int middleMax1 = inExt[3] - (kernelSize[1] - 1) + kernelMiddle[1];
  int middleMin2 = inExt[4] + kernelMiddle[2];
  int middleMax2 = inExt[5] - (kernelSize[2] - 1) + kernelMiddle[2];

  unsigned long target = static_cast<unsigned long>((outExt[5] - outExt[4] + 1) * (outExt[3] - outExt[2] + 1) / 50.0) + 1;
  unsigned long count = 0;

  T backgroundValue = static_cast<T>(self->GetBackgroundValue());

  // loop through pixel of output
  inPtr = static_cast<T*>(inArray->GetVoidPointer((hoodMin0 - inExt[0]) * inInc0 +
    (hoodMin1 - inExt[2]) * inInc1 + (hoodMin2 - inExt[4]) * inInc2));
  T* inPtr2 = inPtr;
  for (int outIdx2 = outExt[4]; outIdx2 <= outExt[5]; ++outIdx2)
    {
    T* inPtr1 = inPtr2;
    hoodMin1 = hoodStartMin1;
    hoodMax1 = hoodStartMax1;
    vtkIdType hoodCenterOffset2 = (hoodMax2 - hoodMin2 ) / 2 * inInc2;
    for (int outIdx1 = outExt[2]; !self->AbortExecute && outIdx1 <= outExt[3]; ++outIdx1)
      {
      if (!id)
        {
        if (!(count % target))
          {
          self->UpdateProgress(count / (50.0 * target));
          }
        count++;
        }
      T* inPtr0 = inPtr1;
      hoodMin0 = hoodStartMin0;
      hoodMax0 = hoodStartMax0;
      vtkIdType hoodCenterOffset1 = (hoodMax1 - hoodMin1) / 2 * inInc1;
      T *tmpPtr0, *tmpPtr1, *tmpPtr2;
      std::vector<T> workArray;
      for (int outIdx0 = outExt[0]; outIdx0 <= outExt[1]; ++outIdx0)
        {
        vtkIdType hoodCenterOffset = hoodCenterOffset2 + hoodCenterOffset1 + (hoodMax0 - hoodMin0) / 2 * inInc0;

        for (int outIdxC = 0; outIdxC < numComp; outIdxC++)
          {
          T centerVoxelValue = *(inPtr0 + hoodCenterOffset + outIdxC);
          if (centerVoxelValue != backgroundValue)
            {
            // Center voxel is not background voxel, leave it unchanged
            *outPtr++ = centerVoxelValue;
            continue;
            }

          // Center voxel is a background voxel, replace it with the most frequent non-background value
          // in the neighborhood.

          workArray.clear();
          tmpPtr2 = inPtr0 + outIdxC;
          for (int hoodIdx2 = hoodMin2; hoodIdx2 <= hoodMax2; ++hoodIdx2)
            {
            tmpPtr1 = tmpPtr2;
            for (int hoodIdx1 = hoodMin1; hoodIdx1 <= hoodMax1; ++hoodIdx1)
              {
              tmpPtr0 = tmpPtr1;
              for (int hoodIdx0 = hoodMin0; hoodIdx0 <= hoodMax0; ++hoodIdx0)
                {
                if (*tmpPtr0 != backgroundValue)
                  {
                  workArray.push_back(*tmpPtr0);
                  }
                tmpPtr0 += inInc0;
                }
              tmpPtr1 += inInc1;
              }
            tmpPtr2 += inInc2;
            }

          // Replace this pixel with the mode of the neighborhood
          if (workArray.empty())
            {
            *outPtr++ = backgroundValue;
            }
          else
            {
            *outPtr++ = vtkComputeModeOfArray(workArray);
            }
        }

        // shift neighborhood considering boundaries
        if (outIdx0 >= middleMin0)
          {
          inPtr0 += inInc0;
          ++hoodMin0;
          }
        if (outIdx0 < middleMax0)
          {
          ++hoodMax0;
          }
      }
      // shift neighborhood considering boundaries
      if (outIdx1 >= middleMin1)
        {
        inPtr1 += inInc1;
        ++hoodMin1;
        }
      if (outIdx1 < middleMax1)
        {
        ++hoodMax1;
        }
      outPtr += outIncY;
      }
    // shift neighborhood considering boundaries
    if (outIdx2 >= middleMin2)
      {
      inPtr2 += inInc2;
      ++hoodMin2;
      }
    if (outIdx2 < middleMax2)
      {
      ++hoodMax2;
      }
    outPtr += outIncZ;
    }
}

//------------------------------------------------------------------------------
// This method contains the first switch statement that calls the correct
// templated function for the input and output region types.
void vtkImageLabelDilate3D::ThreadedRequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* vtkNotUsed(outputVector),
  vtkImageData*** inData, vtkImageData** outData, int outExt[6], int id)
{
  void* inPtr;
  void* outPtr = outData[0]->GetScalarPointerForExtent(outExt);

  vtkDataArray* inArray = this->GetInputArrayToProcess(0, inputVector);
  if (id == 0)
    {
    outData[0]->GetPointData()->GetScalars()->SetName(inArray->GetName());
    }

  inPtr = inArray->GetVoidPointer(0);

  // this filter expects that input is the same type as output.
  if (inArray->GetDataType() != outData[0]->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: input data type, " << inArray->GetDataType()
                  << ", must match out ScalarType " << outData[0]->GetScalarType());
    return;
    }

  switch (inArray->GetDataType())
    {
    vtkTemplateMacro(vtkImageLabelDilate3DExecute(this, inData[0][0], static_cast<VTK_TT*>(inPtr),
      outData[0], static_cast<VTK_TT*>(outPtr), outExt, id, inArray));
    default:
      vtkErrorMacro(<< "Execute: Unknown input ScalarType");
      return;
    }
}
