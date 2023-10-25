/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageLabelDilate3D.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkImageLabelDilate3D
 * @brief   Median Filter
 *
 * vtkImageLabelDilate3D dilates a labelmap image by replacing each background voxel with the
 * most dominant label voxel in its neighborhood.
 */

#ifndef vtkImageLabelDilate3D_h
#define vtkImageLabelDilate3D_h

#include "vtkAddonExport.h" // For export macro
#include "vtkImageSpatialAlgorithm.h"

class VTK_ADDON_EXPORT vtkImageLabelDilate3D : public vtkImageSpatialAlgorithm
{
public:
  static vtkImageLabelDilate3D* New();
  vtkTypeMacro(vtkImageLabelDilate3D, vtkImageSpatialAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Size of the neighborhood where the most dominant label value will be searched for.
   */
  void SetKernelSize(int size0, int size1, int size2);

  ///@{
  /**
   * Set/Get the background voxel value that is excluded from the median computation.
   */
  vtkSetMacro(BackgroundValue, double);
  vtkGetMacro(BackgroundValue, double);
  ///@}

protected:
  vtkImageLabelDilate3D();
  ~vtkImageLabelDilate3D() override;

  double BackgroundValue;

  void ThreadedRequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector, vtkImageData*** inData, vtkImageData** outData,
    int outExt[6], int id) override;

private:
  vtkImageLabelDilate3D(const vtkImageLabelDilate3D&) = delete;
  void operator=(const vtkImageLabelDilate3D&) = delete;
};

#endif
