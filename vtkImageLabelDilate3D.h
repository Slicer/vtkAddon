// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @class   vtkImageLabelDilate3D
 * @brief   Label image dilation filter
 *
 * vtkImageLabelDilate3D dilates a labelmap image by replacing each background voxel with the
 * most dominant label voxel in its neighborhood.
 *
 * @par Acknowledgments:
 * This class was developed by Andras Lasso PerkLab, Queen's University
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
   * Set/Get the background voxel value that label values are dilated into.
   * The filter will only change voxels that are of this value.
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
