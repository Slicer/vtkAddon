/*==============================================================================

  Copyright (c) Laboratory for Percutaneous Surgery (PerkLab)
  Queen's University, Kingston, ON, Canada. All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Kyle Sunderland, PerkLab, Queen's University
  and was supported through CANARIE's Research Software Program, Cancer
  Care Ontario, OpenAnatomy, and Brigham and Women's Hospital through NIH grant R01MH112748.

==============================================================================*/
/**
 * @class   vtkImageMathematicsAddon
 * @brief   implicit function consisting of boolean combinations of implicit functions, with invert option
 *
 * This class is a subclass of vtkImageMathematics that adds an option to multiply an image by the value of
 * a second image, normalized using the specified range.
 */

#ifndef vtkImageMathematicsAddon_h
#define vtkImageMathematicsAddon_h

#define VTK_NORMALIZE_MULTIPLY 100

// VTK includes
#include "vtkImageMathematics.h"

// vtkAddon includes
#include "vtkAddon.h"

class VTK_ADDON_EXPORT vtkImageMathematicsAddon : public vtkImageMathematics
{
public:
  vtkTypeMacro(vtkImageMathematicsAddon, vtkImageMathematics);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkImageMathematicsAddon* New();

  void SetOperationToNormalizeMultiply() { this->SetOperation(VTK_NORMALIZE_MULTIPLY); }

  vtkSetVector2Macro(NormalizeRange, double);
  vtkGetVector2Macro(NormalizeRange, double);

protected:
  vtkImageMathematicsAddon();
  ~vtkImageMathematicsAddon() override;

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  void ThreadedRequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector, vtkImageData*** inData, vtkImageData** outData,
    int outExt[6], int threadId) override;

  double NormalizeRange[2] = { 0.0, 1.0 };

private:
  vtkImageMathematicsAddon(const vtkImageMathematicsAddon&) = delete;
  void operator=(const vtkImageMathematicsAddon&) = delete;
};
#endif
