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
 * @class   vtkCapPolyData
 * @brief   implicit function consisting of boolean combinations of implicit functions, with invert option
 *
 * This class is a subclass of vtkPolyDataAlgorithm that will generate an end cap for a polydata cut with the specified function
 *
 */

#ifndef vtkCapPolyData_h
#define vtkCapPolyData_h

 // VTK includes
#include <vtkAppendPolyData.h>
#include <vtkImplicitFunction.h>
#include <vtkPolyDataAlgorithm.h>

// vtkAddon includes
#include "vtkAddon.h"

class vtkPlaneCollection;

class VTK_ADDON_EXPORT vtkCapPolyData : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkCapPolyData, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkCapPolyData* New();

  ///@{
  /**
   * Specify the implicit function with which to perform the
   * clipping. If you do not define an implicit function, then the input
   * scalar data will be used for clipping.
   */
  vtkSetMacro(ClipFunction, vtkSmartPointer<vtkImplicitFunction>);
  vtkGetMacro(ClipFunction, vtkSmartPointer<vtkImplicitFunction>);
  ///@}

  /// Return the mtime also considering the locator and clip function.
  vtkMTimeType GetMTime() override;

  /// Get the list of planes from the implicit function
  static void GetPlanes(vtkImplicitFunction* function, vtkPlaneCollection* planes, vtkAbstractTransform* parentTransform = nullptr);

  ///@{
  /**
   * Set/Get whether to generate an outline of the cap.
   * Default is on.
   */
  vtkSetMacro(GenerateOutline, bool);
  vtkGetMacro(GenerateOutline, bool);
  ///@}

  ///@{
  /**
   * Set/Get whether to generate cell type scalars.
   * Default is on.
   */
  vtkSetMacro(GenerateCellTypeScalars, bool);
  vtkGetMacro(GenerateCellTypeScalars, bool);
  ///@}

protected:
  vtkCapPolyData();
  ~vtkCapPolyData() override;

  /// Generate the end cap for the input polydata cut using planes in the cutFunction.
  void CreateEndCap(vtkPlaneCollection* planes, vtkPolyData* originalPolyData,
    vtkImplicitFunction* cutFunction, vtkPolyData* outputEndCap);

  /// Updates the polydata cell scalar array to reflect the cell type.
  void UpdateCellTypeArray(vtkPolyData* polyData);

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

protected:

  vtkSmartPointer<vtkImplicitFunction> ClipFunction{ nullptr };

  bool GenerateOutline{ true };
  bool GenerateCellTypeScalars{ true };

private:
  vtkCapPolyData(const vtkCapPolyData&) = delete;
  void operator=(const vtkCapPolyData&) = delete;
};
#endif
