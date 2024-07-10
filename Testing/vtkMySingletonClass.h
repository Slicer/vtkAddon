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

// VTK includes
#include <vtkObject.h>

// vtkAddon includes
#include <vtkSingleton.h>

/// \brief Singleton VTK object containing a single integer.
class vtkMySingletonClass : public vtkObject
{
public:
  static vtkMySingletonClass* New();
  vtkTypeMacro(vtkMySingletonClass, vtkObject);

  vtkSetMacro(Value, int);
  vtkGetMacro(Value, int);

  VTK_SINGLETON_DECLARE_GETINSTANCE(vtkMySingletonClass);

protected:
  vtkMySingletonClass() = default;
  ~vtkMySingletonClass() override = default;

  int Value{ 0 };

  VTK_SINGLETON_DECLARE(vtkMySingletonClass);
};

VTK_SINGLETON_DECLARE_INITIALIZER(, vtkMySingletonClass);
