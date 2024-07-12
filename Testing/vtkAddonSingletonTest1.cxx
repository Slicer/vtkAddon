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

#include "vtkAddonTestingMacros.h"
#include "vtkMySingletonClass.h"

//----------------------------------------------------------------------------
int TestSingletonGetValue(int value)
{
  vtkMySingletonClass* singleton = vtkMySingletonClass::GetInstance();
  CHECK_NOT_NULL(singleton);
  CHECK_INT(singleton->GetValue(), value);
  return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------
int TestSingletonSetValue(int value)
{
  vtkMySingletonClass* singleton = vtkMySingletonClass::GetInstance();
  CHECK_NOT_NULL(singleton);
  singleton->SetValue(value);
  CHECK_INT(singleton->GetValue(), value);
  return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------
int vtkAddonSingletonTest1(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  int value = 0;
  CHECK_EXIT_SUCCESS(TestSingletonGetValue(0));

  value = 1;
  CHECK_EXIT_SUCCESS(TestSingletonSetValue(value));
  CHECK_EXIT_SUCCESS(TestSingletonGetValue(value));

  value = -512;
  CHECK_EXIT_SUCCESS(TestSingletonSetValue(value));
  CHECK_EXIT_SUCCESS(TestSingletonGetValue(value));

  return EXIT_SUCCESS;
}
