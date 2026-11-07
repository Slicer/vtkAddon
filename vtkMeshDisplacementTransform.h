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

#ifndef vtkMeshDisplacementTransform_h
#define vtkMeshDisplacementTransform_h

// vtkAddon includes
#include "vtkAddon.h"

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkWarpTransform.h>

// STD includes
#include <atomic>

class vtkDataArray;
class vtkIdTypeArray;
class vtkPointSet;
class vtkPolyData;
class vtkStaticCellLocator;
class vtkStaticPointLocator;

/// \brief Transform that moves each position by the displacement (3-component point data) of a mesh.
///
/// It is a transform that is created from displacements at the positions of the mesh points, but it defines the
/// displacement at any position in the cells of the mesh and, with a border falloff, slightly outside them. So it can
/// deform other data (images when resampled, other meshes, points), it can be inverted and concatenated with other
/// transforms, and it can be sampled into a displacement field. A typical use is a displacement computed on a
/// volumetric mesh, such as the result of a finite element simulation, applied to images or other models.
///
/// How it moves a position: a mesh point moves by its own vector; any other position by the vector interpolated in the
/// cell that contains it (as vtkProbeFilter does); outside the cells (for example around a surface mesh or a point
/// cloud) positions do not move, unless BorderFalloffDistance is set. The inverse is iterative (Newton's method),
/// reliable where the field does not fold.
///
/// Compared to vtkWarpVector: that filter only outputs a warped copy of the mesh that holds the vectors. Transforming
/// the mesh itself with this transform gives the same result (its points move exactly by their own vectors, without a
/// cell search), but as a transform it also applies to anything else, as described above.
///
/// Compared to vtkScatteredBSplineTransform, which takes the same input: that one fits a smooth B-spline to the
/// displacements at the mesh points, ignoring the cells. Its field is smooth everywhere (continuous derivatives, also
/// across cell boundaries and outside the mesh) and is a B-spline transform that can be saved as one, but it only
/// approximates the displacements, up to its fitting tolerance, and the fit takes time and memory that grow with the
/// detail it has to follow. This transform uses the displacements as they are, interpolated in the cells: exact at the
/// mesh points and as accurate as the mesh in between, with no fitting, but its derivative is only piecewise smooth
/// (it changes across cell boundaries), and outside the cells it is defined by the border falloff. Use this one when the
/// displacements must be followed exactly (warping the mesh itself, or data inside a volumetric mesh), the B-spline fit
/// when a smooth field matters more (a point cloud, or a surface whose inside should move with it).
///
/// Other ways of making a transform from such displacements: vtkThinPlateSplineTransform also fits a smooth field, but
/// is slow for many points; probing the mesh onto an image (vtkProbeFilter) for a vtkGridTransform loses accuracy
/// depending on the grid spacing.
///
/// Border treatment. Outside the cells, within BorderFalloffDistance, a position moves by the displacement of the
/// closest point of the mesh, weighted by a cubic function of the distance that has zero slope at both ends. The
/// closest point is found on the boundary surface of the mesh (its faces), where it is exact and interpolated as at the
/// border of the cell; the closest point that VTK finds on a 3D cell from outside is an estimate. The remaining
/// limitations: the derivative is discontinuous at the border of the mesh, and around a non-convex mesh the
/// displacement jumps where two parts of the mesh are equally close. The jump is smaller, and stays closer to the mesh,
/// with a smaller BorderFalloffDistance, so the best choice is the smallest distance that keeps the transform
/// invertible, which AutomaticBorderFalloffDistance provides.
class VTK_ADDON_EXPORT vtkMeshDisplacementTransform : public vtkWarpTransform
{
public:
  static vtkMeshDisplacementTransform* New();
  vtkTypeMacro(vtkMeshDisplacementTransform, vtkWarpTransform);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /// Make another transform of the same type.
  vtkAbstractTransform* MakeTransform() override;

  /// The modification time includes that of the input mesh, so that editing the mesh is noticed.
  vtkMTimeType GetMTime() override;

  ///@{
  /// The mesh whose points carry the displacements. The mesh is not modified.
  virtual void SetInputPointSet(vtkPointSet* pointSet);
  vtkGetObjectMacro(InputPointSet, vtkPointSet);
  ///@}

  ///@{
  /// Name of the 3-component point data array that holds the displacements. If it is not set,
  /// the active vectors are used.
  vtkSetStringMacro(DisplacementArrayName);
  vtkGetStringMacro(DisplacementArrayName);
  ///@}

  ///@{
  /// Multiplies the displacements. Default is 1.
  vtkSetMacro(DisplacementScale, double);
  vtkGetMacro(DisplacementScale, double);
  ///@}

  ///@{
  /// Distance from the cells of the mesh over which the displacement falls smoothly to zero. Default is 0: positions outside
  /// the cells do not move, so the displacement jumps to zero at the boundary of the mesh, which can make the iterative
  /// inverse fail near the boundary. If it is positive, a position outside the cells gets the displacement at the closest
  /// point of the mesh, multiplied by a weight that falls from 1 at the mesh to 0 at this distance (with zero slope at both
  /// ends). The transform remains invertible if the distance is larger than about 1.5 times the largest displacement.
  /// Set it to VTK_DOUBLE_MAX to keep the displacement of the closest mesh point everywhere outside the mesh (as
  /// vtkBSplineTransform does in Edge border mode beyond its grid).
  /// Positions inside a closed surface mesh are outside its cells, so they get the falloff too. Not used for a point cloud
  /// without cells, where the displacement of the closest point would jump between points.
  vtkSetMacro(BorderFalloffDistance, double);
  vtkGetMacro(BorderFalloffDistance, double);
  ///@}

  ///@{
  /// If on, BorderFalloffDistance is ignored and the border falloff distance is computed from the displacements as twice
  /// the largest displacement (|DisplacementScale| times the longest vector of the array), which keeps the transform
  /// invertible and follows the mesh, the array and the scale when they change. Default is off.
  /// The transform folds, and cannot be inverted, if the distance is less than 1.5 times the largest displacement (1.5 is
  /// the largest slope of the falloff weight); twice leaves a margin, so the falloff compresses space by at most a factor
  /// of 4 and the inverse converges reliably (vtkMeshDisplacementTransformTest3). A larger distance makes the falloff
  /// smoother (the gradient scales with its inverse) but moves positions farther from the mesh and makes the jump around
  /// a non-convex mesh larger.
  vtkSetMacro(AutomaticBorderFalloffDistance, bool);
  vtkGetMacro(AutomaticBorderFalloffDistance, bool);
  vtkBooleanMacro(AutomaticBorderFalloffDistance, bool);
  ///@}

  /// The border falloff distance in use: computed from the displacements if AutomaticBorderFalloffDistance is on,
  /// BorderFalloffDistance otherwise (0 if there is no usable input). Update() must have been called.
  double GetEffectiveBorderFalloffDistance();

  /// Displacement at a mesh point, interpolated in the containing cell, or within the border falloff distance of the
  /// cells, times DisplacementScale. Update() must have been called. Returns false, and a zero displacement, if the
  /// position is not in any of these.
  bool GetDisplacement(const double position[3], double displacement[3]);

  ///@{
  /// The inverse is computed iteratively at each position. If it does not converge (where the transform folds, for
  /// example around a non-convex mesh), the closest estimate is used and the failure is counted. Only the first failure
  /// after the transform is modified is logged as a warning, as in vtkOrientedGridTransform, unless this is on, which
  /// logs every one. Default is off.
  vtkSetMacro(LogAllInverseConvergenceFailures, bool);
  vtkGetMacro(LogAllInverseConvergenceFailures, bool);
  vtkBooleanMacro(LogAllInverseConvergenceFailures, bool);
  ///@}

  /// Number of positions where the inverse did not converge since this transform object was last updated. For the
  /// inverse of a transform, ask the object that GetInverse() returns.
  vtkIdType GetNumberOfInverseConvergenceFailures() { return this->NumberOfInverseConvergenceFailures.load(); }

protected:
  vtkMeshDisplacementTransform();
  ~vtkMeshDisplacementTransform() override;

  void InternalUpdate() override;
  void InternalDeepCopy(vtkAbstractTransform* transform) override;

  void ForwardTransformPoint(const float in[3], float out[3]) override;
  void ForwardTransformPoint(const double in[3], double out[3]) override;
  void ForwardTransformDerivative(const float in[3], float out[3], float derivative[3][3]) override;
  void ForwardTransformDerivative(const double in[3], double out[3], double derivative[3][3]) override;

  /// Newton's method with a Jacobian that is only updated when the error does not decrease fast enough, steps that are
  /// shortened instead of increasing the error, and an early stop once a fresh Jacobian cannot decrease it any more.
  /// This needs a fraction of the evaluations of vtkWarpTransform's inverse, which costs most in the falloff zone.
  void InverseTransformPoint(const float in[3], float out[3]) override;
  void InverseTransformPoint(const double in[3], double out[3]) override;
  void InverseTransformDerivative(const float in[3], float out[3], float derivative[3][3]) override;
  void InverseTransformDerivative(const double in[3], double out[3], double derivative[3][3]) override;

  /// Displacement at the closest point of the boundary surface of the mesh, within the given distance, weighted by the
  /// falloff. Returns false if no boundary is that close.
  bool GetBorderFalloffDisplacement(const double position[3], double borderFalloffDistance, double displacement[3]);

  /// Count, and log as set by LogAllInverseConvergenceFailures, a position where the inverse did not converge
  void ReportInverseConvergenceFailure(const double position[3], double error, int numberOfIterations);

  vtkPointSet* InputPointSet{ nullptr };
  char* DisplacementArrayName{ nullptr };
  double DisplacementScale{ 1.0 };
  double BorderFalloffDistance{ 0.0 };
  bool AutomaticBorderFalloffDistance{ false };
  bool LogAllInverseConvergenceFailures{ false };

  /// Inverse failures since the last update, and the time of the last one that was logged; the inverse may be
  /// computed from several threads at once
  std::atomic<vtkIdType> NumberOfInverseConvergenceFailures{ 0 };
  std::atomic<vtkMTimeType> LastInverseConvergenceWarningTime{ 0 };

  // Set up by InternalUpdate()
  vtkSmartPointer<vtkDataArray> Displacements;
  /// BorderFalloffDistance, or twice the largest displacement if AutomaticBorderFalloffDistance is on
  double EffectiveBorderFalloffDistance{ 0.0 };
  vtkSmartPointer<vtkStaticPointLocator> PointLocator;
  vtkSmartPointer<vtkStaticCellLocator> CellLocator;
  /// Boundary surface of the mesh (its faces, or the mesh itself if it is a vtkPolyData), for the closest point of the
  /// mesh in the falloff zone: on a face it is found exactly and its interpolation matches that of the cell at the border
  vtkSmartPointer<vtkPolyData> BoundarySurface;
  vtkSmartPointer<vtkStaticCellLocator> BoundaryLocator;
  /// Point of the mesh of each boundary surface point (nullptr if the boundary surface is the mesh itself)
  vtkSmartPointer<vtkIdTypeArray> BoundaryPointIds;
  /// Largest number of points in a face of the boundary surface
  int MaximumBoundaryCellSize{ 0 };
  /// Squared distance within which a position is the same as a mesh point
  double PointTolerance2{ 0.0 };
  /// Squared distance within which a position is inside a cell
  double CellTolerance2{ 0.0 };
  /// Step of the finite differences that the derivative is computed with
  double DerivativeStep{ 1.0 };
  /// Largest number of points in a cell of the mesh, for the interpolation weights
  int MaximumCellSize{ 0 };
  /// Bounds of the mesh, for quickly skipping positions far from it
  double MeshBounds[6]{ 0.0, -1.0, 0.0, -1.0, 0.0, -1.0 };

private:
  vtkMeshDisplacementTransform(const vtkMeshDisplacementTransform&) = delete;
  void operator=(const vtkMeshDisplacementTransform&) = delete;
};

#endif
