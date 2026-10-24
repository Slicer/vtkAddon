/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

=========================================================================auto=*/

/// \brief vtkScatteredBSplineTransform - cubic b-spline transform fitted to
/// displacements defined at scattered points.
///
/// This transform takes a surface or volumetric mesh (any vtkPointSet, such as
/// vtkPolyData or vtkUnstructuredGrid) whose points carry a 3-component
/// displacement vector as point data, and fits a cubic b-spline transform to
/// those displacements.
///
/// The b-spline coefficients are computed with the multilevel b-spline
/// approximation (MBA) algorithm: a coarse control point lattice is fitted
/// first, then the lattice is repeatedly refined (halving the control point
/// spacing) and fitted to the remaining residual, until the residual drops
/// below Tolerance or a refinement limit is reached. The result is smooth
/// everywhere, so it can be evaluated (and inverted) at any point in space,
/// not only at the input mesh points.
///
/// S. Lee, G. Wolberg, and S. Y. Shin, "Scattered data interpolation with
/// multilevel B-splines", IEEE Transactions on Visualization and Computer
/// Graphics, 3(3):228-244, 1997.
///
/// The class derives from vtkOrientedBSplineTransform, therefore the fitted
/// transform can be used anywhere a b-spline transform is expected (for
/// example, it can be stored in a vtkMRMLBSplineTransformNode and written to
/// an ITK transform file). The fitted control point lattice is axis-aligned,
/// so GridDirectionMatrix is left unset (identity).
///
/// Usage:
///   transform->SetInputPointSet(mesh);
///   transform->SetDisplacementArrayName("Displacement");
///   transform->SetTolerance(0.01);
///   transform->Update();
///   // transform->GetResidual() reports the achieved accuracy

#ifndef vtkScatteredBSplineTransform_h
#define vtkScatteredBSplineTransform_h

#include "vtkAddon.h"

#include "vtkOrientedBSplineTransform.h"

class vtkPointSet;

class VTK_ADDON_EXPORT vtkScatteredBSplineTransform : public vtkOrientedBSplineTransform
{
public:
  static vtkScatteredBSplineTransform *New();
  vtkTypeMacro(vtkScatteredBSplineTransform,vtkOrientedBSplineTransform);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /* Explicitly deleted functions belong in the public interface */
  vtkScatteredBSplineTransform(const vtkScatteredBSplineTransform&) = delete;
  void operator=(const vtkScatteredBSplineTransform&) = delete;

  // Description:
  // Make another transform of the same type.
  vtkAbstractTransform *MakeTransform() override;

  // Description:
  // Get the MTime. The modification time of the input point set is taken into
  // account, so that editing the mesh triggers a new fit.
  vtkMTimeType GetMTime() override;

  // Description:
  // Set/Get the mesh (vtkPolyData, vtkUnstructuredGrid, or any other
  // vtkPointSet) whose points define the displacement samples.
  // This transform class will never modify the data.
  // While it is not set, the transform is the identity.
  virtual void SetInputPointSet(vtkPointSet*);
  vtkGetObjectMacro(InputPointSet, vtkPointSet);

  // Description:
  // Set/Get the name of the 3-component point data array that contains the
  // displacement vectors. If not set (default), the active vectors array of
  // the input point set is used.
  vtkSetStringMacro(DisplacementArrayName);
  vtkGetStringMacro(DisplacementArrayName);

  // Description:
  // Set/Get the control point spacing of the coarsest (first) level, in the
  // coordinate system of the input points.
  // If an axis is set to 0 (default) then one tenth of the size of the fitting
  // domain along that axis is used.
  vtkSetVector3Macro(ControlPointSpacing, double);
  vtkGetVector3Macro(ControlPointSpacing, double);

  // Description:
  // Set/Get the smallest control point spacing that refinement may produce.
  // Refinement stops before the spacing would drop below this value along any
  // axis. Default is 0.1.
  vtkSetMacro(MinimumControlPointSpacing, double);
  vtkGetMacro(MinimumControlPointSpacing, double);

  // Description:
  // Set/Get the maximum number of levels. Each level after the first halves
  // the control point spacing, so the memory need grows by about 8x per level.
  // Default is 8.
  vtkSetClampMacro(MaximumNumberOfLevels, int, 1, 16);
  vtkGetMacro(MaximumNumberOfLevels, int);

  // Description:
  // Set/Get the fitting tolerance. Refinement stops when the largest residual
  // displacement error at the input points drops below this value (in the
  // physical units of the input points). Default is 0.01.
  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  // Description:
  // Set/Get if the coarsest level is seeded with a least-squares linear fit of
  // the displacements. This makes the transform extrapolate the overall linear
  // trend outside the sampled region instead of flattening out.
  // Default is off.
  vtkSetMacro(UseLinearApproximation, bool);
  vtkGetMacro(UseLinearApproximation, bool);
  vtkBooleanMacro(UseLinearApproximation, bool);

  // Description:
  // Set/Get the fitting domain (xmin, xmax, ymin, ymax, zmin, zmax).
  // Ignored while AutomaticGridBounds is on.
  vtkSetVector6Macro(GridBounds, double);
  vtkGetVector6Macro(GridBounds, double);

  // Description:
  // Set/Get if the fitting domain is computed from the bounds of the input
  // points (default) or taken from GridBounds.
  vtkSetMacro(AutomaticGridBounds, bool);
  vtkGetMacro(AutomaticGridBounds, bool);
  vtkBooleanMacro(AutomaticGridBounds, bool);

  // Description:
  // Set/Get the distance that the automatic fitting domain extends beyond the
  // bounds of the input points along each axis (in the physical units of the
  // input points). Within this margin the transform follows the displacements
  // smoothly; beyond it, it behaves as set by BorderMode. As the domain follows
  // the points, so does the margin, for example along a sequence of meshes.
  // Ignored while AutomaticGridBounds is off. Default is 0.
  vtkSetClampMacro(GridBoundsMargin, double, 0.0, VTK_DOUBLE_MAX);
  vtkGetMacro(GridBoundsMargin, double);

  // Description:
  // Get the largest displacement error at the input points after the last fit,
  // or -1.0 if no successful fit has been performed.
  // Valid after Update().
  vtkGetMacro(Residual, double);

  // Description:
  // Get the number of refinement levels used by the last fit, or 0 if no
  // successful fit has been performed.
  // Valid after Update().
  vtkGetMacro(NumberOfLevels, int);

protected:
  vtkScatteredBSplineTransform();
  ~vtkScatteredBSplineTransform() override;

  // Description:
  // Fit the b-spline coefficients to the input point set, then let the parent
  // class set up the transformation from the fitted coefficients.
  void InternalUpdate() override;

  // Description:
  // Copy this transform from another of the same type.
  void InternalDeepCopy(vtkAbstractTransform *transform) override;

  // Description:
  // Run the multilevel b-spline approximation and store the resulting control
  // point lattice as the coefficient data of the parent class.
  // Returns true on success. On failure the coefficient data is cleared, which
  // makes the transform behave as identity.
  bool FitCoefficients();

  // Description:
  // Discard the fitted coefficients and reset the grid geometry of the parent
  // class, so that the transform behaves as identity.
  void ClearCoefficients();

  vtkPointSet* InputPointSet;
  char* DisplacementArrayName;

  double ControlPointSpacing[3];
  double MinimumControlPointSpacing;
  int MaximumNumberOfLevels;
  double Tolerance;
  bool UseLinearApproximation;
  double GridBounds[6];
  bool AutomaticGridBounds;
  double GridBoundsMargin;

  double Residual;
  int NumberOfLevels;
};

#endif
