/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

=========================================================================auto=*/

/// \brief Compute orthonormal basis along a curve with minimal torsion.
/// 
/// Parallel transport frame provides a smoothly changing orthonormal basis
/// along a curve. The basic idea is to sweep an initial orthonormal coordinate
/// frame along the line without being affected by torsion: at each point, the unit
/// vectors making up the frame are rotated of the exact amount required by curvature,
///
/// Vectors of the basis: normal (x), binormal (y), tangent (z).
/// Tangent vector direction is always the curve's tangent.
/// Normal vector starts pointing towards a preferred orientation.
/// Binormal is the cross product of tangent and normal.
///
/// Note that the classic Frenet-Serret method can be used to compute othonormal basis along
/// a curve, too. However, in a Frenet-Serret frame the normal is not defined when the
/// curvature vanishes and abruptly changes orientation when the direction of concavity of the curve
/// changes. Therefore, Frenet-Serret frame is not usable for arbitrary curves.
///
/// References:
/// - Parallel transport theory: R. Bishop, "There is more than one way to frame a curve",
///   American Mathematical Monthly, vol. 82, no. 3, pp. 246�251, 1975
/// - Parallel transport implementation: Piccinelli M, Veneziani A, Steinman DA, Remuzzi A, Antiga L.
///   "A framework for geometric analysis of vascular structures: application to cerebral aneurysms.",
///   IEEE Trans Med Imaging. 2009 Aug;28(8):1141-55. doi: 10.1109/TMI.2009.2021652.
/// 
/// The initial implementation was based on VMTK (vtkvmtkCenterlineAttributesFilter) which was optimized
/// and enhanced with more predictable initial normal vector direction. In the future, support for closed
/// curves may be added (where constant torsion need to be added or subtracted to make the normal directions
/// match between the start and end point of the curve).

#ifndef vtkParallelTransportFrame_h
#define vtkParallelTransportFrame_h

#include "vtkAddon.h"  // For export macro
#include "vtkPolyDataAlgorithm.h"

class VTK_ADDON_EXPORT vtkParallelTransportFrame : public vtkPolyDataAlgorithm
{
public:
  void PrintSelf(ostream& os, vtkIndent indent) override;
  vtkTypeMacro(vtkParallelTransportFrame, vtkPolyDataAlgorithm);
  static vtkParallelTransportFrame* New();

  ///@{
  /// Get/set the point array name that contains the tangent (z) axis.
  /// Default value is "Tangents"
  vtkSetStringMacro(TangentsArrayName);
  vtkGetStringMacro(TangentsArrayName);
  ///@} 

  ///@{
  /// Get/set the point array name that contains the normals (x) axis.
  /// Default value is "Normals"
  vtkSetStringMacro(NormalsArrayName);
  vtkGetStringMacro(NormalsArrayName);
  ///@} 

  ///@{
  /// Get/set the point array name that contains the binormal (y) axis.
  /// Default value is "Binormals"
  vtkSetStringMacro(BinormalsArrayName);
  vtkGetStringMacro(BinormalsArrayName);
  ///@} 

  ///@{
  /// Use rotation minimizing frames, otherwise use Bishop frames
  /// By default is False
  vtkSetMacro(RotationMinimizingFrames, vtkTypeBool);
  vtkGetMacro(RotationMinimizingFrames, vtkTypeBool);
  ///@} 

  /// Define the preferred direction of the normal vector at the first point of the curve.
  /// It is just "preferred" because the direction has to be orhogonal to the tangent,
  /// so in general the normal vector cannot point into exactly to a required direction.
  /// By default it is (1, 0, 0).
  /// \sa PreferredInitialBinormalVector
  vtkSetVectorMacro(PreferredInitialNormalVector, double, 3);
  vtkGetVectorMacro(PreferredInitialNormalVector, double, 3);
  ///@}

  /// Define the preferred direction of the binormal vector at the first point.
  /// It is used only if the curve's tangent at the first point is parallel to the
  /// preferred initial normal vector.
  /// By default it is (0, 1, 0).
  /// \sa PreferredInitialNormalVector
  vtkSetVectorMacro(PreferredInitialBinormalVector, double, 3);
  vtkGetVectorMacro(PreferredInitialBinormalVector, double, 3);
  ///@}

protected:
  vtkParallelTransportFrame();
  ~vtkParallelTransportFrame() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

  void ComputeAxisDirections(vtkPolyData* input, vtkIdType cellIndex, vtkDoubleArray* tangentsArray, vtkDoubleArray* normalsArray, vtkDoubleArray* binormalsArray);
  void ComputeAxisDirections2(vtkPolyData* input, vtkIdType cellIndex, vtkDoubleArray* tangentsArray, vtkDoubleArray* normalsArray, vtkDoubleArray* binormalsArray);

  /// Rotate a vector around an axis
  static void RotateVector(double* inVector, double* outVector, const double* axis, double angle);

private:
  vtkParallelTransportFrame(const vtkParallelTransportFrame&) = delete;
  void operator=(const vtkParallelTransportFrame&) = delete;

  char* TangentsArrayName = nullptr;
  char* NormalsArrayName = nullptr;
  char* BinormalsArrayName = nullptr;

  vtkTypeBool RotationMinimizingFrames = false;

  /// Tolerance value used for checking that a value is non-zero.
  double Tolerance = 1e-6;
  /// Minimum distance for comuting initial tangent direction
  double MinimumDistance = 1e-3;
  
  double PreferredInitialNormalVector[3] = { 1.0, 0.0, 0.0 };
  double PreferredInitialBinormalVector[3] = { 0.0, 0.0, 0.0 };
};

#endif // vtkParallelTransportFrame_h
