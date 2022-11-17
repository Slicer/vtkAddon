
#
# TODO
# * Support vtkAddon_USE_UTF8. See https://gitlab.kitware.com/vtk/vtk/-/issues/18542
# * Support testing and use of vtkAddon_LAUNCH_COMMAND
#

set(classes
  vtkAddonMathUtilities
  vtkAddonTestingUtilities
  vtkCurveGenerator
  vtkErrorSink
  vtkLinearSpline
  vtkOrientedBSplineTransform
  vtkOrientedGridTransform
  vtkParallelTransportFrame
  vtkParametricPolynomialApproximation
  vtkPersonInformation
  vtkRawRGBVolumeCodec
  vtkSlicerDijkstraGraphGeodesicPath
  vtkStreamingVolumeCodec
  vtkStreamingVolumeCodecFactory
  vtkStreamingVolumeFrame
  vtkTestingOutputWindow
)

if(TARGET VTK::RenderingOpenGL2)
  list(APPEND classes
    vtkOpenGLShaderComputation
    vtkOpenGLTextureImage
  )
endif()

set(sources
  vtkAddonTestingUtilities.txx
)

set(headers
  vtkAddonSetGet.h
  vtkLoggingMacros.h
)

vtk_module_add_module(VTK::Addon
  EXPORT_MACRO_PREFIX VTK_ADDON
  CLASSES ${classes}
  SOURCES ${sources}
  HEADERS ${headers}
)

# Install vtkAddon CMake files
if(DEFINED vtk_cmake_destination)
  install(
    FILES
      ${CMAKE_SOURCE_DIR}/CMake/vtkAddonFunctionAddExecutable.cmake
      ${CMAKE_SOURCE_DIR}/CMake/vtkAddonFunctionAddExecutable.cmake
      ${CMAKE_SOURCE_DIR}/CMake/vtkMacroKitPythonWrap.cmake
      ${CMAKE_SOURCE_DIR}/CMake/vtkWrapHierarchy.cmake
      ${CMAKE_SOURCE_DIR}/CMake/vtkWrapperInit.data.in
      ${CMAKE_SOURCE_DIR}/CMake/vtkWrapPython.cmake
      ${CMAKE_SOURCE_DIR}/CMake/WindowsApplicationUseUtf8.manifest
    DESTINATION ${vtk_cmake_destination}
    COMPONENT "development")
endif()
