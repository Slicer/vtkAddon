#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.13.4)
#-----------------------------------------------------------------------------

project(vtkAddon
  LANGUAGES CXX C
  DESCRIPTION "General-purpose features that may be integrated into VTK library in the future."
)


set(vtkAddon_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
set(CMAKE_MODULE_PATH ${vtkAddon_CMAKE_DIR} ${CMAKE_MODULE_PATH})

include(CMakeDependentOption)
include(vtkAddonFunctionAddExecutable)

# --------------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------------

#
# VTK
#
find_package(VTK REQUIRED)
if(${VTK_VERSION} VERSION_LESS "8.90")
  include(${VTK_USE_FILE})
endif()
set(vtkAddon_LIBS ${VTK_LIBRARIES})

# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------
if(NOT DEFINED BUILD_SHARED_LIBS)
  option(BUILD_SHARED_LIBS "Build with shared libraries." ON)
endif()

if(NOT DEFINED BUILD_TESTING)
  option(BUILD_TESTING "Test the project" OFF)
endif()

if(NOT DEFINED vtkAddon_USE_UTF8)
  option(vtkAddon_USE_UTF8 "Make applications use UTF-8 as code page." OFF)
endif()

if(NOT DEFINED vtkAddon_WRAP_PYTHON)
  cmake_dependent_option(vtkAddon_WRAP_PYTHON "Build python wrapping" OFF "VTK_WRAP_PYTHON" OFF)
endif()

# --------------------------------------------------------------------------
# Configure headers
# --------------------------------------------------------------------------
set(configure_header_file vtkAddonConfigure.h)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${configure_header_file}.in
  ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
  )

# --------------------------------------------------------------------------
# Install headers
# --------------------------------------------------------------------------

if(NOT DEFINED ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  set(${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT ON)
endif()

if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)

  if(NOT DEFINED ${PROJECT_NAME}_INSTALL_INCLUDE_DIR)
    set(${PROJECT_NAME}_INSTALL_INCLUDE_DIR include/${PROJECT_NAME})
  endif()

  file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  install(
    FILES ${headers} ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
    DESTINATION ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR} COMPONENT Development)

  file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
  install(
    FILES ${headers} ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
    DESTINATION ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR} COMPONENT Development)
endif()

# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------
set(vtkAddon_SRCS
  vtkAddonTestingUtilities.cxx
  vtkAddonTestingUtilities.h
  vtkAddonTestingUtilities.txx
  vtkErrorSink.cxx
  vtkErrorSink.h
  vtkLoggingMacros.h
  vtkTestingOutputWindow.cxx
  vtkTestingOutputWindow.h
  vtkOrientedBSplineTransform.cxx
  vtkOrientedBSplineTransform.h
  vtkOrientedGridTransform.cxx
  vtkOrientedGridTransform.h
  vtkPersonInformation.cxx
  vtkPersonInformation.h
  vtkAddonMathUtilities.h
  vtkAddonMathUtilities.cxx
  vtkAddonSetGet.h
  vtkStreamingVolumeCodec.cxx
  vtkStreamingVolumeCodec.h
  vtkStreamingVolumeFrame.cxx
  vtkStreamingVolumeFrame.h
  vtkStreamingVolumeCodecFactory.cxx
  vtkStreamingVolumeCodecFactory.h
  vtkRawRGBVolumeCodec.cxx
  vtkRawRGBVolumeCodec.h
)

if(VTK_RENDERING_BACKEND STREQUAL "OpenGL2")
  set(vtkAddon_SRCS
    ${vtkAddon_SRCS}
    vtkOpenGLTextureImage.cxx
    vtkOpenGLTextureImage.h
    vtkOpenGLShaderComputation.cxx
    vtkOpenGLShaderComputation.h
  )
  find_package(OpenGL)
  set(vtkAddon_LIBS
    ${vtkAddon_LIBS}
    ${OPENGL_LIBRARIES}
  )
endif()

# Abstract/pure virtual classes

#set_source_files_properties(
#  ABSTRACT
#  )

# Helper classes

set_source_files_properties(
  vtkAddonTestingUtilities.h
  vtkLoggingMacros.h 
  vtkAddonSetGet.h
  WRAP_EXCLUDE
  )
# --------------------------------------------------------------------------
# Include dirs
# --------------------------------------------------------------------------
set(include_dirs
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )
include_directories(BEFORE ${include_dirs})

# --------------------------------------------------------------------------
# Build library
# --------------------------------------------------------------------------
set(lib_name ${PROJECT_NAME})

set(srcs ${vtkAddon_SRCS})
add_library(${lib_name} ${srcs})

target_link_libraries(${lib_name} ${vtkAddon_LIBS})

# --------------------------------------------------------------------------
# Folder
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_FOLDER)
  set(${PROJECT_NAME}_FOLDER ${PROJECT_NAME})
endif()
if(NOT "${${PROJECT_NAME}_FOLDER}" STREQUAL "")
  set_target_properties(${lib_name} PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
endif()

# --------------------------------------------------------------------------
# Export target
# --------------------------------------------------------------------------

# Export targets for a build tree
if(NOT DEFINED ${PROJECT_NAME}_EXPORT_FILE)
  set(${PROJECT_NAME}_EXPORT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)
endif()
export(TARGETS ${lib_name} APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})

# --------------------------------------------------------------------------
# Install library
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_BIN_DIR)
  set(${PROJECT_NAME}_INSTALL_BIN_DIR bin)
endif()
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_LIB_DIR)
  set(${PROJECT_NAME}_INSTALL_LIB_DIR lib/${PROJECT_NAME})
endif()

install(TARGETS ${lib_name}
  RUNTIME DESTINATION ${${PROJECT_NAME}_INSTALL_BIN_DIR} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${${PROJECT_NAME}_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${${PROJECT_NAME}_INSTALL_LIB_DIR} COMPONENT Development
  )

# --------------------------------------------------------------------------
# Python Wrapping
# --------------------------------------------------------------------------
if(vtkAddon_WRAP_PYTHON)

  if(NOT DEFINED ${PROJECT_NAME}_INSTALL_PYTHON_MODULE_LIB_DIR)
    set(${PROJECT_NAME}_INSTALL_PYTHON_MODULE_LIB_DIR ${${PROJECT_NAME}_INSTALL_LIB_DIR})
  endif()
  if(NOT DEFINED ${PROJECT_NAME}_INSTALL_PYTHON_LIB_DIR)
    set(${PROJECT_NAME}_INSTALL_PYTHON_LIB_DIR ${${PROJECT_NAME}_INSTALL_LIB_DIR})
  endif()

  include(vtkMacroKitPythonWrap)

  vtkMacroKitPythonWrap(
    KIT_NAME ${lib_name}
    KIT_SRCS ${vtkAddon_SRCS}
    KIT_INSTALL_BIN_DIR ${${PROJECT_NAME}_INSTALL_BIN_DIR}
    KIT_INSTALL_LIB_DIR ${${PROJECT_NAME}_INSTALL_PYTHON_LIB_DIR}
    KIT_MODULE_INSTALL_LIB_DIR ${${PROJECT_NAME}_INSTALL_PYTHON_MODULE_LIB_DIR}
    )
  # Export target
  export(TARGETS ${lib_name}Python APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})
  if(${VTK_VERSION} VERSION_LESS "8.90")
    export(TARGETS ${lib_name}PythonD APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})
  endif()
  # Folder
  if(NOT "${${PROJECT_NAME}_FOLDER}" STREQUAL "")
    set_target_properties(${lib_name}Python PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
    if(${VTK_VERSION} VERSION_LESS "8.90")
      set_target_properties(${lib_name}PythonD PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
    endif()
    if(TARGET ${lib_name}Hierarchy)
      set_target_properties(${lib_name}Hierarchy PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
    endif()
  endif()
endif()

# --------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(Testing)
endif()

# --------------------------------------------------------------------------
# Set INCLUDE_DIRS variable
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include dirs" FORCE)

# --------------------------------------------------------------------------
# Configure support file
# --------------------------------------------------------------------------

# Configure 'vtkAddonConfig.cmake' for a build tree
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkAddonConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/vtkAddonConfig.cmake
  )

# Configuret vtkAddonConfig.cmake for an install tree
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkAddonInstallConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/install/vtkAddonConfig.cmake
  )


# --------------------------------------------------------------------------
# Install configuration files
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_CMAKE_DIR)
  set(${PROJECT_NAME}_INSTALL_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMake)
endif()

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/install/vtkAddonConfig.cmake
  DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKE_DIR}
  )

