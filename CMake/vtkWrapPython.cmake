#
# CMake implementation of the Wrap Python command.
#
macro(VTK_WRAP_PYTHON2 TARGET SOURCE_LIST_NAME)
  # convert to the WRAP3 signature
  vtk_wrap_python3(${TARGET} ${SOURCE_LIST_NAME} "${ARGN}")
endmacro()

macro(VTK_WRAP_PYTHON3 TARGET SRC_LIST_NAME SOURCES)
  if(NOT VTK_WRAP_PYTHON_INIT_EXE)
    if(TARGET vtkWrapPythonInit)
      set(VTK_WRAP_PYTHON_INIT_EXE vtkWrapPythonInit)
    elseif(TARGET VTK::WrapPythonInit)
      set(VTK_WRAP_PYTHON_INIT_EXE VTK::WrapPythonInit)
    else()
      message(SEND_ERROR
        "VTK_WRAP_PYTHON_INIT_EXE not specified when calling vtk_wrap_python")
    endif()
  endif()
  if(NOT VTK_WRAP_PYTHON_EXE)
    if(TARGET vtkWrapPython)
      set(VTK_WRAP_PYTHON_EXE vtkWrapPython)
    elseif(TARGET VTK::WrapPython)
      set(VTK_WRAP_PYTHON_EXE VTK::WrapPython)
    else()
      message(SEND_ERROR
        "VTK_WRAP_PYTHON_EXE not specified when calling vtk_wrap_python")
    endif()
  endif()

  # start writing the input file for the init file
  set(VTK_WRAPPER_INIT_DATA "${TARGET}")

  # collect the common wrapper-tool arguments
  set(_common_args)
  foreach(file IN LISTS VTK_WRAP_HINTS)
    set(_common_args "${_common_args}--hints \"${file}\"\n")
  endforeach()

  # write wrapper-tool arguments to a file
  set(_args_file ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}.$<CONFIGURATION>.args)
  file(GENERATE OUTPUT ${_args_file} CONTENT "${_common_args}
$<$<BOOL:${KIT_HIERARCHY_FILE}>:
--types \"$<JOIN:${KIT_HIERARCHY_FILE},\"
--types \">\">
$<$<BOOL:${OTHER_HIERARCHY_FILES}>:
--types \"$<JOIN:${OTHER_HIERARCHY_FILES},\"
--types \">\">
$<$<BOOL:$<TARGET_PROPERTY:${TARGET},COMPILE_DEFINITIONS>>:
-D\"$<JOIN:$<TARGET_PROPERTY:${TARGET},COMPILE_DEFINITIONS>,\"
-D\">\">
$<$<BOOL:$<TARGET_PROPERTY:${TARGET},INCLUDE_DIRECTORIES>>:
-I\"$<JOIN:$<TARGET_PROPERTY:${TARGET},INCLUDE_DIRECTORIES>,\"
-I\">\">
")

  if (CMAKE_GENERATOR MATCHES "Ninja")
    set(hierarchy_depend ${KIT_HIERARCHY_FILE})
  else ()
    string(LENGTH "${TARGET}" target_length)
    math(EXPR target_length "${target_length} - 6")
    string(SUBSTRING "${TARGET}" 0 "${target_length}" target_basename)
    get_property(is_kit GLOBAL
      PROPERTY "_vtk_${target_basename}_is_kit")
    if (is_kit)
      set(hierarchy_depend)
      get_property(kit_modules GLOBAL
        PROPERTY "_vtk_${target_basename}_kit_modules")
      message("depends for ${TARGET}: ${kit_modules}")
      foreach (depend IN LISTS kit_modules)
        if(TARGET "${depend}Hierarchy")
          list(APPEND hierarchy_depend
            "${depend}Hierarchy")
        endif()
      endforeach ()
    else ()
      set(hierarchy_depend "${target_basename}Hierarchy")
    endif ()
  endif ()

  # for each class
  foreach(FILE ${SOURCES})
    # what is the filename without the extension
    get_filename_component(TMP_FILENAME ${FILE} NAME_WE)

    # the input file might be full path so handle that
    get_filename_component(TMP_FILEPATH ${FILE} PATH)

    # compute the input filename
    if(TMP_FILEPATH)
      set(TMP_INPUT ${TMP_FILEPATH}/${TMP_FILENAME}.h)
    else()
      set(TMP_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${TMP_FILENAME}.h)
    endif()

    # add the info to the init file
    set(VTK_WRAPPER_INIT_DATA
      "${VTK_WRAPPER_INIT_DATA}\n${TMP_FILENAME}")

    # new source file is namePython.cxx, add to resulting list
    list(APPEND ${SRC_LIST_NAME} ${TMP_FILENAME}Python.cxx)

    # add custom command to output
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}Python.cxx
      DEPENDS ${VTK_WRAP_PYTHON_EXE}
              ${VTK_WRAP_HINTS}
              ${TMP_INPUT}
              ${_args_file}
              ${hierarchy_depend}
      IMPLICIT_DEPENDS CXX ${TMP_INPUT}
      COMMAND ${VTK_WRAP_PYTHON_EXE}
              @${_args_file}
              -o ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}Python.cxx
              ${TMP_INPUT}
      COMMENT "Generating Python wrapper sources for ${TMP_FILENAME}"
      VERBATIM
      )
  endforeach()

  if(${VTK_VERSION} VERSION_GREATER_EQUAL "8.90")
    # find the python modules needed by this python module
    string(REGEX REPLACE "Python$" "" _basename "${TARGET}")
    string(REGEX REPLACE "Kit$" "" kit_basename "${_basename}")
    if(_${kit_basename}_is_kit)
      set(${_basename}_WRAP_DEPENDS)
      foreach(kit_module IN LISTS _${kit_basename}_modules)
        list(APPEND ${_basename}_WRAP_DEPENDS
          ${${kit_module}_WRAP_DEPENDS})
      endforeach()
    endif()
    set(_python_module_depends)
    foreach(dep IN LISTS ${_basename}_WRAP_DEPENDS)
      if(dep STREQUAL "${_basename}" OR dep STREQUAL "${kit_basename}")
        continue()
      endif()
      if(VTK_ENABLE_KITS AND ${dep}_KIT)
        if(NOT ${dep}_KIT STREQUAL kit_basename)
          list(APPEND _python_module_depends ${${dep}_KIT}KitPython)
        endif()
      elseif(dep MATCHES "^VTK::")
        # Import the Python module of each VTK dependency so that base classes
        # defined in VTK (for example vtkGridTransform in vtkFiltersHybrid) are
        # registered before this library's wrapped subclasses are built. The
        # eager "import vtk" used to load every VTK module up front; with lazy
        # VTK loading it does not, so a wrapped subclass such as
        # vtkOrientedGridTransform would otherwise be built without its VTK base
        # class and silently lose the inherited methods. Follow VTK's own module
        # wrapper (vtkModuleWrapPython): skip modules excluded from wrapping and
        # name the rest "<python_package>.<library_name>".
        _vtk_module_get_module_property("${dep}"
          PROPERTY  "exclude_wrap"
          VARIABLE  _dep_exclude_wrap)
        if(NOT _dep_exclude_wrap)
          _vtk_module_get_module_property("${dep}"
            PROPERTY  "python_package"
            VARIABLE  _dep_python_package)
          _vtk_module_get_module_property("${dep}"
            PROPERTY  "library_name"
            VARIABLE  _dep_library_name)
          if(_dep_python_package AND _dep_library_name)
            list(APPEND _python_module_depends "${_dep_python_package}.${_dep_library_name}")
          endif()
        endif()
      elseif(TARGET ${dep}Python)
        list(APPEND _python_module_depends ${dep}Python)
      endif()
    endforeach()
    if(_python_module_depends)
      # add a DEPENDS section to the Init.data file for this module
      set(VTK_WRAPPER_INIT_DATA "${VTK_WRAPPER_INIT_DATA}\nDEPENDS")
      list(REMOVE_DUPLICATES _python_module_depends)
      foreach(dep IN LISTS _python_module_depends)
        set(VTK_WRAPPER_INIT_DATA "${VTK_WRAPPER_INIT_DATA}\n${dep}")
      endforeach()
    endif()
  endif()

  # finish the data file for the init file
  configure_file(
    ${vtkAddon_CMAKE_DIR}/vtkWrapperInit.data.in
    ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.data
    @ONLY
    )

  # The DEPENDS just written lists every VTK python module on the link line,
  # because the wrapping dependencies are seeded from ${VTK_LIBRARIES}. Only the
  # modules that own a direct base class of a wrapped class actually need to be
  # imported first (each VTK module imports its own bases transitively). Prune
  # the DEPENDS to that minimal set so lazy VTK loading is not defeated by
  # importing all of VTK when this module is imported. Cross-module Slicer
  # dependencies are preserved unchanged. Requires the merged hierarchy file and
  # a Python interpreter; fall back to the unpruned data file when unavailable.
  set(_init_data_file "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.data")
  if(KIT_HIERARCHY_FILE AND PYTHON_EXECUTABLE AND ${VTK_VERSION} VERSION_GREATER_EQUAL "8.90")
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.pruned.data
      DEPENDS ${vtkAddon_CMAKE_DIR}/SlicerPrunePythonModuleDepends.py
              ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.data
              ${KIT_HIERARCHY_FILE}
      COMMAND ${PYTHON_EXECUTABLE}
              ${vtkAddon_CMAKE_DIR}/SlicerPrunePythonModuleDepends.py
              ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.data
              ${KIT_HIERARCHY_FILE}
              ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.pruned.data
      COMMENT "Pruning Python module dependencies for ${TARGET}"
      VERBATIM
      )
    set(_init_data_file "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.pruned.data")
  endif()

  set(_init_impl_src "")
  if(${VTK_VERSION} VERSION_LESS "9.5.0")
    set(_init_impl_src "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}InitImpl.cxx")
  endif()
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.cxx
           ${_init_impl_src}
    DEPENDS ${VTK_WRAP_PYTHON_INIT_EXE}
            ${_init_data_file}
    COMMAND ${VTK_WRAP_PYTHON_INIT_EXE}
            ${_init_data_file}
            ${CMAKE_CURRENT_BINARY_DIR}/${TARGET}Init.cxx
            ${_init_impl_src}
    COMMENT "Generating the Python module initialization sources for ${TARGET}"
    VERBATIM
    )
  if(${VTK_VERSION} VERSION_LESS "9.5.0")
    # Create the Init File
    set(${SRC_LIST_NAME} ${${SRC_LIST_NAME}} ${TARGET}InitImpl.cxx)
  endif()

endmacro()

if(VTK_WRAP_PYTHON_FIND_LIBS)
  find_package(Python3 COMPONENTS Development.Module REQUIRED)
endif()

# Determine the location of the supplied header in the include_dirs supplied.
macro(vtk_find_header header include_dirs full_path)
  unset(${full_path})
  foreach(_dir ${include_dirs})
    if(EXISTS "${_dir}/${header}")
      set(${full_path} "${_dir}/${header}")
      break()
    endif()
  endforeach()
endmacro()

# Macro that just takes the a list of module names, figure the rest out from there.
macro(vtk_wrap_python TARGET SRC_LIST_NAME)
  # List of all headers to wrap.
  set(headers_to_wrap)

  foreach(module ${ARGN})

  # Decide what to do for each header.
  foreach(header ${${module}_HEADERS})
    # Find the full path to the header file to be wrapped.
    vtk_find_header(${header}.h "${${module}_INCLUDE_DIRS}" class_header_path)
    if(NOT class_header_path)
      message(FATAL_ERROR "Could not find ${header}.h for Python wrapping!")
    endif()

    # The new list of headers has the full path to each file.
    list(APPEND headers_to_wrap ${class_header_path})
  endforeach()

  endforeach() # end ARGN loop

  # Delegate to vtk_wrap_python3
  vtk_wrap_python3(${TARGET} ${SRC_LIST_NAME} "${headers_to_wrap}" ${ARGN})
endmacro()
