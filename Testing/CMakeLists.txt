set(KIT vtkAddon)

create_test_sourcelist(Tests ${KIT}CxxTests.cxx
  vtkAddonMathUtilitiesTest1.cxx
  vtkAddonTestingUtilitiesTest1.cxx
  vtkLoggingMacrosTest1.cxx
  vtkParallelTransportTest1.cxx
  vtkPersonInformationTest1.cxx
  )

set(LIBRARY_NAME ${PROJECT_NAME})

vtkaddon_add_executable(${KIT}CxxTests ${Tests})
target_link_libraries(${KIT}CxxTests ${lib_name})

set_target_properties(${KIT}CxxTests PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})

function(vtkaddon_add_test testname)
  add_test(NAME ${testname} COMMAND ${vtkAddon_LAUNCH_COMMAND} $<TARGET_FILE:${KIT}CxxTests> ${testname} ${ARGN})
  set_property(TEST ${testname} PROPERTY LABELS ${KIT})
endfunction()

vtkaddon_add_test( vtkAddonMathUtilitiesTest1 )
vtkaddon_add_test( vtkAddonTestingUtilitiesTest1 )
vtkaddon_add_test( vtkLoggingMacrosTest1 )
vtkaddon_add_test( vtkPersonInformationTest1 )
