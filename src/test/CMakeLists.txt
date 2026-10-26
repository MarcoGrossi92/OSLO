file(GLOB TESTS CONFIGURE_DEPENDS *.f90 *.F90)

# robertsonOMP runs the solvers from an OpenMP parallel region: without OpenMP it has nothing to test
if(NOT USE_OPENMP)
  list(FILTER TESTS EXCLUDE REGEX "/robertsonOMP[.]f90$")
endif()

set(TESTS ${TESTS} PARENT_SCOPE)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")

# The problems, reference solutions and checks shared by the tests
add_library(oslo_test_support STATIC support/oslo_test_support.f90)
set_target_properties(oslo_test_support PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/support")
target_include_directories(oslo_test_support PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/modules/support")
# OSLO's OpenMP flag: the right-hand sides and Jacobians here are called by concurrent solvers, and without the flag
# ifx allocates their temporaries (the array constructors of robertson_jac) statically, shared by all threads
target_link_libraries(oslo_test_support PRIVATE OSLO)

foreach(TEST ${TESTS})
  get_filename_component(EXE_NAME ${TEST} NAME_WE)
  add_executable(${EXE_NAME} ${TEST})
  target_link_libraries(${EXE_NAME} PRIVATE OSLO oslo_test_support)
  # Each program has its own module folder, so that the modules of one never clash with those of another in a
  # parallel build
  set_target_properties(${EXE_NAME} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/${EXE_NAME}")
  target_include_directories(${EXE_NAME}
      PRIVATE
          $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  )
  # Each program checks its results and stops with a nonzero exit status when one is wrong
  add_test(NAME OSLO.${EXE_NAME} COMMAND ${EXE_NAME})
endforeach()

# Four threads whatever the machine, so that the solvers do run concurrently
if(USE_OPENMP)
  set_tests_properties(OSLO.robertsonOMP PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
