file(GLOB TESTS CONFIGURE_DEPENDS *.f90 *.F90)

# robertsonOMP runs the solvers from an OpenMP parallel region: without OpenMP it has nothing to test
if(NOT USE_OPENMP)
  list(FILTER TESTS EXCLUDE REGEX "/robertsonOMP[.]f90$")
endif()
# fatode_ros_disabled checks that setup_odesolver refuses FATODE's Rosenbrock solvers, in the builds that disable them
if(NOT FATODE_ROS_DISABLED)
  list(FILTER TESTS EXCLUDE REGEX "/fatode_ros_disabled[.]f90$")
endif()

set(TESTS ${TESTS} PARENT_SCOPE)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")

# The problems, reference solutions and checks shared by the tests
add_library(oslo_test_support STATIC support/oslo_test_support.f90)
set_target_properties(oslo_test_support PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/support")
target_include_directories(oslo_test_support PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/modules/support")
# The OpenMP flag: the right-hand sides and Jacobians here are called by concurrent solvers, and without the flag ifx
# allocates their temporaries (the array constructors of robertson_jac) statically, shared by all threads
if(USE_OPENMP)
  target_link_libraries(oslo_test_support PRIVATE OpenMP::OpenMP_Fortran)
endif()

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
  if(NOT EXE_NAME STREQUAL "fatode_ros_disabled")
    add_test(NAME OSLO.${EXE_NAME} COMMAND ${EXE_NAME})
  endif()
endforeach()

# Once for each solver disabled: the test passes when setup_odesolver stops with its message
if(FATODE_ROS_DISABLED)
  foreach(SOLVER ros2 ros3 ros4 rodas3 rodas4)
    add_test(NAME OSLO.fatode_ros_disabled.${SOLVER} COMMAND fatode_ros_disabled ${SOLVER})
    set_tests_properties(OSLO.fatode_ros_disabled.${SOLVER} PROPERTIES PASS_REGULAR_EXPRESSION "is disabled in this build")
  endforeach()
endif()

# Four threads whatever the machine, so that the solvers do run concurrently
if(USE_OPENMP)
  set_tests_properties(OSLO.robertsonOMP PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
