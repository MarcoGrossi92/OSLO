function(gather_oslo_files out_var)

  # FATODE
  file(GLOB FATODE_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/lib/FATODE/*.F90")

  if (NOT LAPACK_FOUND)
    file(GLOB LAPACK_FILES "${CMAKE_CURRENT_SOURCE_DIR}/lib/FATODE/lapack/*.f")
    list(APPEND FATODE_FILES ${LAPACK_FILES})
  endif()

  # HAIRER
  file(GLOB HAIRER_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/lib/hairer/*.f")

  # MAIN FILE
  set(OSLO_MAIN "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/oslo.f90")

  # Combine everything
  set(all_files ${OSLO_MAIN} ${FATODE_FILES} ${HAIRER_FILES})

  # Export to caller
  set(${out_var} ${all_files} PARENT_SCOPE)

endfunction()


function(link_oslo_dependencies target)

  target_compile_definitions(${target} PRIVATE FULL_ALGEBRA)

  # The solvers are called from inside the caller's OpenMP parallel regions.
  # Without the OpenMP flag ifx allocates fixed-size local arrays statically,
  # i.e. shared by all threads, and the COMMON-block THREADPRIVATE directives
  # in the Hairer sources are ignored.
  if (USE_OPENMP)
    target_link_libraries(${target} PUBLIC OpenMP::OpenMP_Fortran)
  endif()

  # ------------------------------
  # LINKING
  # ------------------------------
  if(BLAS_FOUND AND LAPACK_FOUND)
    target_link_libraries(${target} PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  endif()

  if (UNIX AND NOT APPLE)
    file(GLOB INTELODE "${CMAKE_CURRENT_SOURCE_DIR}/lib/Intel-ODE/lib/intel64/*")
    target_link_libraries(${target} PRIVATE ${INTELODE})
  endif()

  if (USE_SUNDIALS)
    target_link_libraries(${target}
      PRIVATE
        SUNDIALS::fcvodes_mod
        SUNDIALS::fidas_mod
        SUNDIALS::fcore_mod
        SUNDIALS::fsunnonlinsolnewton_mod
    )
  endif()

  # ------------------------------
  # MODULE DIR
  # ------------------------------
  target_include_directories(${target} INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

endfunction()