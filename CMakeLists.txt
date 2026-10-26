cmake_minimum_required(VERSION 3.23)

project (OSLO LANGUAGES Fortran C)

set(CMAKE_VERBOSE_MAKEFILE FALSE)
# seach path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING. You should review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries. There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

# Include the functions that gather source files and link dependencies for OSLO.
INCLUDE(${CMAKE_MODULE_PATH}/SetOSLO.cmake)

# Set directories
include(GNUInstallDirs)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# ------------------------------
# Define third-party libs
# ------------------------------

# Intel-ODE
if(UNIX AND NOT APPLE)
  add_compile_definitions(INTEL)
  file(GLOB INTELODE "${CMAKE_CURRENT_SOURCE_DIR}/lib/Intel-ODE/lib/intel64/*")
endif()

# SUNDIALS
OPTION(USE_SUNDIALS "Use Sundials suite of solvers" ON)
# The SUNDIALS installation, built from lib/sundials when it is missing or stale: by default sundials-install in these
# sources. A folder outside them keeps SUNDIALS for the next copy of OSLO (a project that fetches OSLO clones it in
# each of its build trees)
set(OSLO_SUNDIALS_DIR "" CACHE PATH "SUNDIALS installation (empty: sundials-install in OSLO's sources)")
if (USE_SUNDIALS)
  if (OSLO_SUNDIALS_DIR)
    set(SUNDIALS_INSTALL_DIR "${OSLO_SUNDIALS_DIR}")
  else()
    set(SUNDIALS_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sundials-install")
  endif()
  # Check if already installed
  set(SUNDIALS_CONFIG_FILE "${SUNDIALS_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/cmake/sundials/SUNDIALSConfig.cmake")
  # The SUNDIALS Fortran interface modules are called from inside the caller's
  # OpenMP parallel regions: build them with the OpenMP flag, otherwise ifx makes
  # their derived-type locals (the c_ptr argument wrappers) static and shared by
  # all threads.
  set(SUNDIALS_Fortran_FLAGS "")
  if (USE_OPENMP)
    set(SUNDIALS_Fortran_FLAGS "${OpenMP_Fortran_FLAGS}")
  endif()
  # The stamp of the install: the SUNDIALS commit OSLO records (readable without the submodule), the compilers and
  # the flags. An install with another stamp (another SUNDIALS, compiler or flags) is rebuilt
  execute_process(COMMAND git rev-parse HEAD:lib/sundials
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE SUNDIALS_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  set(SUNDIALS_STAMP "sundials ${SUNDIALS_COMMIT}\n"
                     "C ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}\n"
                     "Fortran ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} ${SUNDIALS_Fortran_FLAGS}\n")
  string(CONCAT SUNDIALS_STAMP ${SUNDIALS_STAMP})
  set(SUNDIALS_STAMP_FILE "${SUNDIALS_INSTALL_DIR}/oslo-stamp.txt")
  set(SUNDIALS_INSTALLED_STAMP "<none>")
  if(EXISTS "${SUNDIALS_STAMP_FILE}")
    file(READ "${SUNDIALS_STAMP_FILE}" SUNDIALS_INSTALLED_STAMP)
  endif()
  if(NOT EXISTS "${SUNDIALS_CONFIG_FILE}" OR NOT SUNDIALS_INSTALLED_STAMP STREQUAL SUNDIALS_STAMP)
    message(STATUS "SUNDIALS not installed (or built from other sources, compilers or flags), building and installing it...")

    # The sources, only to build them: shallow, without the submodules of SUNDIALS (CI configurations, test answers)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials/CMakeLists.txt")
      message(STATUS "Sundials not found. Attempting to initialize...")
      execute_process(COMMAND git submodule update --init --depth 1 lib/sundials
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      RESULT_VARIABLE GIT_SUBMOD_RESULT)
      if(NOT GIT_SUBMOD_RESULT EQUAL 0)
        execute_process(COMMAND git submodule update --init lib/sundials
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
      endif()
      if(NOT GIT_SUBMOD_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to initialize git submodules.")
      endif()
    endif()

    # Build SUNDIALS in a temporary build directory
    execute_process(COMMAND ${CMAKE_COMMAND} -S ${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials
                                          -B ${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials/build
                                          -DCMAKE_INSTALL_PREFIX=${SUNDIALS_INSTALL_DIR}
                                          -DCMAKE_BUILD_TYPE=Release
                                          -DBUILD_SHARED_LIBS=OFF
                                          -DSUNDIALS_BUILD_WITH_MONITORING=OFF
                                          -DSUNDIALS_BUILD_WITH_PROFILING=OFF
                                          -DBUILD_FORTRAN_MODULE_INTERFACE=ON 
                                          -DEXAMPLES_ENABLE_C=OFF
                                          -DEXAMPLES_ENABLE_CXX=OFF
                                          -DEXAMPLES_ENABLE_F2003=OFF
                                          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                          -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                                          -DCMAKE_Fortran_FLAGS=${SUNDIALS_Fortran_FLAGS}
                                          RESULT_VARIABLE CONFIGURE_RESULT)
    if(NOT CONFIGURE_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to configure SUNDIALS")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_SOURCE_DIR}/lib/sundials/build --target install
                    RESULT_VARIABLE BUILD_RESULT)
    if(NOT BUILD_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to build and install SUNDIALS")
    endif()
    file(WRITE "${SUNDIALS_STAMP_FILE}" "${SUNDIALS_STAMP}")
  endif()
  find_package(SUNDIALS REQUIRED PATHS ${SUNDIALS_INSTALL_DIR} NO_DEFAULT_PATH)
  add_compile_definitions(SUNDIALS)
endif(USE_SUNDIALS)

# -----------------------------
# Build library
# ------------------------------

# Gather source files and link dependencies. The functions that do this are defined in cmake/AddOSLO.cmake
find_package(BLAS)
find_package(LAPACK)
gather_oslo_files(OSLO_SRC)
add_library(OSLO ${OSLO_SRC})
link_oslo_dependencies(OSLO)
if (USE_SUNDIALS)
  set_target_properties(OSLO PROPERTIES OSLO_SUNDIALS_DIR "${SUNDIALS_INSTALL_DIR}")
endif()

# ------------------------------
# Build tests and applications
# ------------------------------

# Build tests and applications only if OSLO build has not been called by another program
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  message("-- Full build allowed")
  enable_testing()
  add_subdirectory(src/test)
endif()

