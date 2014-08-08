project(ThermoElasticity)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
set(IFEM_PRECONFIGURED ${IFEM_CONFIGURED})
if (NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
endif(NOT IFEM_CONFIGURED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

include_directories(${IFEM_INCLUDES} ../Common ../AdvectionDiffusion
                    ../Elasticity  ../Elasticity/Linear)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif(NOT WIN32)

if (NOT IFEM_PRECONFIGURED)
  add_subdirectory(../Common CommonIFEM.dir)
endif(NOT IFEM_PRECONFIGURED)

add_executable(ThermoElasticity ../AdvectionDiffusion/AdvectionDiffusion.C
                                ../AdvectionDiffusion/AdvectionDiffusionBDF.C
                                ../Elasticity/ElasticBase.C
                                ../Elasticity/Elasticity.C
                                ../Elasticity/LinIsotropic.C
                                ../Elasticity/SIMElasticity.C
                                ../Elasticity/Linear/AnalyticSolutions.C
                                ../Elasticity/Linear/LinearElasticity.C
                                ../Elasticity/Linear/SIMLinEl2D.C
                                ../Elasticity/Linear/SIMLinEl3D.C
                                ThermoElasticity.C
                                main_ThermoElasticity.C)

target_link_libraries(ThermoElasticity ${IFEM_LIBRARIES} CommonIFEM)

# Installation
install(TARGETS ThermoElasticity DESTINATION bin)

# Testing
enable_testing()

# Generate regtest script with correct paths
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/doregtests.sh.in
  doregtests.sh @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/check-patchseries.sh.in
               check-patchseries.sh)

if(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
else(IFEM_USE_PARALLEL_PETSC)
  file(GLOB AD_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                         ${PROJECT_SOURCE_DIR}/Test/*.reg)
  foreach(TESTFILE ${AD_TESTFILES})
    ifem_add_test(${TESTFILE} ThermoElasticity)
  endforeach()
endif(IFEM_USE_PARALLEL_PETSC)
