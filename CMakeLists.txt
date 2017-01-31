project(ThermoElasticity)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules
                      $ENV{HOME}/cmake/Modules)

# Required packages
if (NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif(NOT IFEM_CONFIGURED)

set(ELASTICITY_DIR ${PROJECT_SOURCE_DIR}/../IFEM-Elasticity)
if(NOT EXISTS ${ELASTICITY_DIR})
  set(ELASTICITY_DIR ${PROJECT_SOURCE_DIR}/../Elasticity)
endif()

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR}
                    ${ELASTICITY_DIR} ThermoElastic)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT WIN32)
  # Emit position-independent code, suitable for dynamic linking
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  # Enable all warnings
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif(NOT WIN32)

# Testing
enable_testing()
include(IFEMTesting)

if(NOT TARGET Elasticity)
  add_subdirectory(${ELASTICITY_DIR} Elasticity)
endif()

if(NOT TARGET ThermoElastic)
  add_subdirectory(ThermoElastic)
endif()

add_executable(ThermoElasticity main_ThermoElasticity.C)
add_executable(HeatEquation main_HeatEquation.C)

target_link_libraries(ThermoElasticity ThermoElastic Elasticity
                                       IFEMAppCommon ${IFEM_LIBRARIES})
target_link_libraries(HeatEquation ThermoElastic Elasticity
                                   IFEMAppCommon ${IFEM_LIBRARIES})

# Installation
install(TARGETS HeatEquation ThermoElasticity DESTINATION bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} \\
                      ${PROJECT_BINARY_DIR} \\
                      ${PROJECT_SOURCE_DIR}/ThermoElastic \\
                      ${PROJECT_SOURCE_DIR}/../IFEM-Elasticity")
add_doc_target(ThermoElasticity ThermoElasticity)

# Generate regtest script with correct paths
configure_file(${IFEM_REGTEST_SCRIPT} regtest.sh)

if(MPI_FOUND)
  # Add parallel tests here
else()
  ifem_add_test(Annulus.reg ThermoElasticity)
  ifem_add_test(Bar.reg ThermoElasticity)
  ifem_add_test(Square.reg ThermoElasticity)
  ifem_add_test(Square-robin.reg ThermoElasticity)
  ifem_add_test(Annulus-heat.reg HeatEquation)
  ifem_add_test(Annulus-heat-be.reg HeatEquation)
  ifem_add_test(Square-heat.reg HeatEquation)
endif()
list(APPEND TEST_APPS HeatEquation ThermoElasticity)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
