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
list(APPEND CHECK_SOURCES main_HeatEquation.C main_ThermoElasticity.C)

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
                      ${ELASTICITY_DIR}")
add_doc_target(ThermoElasticity ThermoElasticity)

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
  ifem_add_hdf5_test(Square-heat.hreg HeatEquation)
  ifem_add_hdf5_test(Annulus.hreg ThermoElasticity)
  if(HDF5_FOUND AND CEREAL_FOUND)
    ifem_add_restart_test(Square-heat-be-restart.reg HeatEquation 5)
    ifem_add_restart_test(Square-heat-bdf2-restart.reg HeatEquation 5)
    ifem_add_restart_test(Annulus-restart.reg ThermoElasticity 5)
  endif()

  ifem_add_vtf_test(Square-heat.vreg HeatEquation)
  ifem_add_vtf_test(Square.vreg ThermoElasticity)
endif()
list(APPEND TEST_APPS HeatEquation ThermoElasticity)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
