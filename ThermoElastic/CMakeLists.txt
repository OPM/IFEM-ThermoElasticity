project(ThermoElastic)

cmake_minimum_required(VERSION 2.6)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR})

# Common ThermoElastic sources
add_library(ThermoElastic STATIC HeatEquation.C
                                 ThermoElasticity.C
                                 ../../Elasticity/Linear/AnalyticSolutions.C)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/../Test/*.C
                  ${PROJECT_SOURCE_DIR}/../Test
                  ThermoElasticity
                  ThermoElastic Elasticity
                  IFEMAppCommon ${IFEM_LIBRARIES})
set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
