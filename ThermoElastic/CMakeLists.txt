project(ThermoElastic)

cmake_minimum_required(VERSION 3.0)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR})

# Common ThermoElastic sources
set(TE_SOURCES HeatEquation.C
               ThermoElasticity.C)

add_library(ThermoElastic STATIC ${TE_SOURCES}
                                 ${ELASTICITY_DIR}/Linear/AnalyticSolutions.C)
foreach(src ${TE_SOURCES})
  list(APPEND CHECK_SOURCES ${PROJECT_SOURCE_DIR}/${src})
endforeach()
list(APPEND CHECK_SOURCES ${ELASTICITY_DIR}/Linear/AnalyticSolutions.C)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/../Test/*.C
                  ${PROJECT_SOURCE_DIR}/../Test
                  ThermoElasticity
                  ThermoElastic Elasticity
                  IFEMAppCommon ${IFEM_LIBRARIES})
set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
