cmake_minimum_required(VERSION 3.5)

project(ThermoElastic)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR})

# Common ThermoElastic sources
set(TE_SOURCES HeatEquation.C
               SIMHeatEquation.C
               SIMThermoElasticity.C
               ThermoElasticity.C)

add_library(ThermoElastic STATIC ${TE_SOURCES}
                                 ${ELASTICITY_DIR}/Linear/AnalyticSolutions.C)
foreach(src ${TE_SOURCES})
  list(APPEND CHECK_SOURCES ${PROJECT_SOURCE_DIR}/${src})
endforeach()
list(APPEND CHECK_SOURCES ${ELASTICITY_DIR}/Linear/AnalyticSolutions.C)

# Unit tests
IFEM_add_test_app(SOURCES ${PROJECT_SOURCE_DIR}/../Test/*.C
                  WORKDIR ${PROJECT_SOURCE_DIR}/../Test
                  NAME ThermoElasticity
                  LIBRARIES ThermoElastic Elasticity
                            IFEMAppCommon ${IFEM_LIBRARIES})
set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
