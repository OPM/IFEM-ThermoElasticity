cmake_minimum_required(VERSION 3.22)

project(ThermoElastic LANGUAGES C CXX)

ifem_add_library(
  NAME
    ThermoElastic
  SOURCES
    HeatEquation.C
    SIMHeatEquation.C
    SIMThermoElasticity.C
    ThermoElasticity.C
    ${Elasticity_DIRECTORY}/Linear/AnalyticSolutions.C
  LIBRARIES
    IFEMAppCommon
    Elasticity
)

target_include_directories(ThermoElastic INTERFACE ${PROJECT_SOURCE_DIR})

# Unit tests
ifem_add_test_app(
  NAME
    ThermoElasticity
  SOURCES
    ../Test/TestSIMHeatEquation.C
    ../Test/TestSIMThermalCoupling.C
    ../Test/TestSIMThermoElasticity.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/../Test
  LIBRARIES
    ThermoElastic
)

set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
