cmake_minimum_required(VERSION 3.22)

project(ThermoElasticity LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

ifem_add_dependency_dir(APP Elasticity TARGET Elasticity)

if(NOT TARGET ThermoElastic)
  add_subdirectory(ThermoElastic)
endif()

ifem_add_application(
  NAME
    ThermoElasticity
  SOURCES
    main_ThermoElasticity.C
  LIBRARIES
    ThermoElastic
)

ifem_add_application(
  NAME
    HeatEquation
  SOURCES
    main_HeatEquation.C
  LIBRARIES
    ThermoElastic
)

# Testing
enable_testing()

# Heat equation
ifem_add_hdf5_test(
  TARGET
    HeatEquation
  TEST_FILES
    Square-heat.hreg
)

ifem_add_vtf_test(
  TARGET
    HeatEquation
  TEST_FILES
    Square-heat.vreg
)

ifem_add_restart_test(
  TARGET
    HeatEquation
  TEST_FILES
    Square-heat-be-restart.reg
    Square-heat-bdf2-restart.reg
  RESTART_LEVEL
    5
)

ifem_add_regression_test(
  TARGET
    HeatEquation
  TEST_FILES
    Annulus-heat.reg
    Annulus-heat-be.reg
    Square-heat.reg
)

# Thermo elasticity
ifem_add_hdf5_test(
  TARGET
    ThermoElasticity
  TEST_FILES
    Annulus.hreg
)

ifem_add_vtf_test(
  TARGET
    ThermoElasticity
  TEST_FILES
    Square.vreg
)

ifem_add_restart_test(
  TARGET
    ThermoElasticity
  TEST_FILES
    Annulus-restart.reg
  RESTART_LEVEL
    5
)

ifem_add_regression_test(
  TARGET
    ThermoElasticity
  TEST_FILES
    Annulus.reg
    Bar.reg
    Square.reg
    Square-robin.reg
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(
  TARGET
    ThermoElasticity
  EXTRA_DIRS
    ${Elasticity_DIRECTORY}
)

# Installation
include(GNUInstallDirs)
install(TARGETS HeatEquation ThermoElasticity DESTINATION ${CMAKE_INSTALL_BINDIR})
