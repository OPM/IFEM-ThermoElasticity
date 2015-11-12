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

include_directories(${IFEM_INCLUDES} ../Common ${PROJECT_SOURCE_DIR}
                    ../Elasticity ../Elasticity/Linear)

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

file(GLOB LinEl_SRCS ../Elasticity/*.C
                     ../Elasticity/Linear/LinearElasticity.C
                     ../Elasticity/Linear/AnalyticSolutions.C
                     ../Elasticity/Dynamic/SIMLinEl.C)

set(ThermoEl_SRCS HeatEquation.C
                  ThermoElasticity.C)

add_library(CommonThermo STATIC ${LinEl_SRCS} ${ThermoEl_SRCS})

add_executable(ThermoElasticity main_ThermoElasticity.C)
add_executable(HeatEquation main_HeatEquation.C)

target_link_libraries(ThermoElasticity CommonThermo ${IFEM_LIBRARIES} CommonIFEM)
target_link_libraries(HeatEquation CommonThermo ${IFEM_LIBRARIES} CommonIFEM)

# Installation
install(TARGETS HeatEquation ThermoElasticity DESTINATION bin)

# Testing
enable_testing()

# Generate regtest script with correct paths
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

if(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
else(IFEM_USE_PARALLEL_PETSC)
  ifem_add_test(Annulus.reg ThermoElasticity)
  ifem_add_test(Bar.reg ThermoElasticity)
  ifem_add_test(Square.reg ThermoElasticity)
  ifem_add_test(Square-robin.reg ThermoElasticity)
  ifem_add_test(Annulus-heat.reg HeatEquation)
  ifem_add_test(Annulus-heat-be.reg HeatEquation)
  ifem_add_test(Square-heat.reg HeatEquation)
endif(IFEM_USE_PARALLEL_PETSC)
list(APPEND TEST_APPS HeatEquation ThermoElasticity)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  ThermoElasticity
                  CommonThermo CommonIFEM ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
