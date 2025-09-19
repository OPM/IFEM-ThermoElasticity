//==============================================================================
//!
//! \file TestSIMHeatEquation.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for solution driver for the Heat equation.
//!
//==============================================================================

#include "SIM2D.h"
#include "SIMHeatEquation.h"
#include "SIMThermalCoupling.h"
#include "HeatEquation.h"

#include <catch2/catch_test_macros.hpp>


TEST_CASE("TestSIMThermalCoupling.Dependencies")
{
  using Heat2D = SIMHeatEquation<SIM2D,HeatEquation>;
  Heat2D sim(1);
  sim.initSol();
  Vector dummy;
  sim.registerField("foobar", dummy);
  SIMThermalCoupling<Heat2D,Heat2D> couple1(sim, sim, {CouplingDef("foobar", 2, 1)});
  couple1.setupDependencies();

  REQUIRE(sim.getDependentField("temperature1") != nullptr);
  REQUIRE(sim.getDependentField("foobar") != nullptr);
}
