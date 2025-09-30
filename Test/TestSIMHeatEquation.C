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

#include "SIMHeatEquation.h"
#include "SIM2D.h"
#include "HeatEquation.h"
#include "LinIsotropic.h"

#include "Catch2Support.h"


TEST_CASE("TestSIMHeatEquation.Parse")
{
  SIMHeatEquation<SIM2D,HeatEquation> sim(2);
  REQUIRE(sim.read("Square.xinp"));

  const HeatEquation& heat = static_cast<const HeatEquation&>(*sim.getProblem());
  const LinIsotropic& mat = static_cast<const LinIsotropic&>(*heat.getMaterial());

  REQUIRE_THAT(mat.getThermalExpansion(1.0), WithinRel(1.2e-7));
  REQUIRE_THAT(mat.getHeatCapacity(1.0), WithinRel(1.0));
  REQUIRE_THAT(mat.getThermalConductivity(1.0), WithinRel(0.1));
}
