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

#include "gtest/gtest.h"

TEST(TestSIMHeatEquation, Parse)
{
  SIMHeatEquation<SIM2D,HeatEquation> sim(2);
  EXPECT_TRUE(sim.read("Square.xinp"));

  const HeatEquation& heat = static_cast<const HeatEquation&>(*sim.getProblem());
  const LinIsotropic& mat = static_cast<const LinIsotropic&>(*heat.getMaterial());

  ASSERT_FLOAT_EQ(mat.getThermalExpansion(1.0), 1.2e-7);
  ASSERT_FLOAT_EQ(mat.getHeatCapacity(1.0), 1.0);
  ASSERT_FLOAT_EQ(mat.getThermalConductivity(1.0), 0.1);
}
