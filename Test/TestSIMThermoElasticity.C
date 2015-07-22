//==============================================================================
//!
//! \file TestSIMThermoElasticity.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for wrapper equpping the linear elasticity solver with dummy
//!        time-stepping support and temperature coupling.
//==============================================================================

#include "SIMThermoElasticity.h"
#include "SIM2D.h"

#include "gtest/gtest.h"

TEST(TestSIMThermoElasticity, Parse)
{
  SIMThermoElasticity<SIM2D> sim;
  EXPECT_TRUE(sim.read("Square.xinp"));
}
