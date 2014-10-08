//==============================================================================
//!
//! \file TestSIMElasticityWrap.C
//!
//! \date Oct 7 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for Wrapper equpping the linear elasticity solver with dummy
//! time-stepping support and temperature coupling
//==============================================================================

#include "SIMElasticityWrap.h"
#include "SIMThermoElasticity.h"

#include "gtest/gtest.h"

TEST(TestSIMElasticityWrap, Parse)
{
  SIMElasticityWrap<SIM2D> sim;
  EXPECT_TRUE(sim.read("Square.xinp"));
}
