// $Id$
//==============================================================================
//!
//! \file SIMThermoElasticity.h
//!
//! \date Aug 05 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for thermo-elastic simulators.
//!
//==============================================================================

#ifndef _SIM_THERMO_ELASTICITY_H_
#define _SIM_THERMO_ELASTICITY_H_

#include "SIMCoupled.h"


/*!
  \brief Driver class for thermo-elastic simulators.
  \details A thermo-elastic simulator is a coupling between a thermal solver
  and an elasticity solver.
*/

template<class TempSolver, class SolidSolver>
class SIMThermoElasticity : public SIMCoupled<TempSolver,SolidSolver>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMThermoElasticity(TempSolver& s1, SolidSolver& s2)
    : SIMCoupled<TempSolver,SolidSolver>(s1,s2) {}

  //! \brief Empty destructor.
  virtual ~SIMThermoElasticity() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S2.registerDependency(&this->S1, "temperature1", 1);
  }

  //! \brief Computes initial nodal temperatures for the head solver.
  bool init(const TimeStep&)
  {
    const RealFunc* f = this->S2.getInitialTemperature();
    if (!f)
      return true;

    for (int i = 0; i < this->S1.getNoPatches(); i++) {
      Vector locvec;
      this->S1.getPatch(i+1)->evaluate(f, locvec);
      this->S1.injectPatchSolution(this->S1.getSolution(), locvec, i, 1);
    }

    return true;
  }
};

#endif
