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
  SIMThermoElasticity(TempSolver& s1, SolidSolver& s2, bool twoway=false)
    : SIMCoupled<TempSolver,SolidSolver>(s1,s2), m_twoway(twoway) {}

  //! \brief Empty destructor.
  virtual ~SIMThermoElasticity() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S2.registerDependency(&this->S1, "temperature1", 1, this->S1.getFEModel(), 1);
    if (m_twoway) {
      this->S1.registerDependency(&this->S2, "displacement1",
                                  this->S1.dimension, this->S2.getFEModel(), 2);
      this->S1.registerDependency(&this->S2, "pressure1", 1,
                                  this->S2.getFEModel(), 2);
    }
  }

  //! \brief Computes initial nodal temperatures for the head solver.
  bool init(const TimeStep&)
  {
    const RealFunc* f = this->S2.getInitialTemperature();
    if (!f)
      return true;
    this->S1.setInitialTemperature(f);

    for (int i = 0; i < this->S1.getNoPatches(); i++) {
      Vector locvec;
      int loc = this->S1.getLocalPatchIndex(i+1);
      if (loc > 0) {
        this->S1.getPatch(loc)->evaluate(f, locvec);
        this->S1.getPatch(loc)->injectNodeVec(locvec, this->S1.getSolution(), 1);
      }
    }

    return true;
  }
protected:
  bool m_twoway;
};

#endif
