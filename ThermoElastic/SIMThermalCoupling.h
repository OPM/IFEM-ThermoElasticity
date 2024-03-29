// $Id$
//==============================================================================
//!
//! \file SIMThermalCoupling.h
//!
//! \date Aug 05 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for simulators coupled to a heat equation.
//!
//==============================================================================

#ifndef _SIM_THERMAL_COUPLING_H_
#define _SIM_THERMAL_COUPLING_H_

#include "SIMCoupled.h"

#include "ASMbase.h"
#include "Function.h"


//! \brief Struct describing a coupling.
struct CouplingDef {
  std::string name; //!< Name of field
  int components;   //!< Number of components in field
  int basis;        //!< Basis of field

  //! \brief Default constructor.
  //! \param[in] nname Name of field
  //! \param[in] comp Number of components in field
  //! \param[in] b Basis of field
  CouplingDef(const std::string& nname, int comp, int b) :
    name(nname), components(comp), basis(b) {}
};

typedef std::vector<CouplingDef> Couplings; //!< Convenience declaration


/*!
  \brief Driver class for thermal coupled simulators.
  \details A thermal coupled simulator is a coupling between a thermal solver
  and another solver.
*/

template<class TempSolver, class OtherSolver,
         template<class T1, class T2> class Coupling=SIMCoupled>
class SIMThermalCoupling : public Coupling<TempSolver,OtherSolver>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMThermalCoupling(TempSolver& s1, OtherSolver& s2,
                     const Couplings& back_coupling = Couplings())
    : Coupling<TempSolver,OtherSolver>(s1,s2), m_back_coupling(back_coupling) {}

  //! \brief Empty destructor.
  virtual ~SIMThermalCoupling() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S2.registerDependency(&this->S1, "temperature1", 1,
                                this->S1.getFEModel(), 1);
    for (auto& it : m_back_coupling)
      this->S1.registerDependency(&this->S2, it.name, it.components,
                                  this->S2.getFEModel(), it.basis);
  }

  //! \brief Computes initial nodal temperatures for the head solver.
  virtual bool init(const TimeStep&)
  {
    const RealFunc* f = this->S2.getInitialTemperature();
    if (!f) return false;

    this->S1.setInitialTemperature(f);

    for (int i = 0; i < this->S1.getNoPatches(); i++) {
      int loc = this->S1.getLocalPatchIndex(i+1);
      if (loc > 0) {
        Vector locvec;
        this->S1.getPatch(loc)->evaluate(f,locvec);
        this->S1.getPatch(loc)->injectNodeVec(locvec,this->S1.getSolution(),1);
      }
    }

    return true;
  }

private:
  Couplings m_back_coupling; //!< Couplings from other solver to thermal solver
};

#endif
