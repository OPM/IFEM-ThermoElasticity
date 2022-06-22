// $Id$
//==============================================================================
//!
//! \file SIMThermoElasticity.h
//!
//! \date Aug 05 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Wrapper equipping the elasticity solver with temperature coupling.
//!
//==============================================================================

#ifndef _SIM_THERMO_ELASTICITY_H_
#define _SIM_THERMO_ELASTICITY_H_

#include "SIMElasticityWrap.h"
#include "SIMconfigure.h"
#include "ElasticityUtils.h"
#include "ThermoElasticity.h"
#include "Linear/AnalyticSolutions.h"
#include "ASMstruct.h"
#include "IFEM.h"
#include "Profiler.h"
#include "TimeStep.h"
#include "Utilities.h"

#include "tinyxml.h"


/*!
  \brief Driver wrapping the elasticity solver with temperature coupling.
*/

template<class Dim> class SIMThermoElasticity : public SIMElasticityWrap<Dim>
{
public:
  typedef bool SetupProps; //!< Dummy declaration (no setup properties required)

  //! \brief Default constructor.
  SIMThermoElasticity();

  //! \brief The destructor clears the VTF-file pointer.
  virtual ~SIMThermoElasticity();

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp);

  //! \brief Returns the initial temperature field.
  const RealFunc* getInitialTemperature() const;

protected:
  using SIMElasticityWrap<Dim>::parse;
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  using SIMElasticityWrap<Dim>::parseAnaSol;
  //! \brief Parses the analytical solution from an XML element.
  virtual bool parseAnaSol(const TiXmlElement* elem);

  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand();

private:
  double startT; //!< Start time for the elasticity solver
};


/*!
  \brief Configuration for a SIMThermoElasticity instance.
*/

template<class Dim> struct SolverConfigurator< SIMThermoElasticity<Dim> >
{
  //! \brief Configures a SIMThermoElasticity instance.
  //! \param elasim The simulator to configure
  //! \param infile The input file to parse
  int setup(SIMThermoElasticity<Dim>& elasim, const bool&, char* infile);
};

#endif
