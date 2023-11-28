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

class RealFunc;
class TimeStep;
namespace tinyxml2 { class XMLElement; }


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
  bool saveStep(const TimeStep& tp, int& nBlock) override;

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp) override;

  //! \brief Returns the initial temperature field.
  const RealFunc* getInitialTemperature() const;

protected:
  using SIMElasticityWrap<Dim>::parse;
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  bool parse(const tinyxml2::XMLElement* elem) override;

  using SIMElasticityWrap<Dim>::parseAnaSol;
  //! \brief Parses the analytical solution from an XML element.
  bool parseAnaSol(const tinyxml2::XMLElement* elem) override;

  //! \brief Returns the actual integrand.
  Elasticity* getIntegrand() override;

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
