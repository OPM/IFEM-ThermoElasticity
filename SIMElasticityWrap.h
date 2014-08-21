// $Id$
//==============================================================================
//!
//! \file SIMElasticityWrap.h
//!
//! \date Aug 05 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Wrapper equpping the linear elasticity solver with dummy
//! time-stepping support and temperature coupling.
//!
//==============================================================================

#ifndef _SIM_ELASTICITY_WRAP_H_
#define _SIM_ELASTICITY_WRAP_H_

#include "SIMLinEl.h"
#include "ThermoElasticity.h"
#include "DataExporter.h"
#include "TimeStep.h"


/*!
  \brief Driver wrapping the linear elasticity solver with an ISolver interface.
*/

template<class Dim> class SIMElasticityWrap : public SIMLinEl<Dim>
{
public:
  //! \brief Setup properties
  struct SetupProps
  {
    bool checkRHS; //!< If true, check for a that model is right-hand oriented.

    //! \brief Default constructor
    SetupProps() : checkRHS(false) {}
  };

  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMElasticityWrap(bool checkRHS = false) : SIMLinEl<Dim>(checkRHS)
  {
    this->myHeading = "Elasticity solver";
  }
  //! \brief Destructor.
  virtual ~SIMElasticityWrap() { this->setVTF(NULL); }

  //! \brief Registers fields for output to a data exporter.
  virtual void registerFields(DataExporter& exporter)
  {
    exporter.registerField("solid displacement", "solid displacement",
                           DataExporter::SIM, DataExporter::PRIMARY |
                                              DataExporter::SECONDARY);
    exporter.setFieldValue("solid displacement", this, &sol);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  virtual bool saveModel(char* fileName, int& geoBlk, int& nBlock) { return true; }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    return this->writeGlvS(sol,iDump,nBlock);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep&) { return true; }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep&)
  {
    return this->assembleSystem() && this->solveSystem(sol,1);
  }

  bool postSolve(const TimeStep&, bool) { return true; }

  const RealFunc* getInitialTemperature() const
  {
    if (!Dim::myProblem)
      return NULL;

    return static_cast<ThermoElasticity*>(Dim::myProblem)->getInitialTemperature();
  }

  void init()
  {
    sol.resize(this->getNoDOFs(), true);
  }
protected:
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"elasticity") &&
        strcasecmp(elem->Value(),"thermoelasticity"))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      this->getIntegrand()->parse(child);

    return this->SIMLinEl<Dim>::parse(elem);
  }

  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
    {
      if (Dim::dimension == 2)
        Dim::myProblem = new ThermoElasticity(2,SIMLinEl2D::axiSymmetry);
      else
        Dim::myProblem = new ThermoElasticity(Dim::dimension);
    }
    return dynamic_cast<Elasticity*>(Dim::myProblem);
  }

private:
  Vector sol;
};


//! \brief Partial specialization for configurator
  template<class Dim>
struct SolverConfigurator< SIMElasticityWrap<Dim> > {
  int setup(SIMElasticityWrap<Dim>& ad,
            const typename SIMElasticityWrap<Dim>::SetupProps& props, char* infile)
  {
    utl::profiler->start("Model input");

    // Reset the global element and node numbers
    ASMstruct::resetNumbering();
    if (!ad.read(infile))
      return 2;

    utl::profiler->stop("Model input");

    // Preprocess the model and establish data structures for the algebraic system
    if (!ad.preprocess())
      return 3;

    // Initialize the linear solvers
    ad.setMode(SIM::STATIC);
    ad.initSystem(ad.opt.solver,1,1);
    ad.setAssociatedRHS(0,0);
    ad.setQuadratureRule(ad.opt.nGauss[0]);
    ad.init();

    return 0;
  }
};

#endif
