// $Id$
//==============================================================================
//!
//! \file SIMThermoElasticity .h
//!
//! \date Aug 05 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Wrapper equipping the linear elasticity solver with
//! dummy time-stepping support and temperature coupling.
//!
//==============================================================================

#ifndef _SIM_THERMO_ELASTICITY_H_
#define _SIM_THERMO_ELASTICITY_H_

#include "SIMElasticity.h"
#include "SIMSolver.h"
#include "ThermoElasticity.h"
#include "Linear/AnalyticSolutions.h"
#include "ASMstruct.h"
#include "DataExporter.h"
#include "Profiler.h"


/*!
  \brief Driver wrapping the linear elasticity solver with an ISolver interface.
*/

template<class Dim> class SIMThermoElasticity : public SIMElasticity<Dim>
{
public:
  typedef bool SetupProps; //!< Dummy declaration (no setup properties required)

  //! \brief Default constructor.
  SIMThermoElasticity(bool checkRHS = false) : SIMElasticity<Dim>(checkRHS)
  {
    Dim::myHeading = "Thermo-Elasticity solver";
    Dim::msgLevel = 1; // prints the solution summary only
    startT = 0.0;
  }

  //! \brief The destructor clears the VTF-file pointer.
  virtual ~SIMThermoElasticity() { this->setVTF(nullptr); }

  //! \brief Registers fields for output to a data exporter.
  void registerFields(DataExporter& exporter)
  {
    int results = DataExporter::PRIMARY;
    if (!Dim::opt.pSolOnly)
      results |= DataExporter::SECONDARY;
    exporter.registerField("solid displacement", "solid displacement",
                           DataExporter::SIM, results);
    exporter.setFieldValue("solid displacement", this, &sol);
  }

  //! \brief Initializes the solution vector.
  void initSol() { sol.resize(this->getNoDOFs(),true); }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (tp.time.t+0.001*tp.time.dt < startT)
      return true;

    PROFILE1("SIMThermoElasticity::saveStep");

    double old_tol = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    bool ok = this->savePoints(sol,tp.time.t,tp.step);
    utl::zero_print_tol = old_tol;

    if (Dim::opt.format < 0 || !ok)
      return ok;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    return this->writeGlvS(sol,iDump,nBlock);
  }

  //! \brief Dummy method.
  bool init(const TimeStep&) { return true; }
  //! \brief Dummy method.
  bool advanceStep(TimeStep&) { return true; }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (tp.time.t+0.001*tp.time.dt < startT)
      return true;

    PROFILE1("SIMThermoElasticity::solveStep");

    this->setMode(SIM::STATIC);
    this->setQuadratureRule(Dim::opt.nGauss[0]);
    if (!this->assembleSystem()) return false;
    if (!this->solveSystem(sol,1)) return false;

    return this->postSolve(tp);
  }

  //! \brief Postprocesses the solution of current time step.
  bool postSolve(const TimeStep& tp, bool = false)
  {
    Vectors gNorm;
    this->setMode(SIM::RECOVERY);
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,Vectors(1,sol),gNorm))
      return false;
    else if (gNorm.empty())
      return true;

    IFEM::cout <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< gNorm[0](1);
    if (gNorm[0](2) != 0.0)
      IFEM::cout <<"\nExternal energy ((f,u^h)+(t,u^h)^0.5 : "<< gNorm[0](2);
    if (this->haveAnaSol() && gNorm[0].size() >= 4)
      IFEM::cout <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< gNorm[0](3)
                 <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< gNorm[0](4)
                 <<"\nExact relative error (%) : "
                 << gNorm[0](4)/gNorm[0](3)*100.0;
    IFEM::cout << std::endl;
    return true;
  }

  //! \brief Returns the initial temperature field.
  const RealFunc* getInitialTemperature() const
  {
    ThermoElasticity* thelp = dynamic_cast<ThermoElasticity*>(Dim::myProblem);
    return thelp ? thelp->getInitialTemperature() : nullptr;
  }

protected:
  using SIMElasticity<Dim>::parse;
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"elasticity") &&
        strcasecmp(elem->Value(),"thermoelasticity"))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"start"))
        utl::getAttribute(child,"time",startT);

      else if (!strcasecmp(child->Value(),"anasol"))
      {
        std::string type;
        utl::getAttribute(child,"type",type,true);
        if (type == "pipe")
        {
          double Ri = 0.0, Ro = 0.0, Ti = 0.0, To = 0.0, T0 = 273.0;
          double E = 0.0, nu = 0.0, alpha = 0.0;
          bool polar = false;
          utl::getAttribute(child,"Ri",Ri);
          utl::getAttribute(child,"Ro",Ro);
          utl::getAttribute(child,"Ti",Ti);
          utl::getAttribute(child,"To",To);
          utl::getAttribute(child,"Tref",T0);
          utl::getAttribute(child,"E",E);
          utl::getAttribute(child,"nu",nu);
          utl::getAttribute(child,"alpha",alpha);
          utl::getAttribute(child,"polar",polar);
          IFEM::cout <<"\tAnalytical solution: Pipe Ri="<< Ri <<" Ro="<< Ro
                     <<" Ti="<< Ti <<" To="<< To << std::endl;
          if (!Dim::mySol)
            Dim::mySol = new AnaSol(new Pipe(Ri,Ro,Ti,To,T0,E,nu,alpha,
                                             Dim::dimension == 3, polar));
        }
      }
      else
        this->getIntegrand()->parse(child);

    return this->SIMElasticity<Dim>::parse(elem);
  }

  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
    {
      if (Dim::dimension == 2)
        Dim::myProblem = new ThermoElasticity(2,this->axiSymmetry);
      else
        Dim::myProblem = new ThermoElasticity(Dim::dimension);
    }
    return dynamic_cast<Elasticity*>(Dim::myProblem);
  }

private:
  Vector sol;    //!< Primary solution vector
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
  int setup(SIMThermoElasticity<Dim>& elasim, const bool&, char* infile)
  {
    utl::profiler->start("Model input");

    ASMstruct::resetNumbering();
    if (!elasim.read(infile))
      return 2;

    utl::profiler->stop("Model input");
    elasim.opt.print(IFEM::cout) << std::endl;

    // Preprocess the model and establish FE data structures
    if (!elasim.preprocess())
      return 3;

    // Initialize the linear equation system solver
    elasim.initSystem(elasim.opt.solver);
    elasim.initSol();

    return 0;
  }
};

#endif
