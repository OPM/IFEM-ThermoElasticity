// $Id$
//==============================================================================
//!
//! \file SIMHeatEquation.h
//!
//! \date Aug 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for the Heat equation
//!
//==============================================================================

#ifndef _SIM_HEATEQUATION_H_
#define _SIM_HEATEQUATION_H_

#include "AnaSol.h"
#include "ASMstruct.h"
#include "DataExporter.h"
#include "Functions.h"
#include "InitialConditionHandler.h"
#include "Profiler.h"
#include "Property.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "HeatEquation.h"
#include "LinIsotropic.h"


/*!
  \brief Driver class for a heat equation simulator.
  \details The class incapsulates data and methods for solving a
  heat equation problem using NURBS-based finite elements.
*/

template<class Dim> class SIMHeatEquation : public Dim
{
public:
  struct SetupProps
  {
    bool shareGrid;
    IntegrandBase* integrand;
    SIMoutput* share;

    SetupProps() : shareGrid(false), integrand(NULL), share(NULL) {}
  };

  //! \brief Default constructor.
  //! \param[in] order Order of temporal integration (1, 2)
  SIMHeatEquation(int order) :
    Dim(1), he(Dim::dimension, order), inputContext("heatequation")
  {
    Dim::myProblem = &he;
    this->myHeading = "Heat equation solver";
  }

  //! \brief Destructor
  virtual ~SIMHeatEquation()
  {
    Dim::myProblem = NULL;
  }

  //! \brief Parses a data section from an input stream (deprecated).
  virtual bool parse(char*, std::istream&) { return false; }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),inputContext.c_str()) &&
        strcasecmp(elem->Value(),"thermoelasticity"))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement()) {
      if (strcasecmp(child->Value(),"anasol") == 0) {
        std::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!Dim::mySol)
          Dim::mySol = new AnaSol(child);

        // Define the analytical boundary traction field
        int code = 0;
        if (utl::getAttribute(child,"code",code)) {
          if (code > 0 && Dim::mySol->getScalarSecSol())
          {
            this->setPropertyType(code,Property::NEUMANN);
            Dim::myVectors[code] = Dim::mySol->getScalarSecSol();
          }
        }
      } else if (!strcasecmp(child->Value(),"isotropic")) {
        int code = this->parseMaterialSet(child,mVec.size());
        std::cout <<"\tMaterial code "<< code <<":";
        mVec.push_back(new LinIsotropic(false,false));
        mVec.back()->parse(child);
      }
      else
        this->Dim::parse(child);
    }

    if (!mVec.empty())
      he.setMaterial(mVec.front());

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "HeatEquation"; }

  void init(const TimeStep& tp)
  {
    // Initialize temperature solution vectors
    size_t n, nSols = this->getNoSolutions();
    temperature.resize(nSols);
    std::string str = "temperature1";
    for (n = 0; n < nSols; n++, str[11]++) {
      temperature[n].resize(this->getNoDOFs(),true);
      this->registerField(str,temperature[n]);
    }
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  virtual bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    // Update temperature vectors between time steps
    const int nNusols = temperature.size();
    for (int n = nNusols-1; n > 0; n--)
      temperature[n] = temperature[n-1];

    he.advanceStep();

    return true;
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMHeatEquation::solveStep");

    if (Dim::myPid == 0 && Dim::msgLevel >= 0)
      std::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    Vector dummy;
    this->updateDirichlet(tp.time.t, &dummy);

    if (!this->assembleSystem(tp.time, temperature))
      return false;

    if (!this->solveSystem(temperature.front(), Dim::msgLevel-1,"temperature "))
      return false;

    if (Dim::msgLevel == 1)
    {
      size_t iMax[1];
      double dMax[1];
      double normL2 = this->solutionNorms(temperature.front(),dMax,iMax,1);
      if (Dim::myPid == 0)
        std::cout <<"Temperature summary: L2-norm        : "<< normL2
                  <<"\n                   Max temperature : "<< dMax[0]
                  << std::endl;
    }

    return true;
  }

  bool postSolve(const TimeStep& tp,bool) {return true;}

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMHeatEquation::saveStep");

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;

    // Write solution fields
    bool result = this->writeGlvS(temperature.front(), iDump, nBlock,
                                  tp.time.t, true, "temperature", 89);

    return result && this->writeGlvStep(iDump, tp.time.t);
  }

  Vector& getSolution(int n=0) { return temperature[n]; }

  const Vector& getSolution(int n=0) const { return temperature[n]; }

  void registerFields(DataExporter& exporter, const std::string& prefix="")
  {
    exporter.registerField("theta","temperature",DataExporter::SIM,
                           DataExporter::PRIMARY|DataExporter::RESTART,
                           prefix);
    exporter.setFieldValue("theta", this, &temperature.front());
  }

  double externalEnergy(const Vectors&) const { return 0.0; }

  //! \brief Sets initial conditions.
  void setInitialConditions() { SIM::setInitialConditions(*this); }

  //! \brief Set context to read from input file
  void setContext(int ctx)
  {
    std::stringstream str;
    str << "heatequation-" << ctx;
    inputContext = str.str();
  }

#ifdef HAS_PETSC
  //! \brief Set MPI communicator for the linear equation solvers
  //! \param comm The communicator to use
  void setCommunicator(const MPI_Comm* comm)
  {
    this->adm.setCommunicator(comm);
  }
#endif

protected:
  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    he.setMaterial(mVec[propInd]);
    return true;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator tit = Dim::myScalars.find(propInd);
    if (tit == Dim::myScalars.end())
      return false;

    he.setFlux(tit->second);
    return true;
  }


private:
  HeatEquation he; //!< Integrand
  std::vector<Material*> mVec; //!< Material data

  Vectors temperature;
  std::string inputContext; //!< Input context
};


//! \brief Partial specialization for configurator
template<class Dim>
struct SolverConfigurator< SIMHeatEquation<Dim> > {
  int setup(SIMHeatEquation<Dim>& ad,
            const typename SIMHeatEquation<Dim>::SetupProps& props, char* infile)
  {
    utl::profiler->start("Model input");

    if (props.shareGrid)
      // Let the turbulence solver use the same grid as the velocity solver
      ad.clonePatches(props.share->getFEModel(), props.share->getGlob2LocMap());

    // Reset the global element and node numbers
    ASMstruct::resetNumbering();
    if (!ad.read(infile))
      return 2;

    utl::profiler->stop("Model input");

    // Preprocess the model and establish data structures for the algebraic system
    if (!ad.preprocess())
      return 3;

    // Initialize the linear solvers
    ad.setMode(SIM::DYNAMIC);
    ad.initSystem(ad.opt.solver,1,1,false);
    ad.setQuadratureRule(ad.opt.nGauss[0]);

    // Time-step loop
    ad.init(TimeStep());
    if (props.share)
      ad.setVTF(props.share->getVTF());
    ad.setInitialConditions();

    return 0;
  }
};

#endif
