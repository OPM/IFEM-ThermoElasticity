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

#ifndef _SIM_HEAT_EQUATION_H_
#define _SIM_HEAT_EQUATION_H_

#include "AnaSol.h"
#include "ASMstruct.h"
#include "DataExporter.h"
#include "ForceIntegrator.h"
#include "Functions.h"
#include "InitialConditionHandler.h"
#include "Profiler.h"
#include "Property.h"
#include "SIMoutput.h"
#include "SIMSolver.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "LinIsotropic.h"
#include <fstream>
#include "HeatQuantities.h"


/*!
  \brief Driver class for a heat equation simulator.
  \details The class incapsulates data and methods for solving a
  heat equation problem using NURBS-based finite elements.
*/

template<class Dim, class Integrand> class SIMHeatEquation : public Dim
{
  //! \brief Struct containing parameters for boundary heat flux calculation.
  struct BoundaryFlux
  {
    std::string file; //!< Name of output file for boundary flux calculation
    std::string set;  //!< Name of topology set for boundary flux calculation
    int code;         //!< Code identifying the boundary for flux calculation
    int timeIncr;     //!< Time level increment for boundary flux calculation
    //! \brief Default constructor.
    BoundaryFlux(int c = 0) : code(c), timeIncr(1) {}
    //! \brief Constructor naming the topology set for force calculation.
    BoundaryFlux(const std::string& s) : set(s), code(0), timeIncr(1) {}
  };

  //! \brief Helper class for searching among BoundaryForce objects.
  class hasCode
  {
    int myCode; //!< The property code to compare with
  public:
    //! \brief Constructor initializing the property code to search for.
    hasCode(int code) : myCode(abs(code)) {}
    //! \brief Returns \e true if the BoundaryForce \a b has the code \a myCode.
    bool operator()(const BoundaryFlux& b) { return abs(b.code) == myCode; }
  };

public:
  struct SetupProps
  {
    bool shareGrid;
    IntegrandBase* integrand;
    SIMoutput* share;

    SetupProps() : shareGrid(false), integrand(NULL), share(NULL) {}
  };

  //! \brief Default constructor.
  //! \param[in] order Order of temporal integration (1 or 2)
  SIMHeatEquation(int order) :
    Dim(1), he(Dim::dimension,order), wdc(Dim::dimension), inputContext("heatequation")
  {
    Dim::myProblem = &he;
    Dim::myHeading = "Heat equation solver";
  }

  //! \brief The destructor zero out the integrand pointer (deleted by parent).
  virtual ~SIMHeatEquation()
  {
    Dim::myProblem = NULL;
    Dim::myInts.clear();
  }

  using Dim::parse;
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
        mVec.push_back(new typename Integrand::MaterialType);
        mVec.back()->parse(child);
      }
      else if (!strcasecmp(child->Value(),"heatflux") ||
               !strcasecmp(child->Value(),"storedenergy")) {
        BoundaryFlux flux;
        utl::getAttribute(child,"set",flux.set);
        utl::getAttribute(child,"file",flux.file);
        utl::getAttribute(child,"stride",flux.timeIncr);
        if (flux.set.empty())
          utl::getAttribute(child,"code",flux.code);
        if (!flux.set.empty()) {
          size_t oldcode = strcasecmp(child->Value(),"heatflux") ? senergy.size() :
                                                                   fluxes.size();
          flux.code = this->getUniquePropertyCode(flux.set,(oldcode+1)*1000);
        }
        strcasecmp(child->Value(),"heatflux") ? senergy.push_back(flux) :
                                                fluxes.push_back(flux);
      } else if (!strcasecmp(child->Value(),"environmentproperties")) {
          double T=273.5, alpha=1.0;
          utl::getAttribute(child,"T",T);
          utl::getAttribute(child,"alpha",alpha);
          wdc.setEnvTemperature(T);
          wdc.setEnvConductivity(alpha);
      } else

        this->Dim::parse(child);
    }

    if (!mVec.empty()) {
      wdc.setMaterial(mVec.front());
      he.setMaterial(mVec.front());
    }

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

    if (!this->solveSystem(temperature.front(),Dim::msgLevel-1,"temperature "))
      return false;

    if (Dim::msgLevel == 1)
    {
      size_t iMax[1];
      double dMax[1];
      double normL2 = this->solutionNorms(temperature.front(),dMax,iMax,1);
      if (Dim::myPid == 0)
        std::cout <<"  Temperature summary: L2-norm         : "<< normL2
                  <<"\n                       Max temperature : "<< dMax[0]
                  << std::endl;
    }

    return true;
  }

  bool postSolve(const TimeStep& tp,bool) {return true;}

  //! \brief Compute and save a boundary heat flux or the stored energy in a volume
  //! \param[in] bf Description of integration domain
  //! \param[in] tp Time stepping information
  //! \param[in] flux True to calculate a heat flux, false for stored energy
  bool saveIntegral(const BoundaryFlux& bf, const TimeStep& tp, bool flux)
  {
    if (bf.code == 0 || bf.timeIncr < 1 || bf.set.empty()) return true;
    if (tp.step < 1 || (tp.step-1)%bf.timeIncr > 0) return true;

    Vector integral;

    if (flux)
      integral = SIM::getBoundaryForce(temperature,this,bf.code,tp.time);
    else {
      HeatEquationStoredEnergy<Integrand> energy(he);
      energy.initBuffer(this->getNoElms());
      SIM::integrate(temperature,this,bf.code,tp.time,&energy);
      energy.assemble(integral);
    }

    if (integral.size() == 0)
      return false;

    std::ostream* os = &std::cout;
    std::stringstream str;

    if (Dim::myPid == 0) {
      if (bf.file.empty())
        std::cout << std::endl;
      else
        os = new std::ofstream(bf.file.c_str(),
                               tp.step == 1 ? std::ios::out : std::ios::app);

      char line[256];
      if (tp.step == 1) {
        if (flux) {
          *os <<"# Heat flux over surface with code "<< bf.code << std::endl;
          sprintf(line,"#%9s %11s\n", "time","Flux");
        } else {
          *os <<"# Stored energy in volumes with code "<< bf.code << std::endl;
          sprintf(line,"#%9s %11s\n", "time","Energy");
        }
        sprintf(line,"#%9s %11s\n", "time","Flux");
        *os << line;
      }

      sprintf(line,"%10.6f", tp.time.t);
      str << line;
      sprintf(line," %11.6g", integral[0]);
      str << line;
      str << '\n';
    }

    if (bf.file.empty())
      utl::printSyncronized(std::cout, str, Dim::myPid);
    else if (Dim::myPid == 0)
    {
      *os << str.str();
      delete os;
    }

    return true;
  }


  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMHeatEquation::saveStep");

    for (size_t i = 0; i < fluxes.size(); ++i)
      saveIntegral(fluxes[i],tp,true);

    for (size_t i = 0; i < senergy.size(); ++i)
      saveIntegral(senergy[i],tp,false);

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
  void setCommunicator(const MPI_Comm* comm) { Dim::adm.setCommunicator(comm); }
#endif

  //! \brief Set the function for the initial temperature field
  //! \param f The function
  void setInitialTemperature(const RealFunc* f) { he.setInitialTemperature(f); }

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to ensure that threading groups are
  //! established for the patch faces subjected to boundary flux integration.
  virtual bool preprocessB()
  {
    PropertyVec::const_iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); p++)
      if (std::find_if(fluxes.begin(),fluxes.end(), hasCode(p->pindx)) != fluxes.end())
        this->generateThreadGroups(*p,SIMinput::msgLevel < 2);

    return true;
  }
  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    he.setMaterial(mVec[propInd]);
    wdc.setMaterial(mVec[propInd]);
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
    wdc.setFlux(tit->second);
    return true;
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the generic Neumann property codes.
  virtual void preprocessA()
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    // Couple the weak Dirichlet integrand to the generic Neumann property codes
    PropertyVec::iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); p++)
      if (p->pcode == Property::NEUMANN_GENERIC ||
          p->pcode == Property::ROBIN)
      {
        if (Dim::myInts.find(p->pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p->pindx,&wdc));
      }
  }
private:
  Integrand he;                 //!< Integrand
  typename Integrand::WeakDirichlet wdc; //!< Weak dirichlet integrand
  std::vector<typename Integrand::MaterialType*> mVec;  //!< Material data

  Vectors temperature;      //!< Temperature solution vectors
  std::string inputContext; //!< Input context

  std::vector<BoundaryFlux> fluxes;  //!< Heat fluxes to calculate
  std::vector<BoundaryFlux> senergy; //!< Stored energies to calculate
};


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
struct SolverConfigurator< SIMHeatEquation<Dim,Integrand> > {
  int setup(SIMHeatEquation<Dim,Integrand>& ad,
            const typename SIMHeatEquation<Dim,Integrand>::SetupProps& props, char* infile)
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
