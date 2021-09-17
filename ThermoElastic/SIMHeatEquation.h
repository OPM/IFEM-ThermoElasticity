// $Id$
//==============================================================================
//!
//! \file SIMHeatEquation.h
//!
//! \date Aug 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for the Heat equation.
//!
//==============================================================================

#ifndef _SIM_HEAT_EQUATION_H_
#define _SIM_HEAT_EQUATION_H_

#include "SIMsolution.h"
#include "SIMoutput.h"
#include "SIMconfigure.h"
#include "LinIsotropic.h"
#include "HeatQuantities.h"
#include "ForceIntegrator.h"
#include "Property.h"
#include "ASMstruct.h"
#include "AnaSol.h"
#include "DataExporter.h"
#include "Functions.h"
#include "Profiler.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <fstream>


/*!
  \brief Driver class for a heat equation simulator.
  \details The class incapsulates data and methods for solving a
  heat equation problem using NURBS-based finite elements.
*/

template<class Dim, class Integrand>
class SIMHeatEquation : public Dim, public SIMsolution
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
  //! \brief Struct containing setup parameters.
  struct SetupProps
  {
    bool shareGrid = false;     //!< \e true to share grid with another simulator
    SIMoutput* share = nullptr; //!< Simulator to share grid with
  };

  //! \brief Default constructor.
  //! \param[in] order Order of temporal integration (1 or 2)
  SIMHeatEquation(int order) :
    Dim(1), heq(Dim::dimension,order), wdc(Dim::dimension)
  {
    Dim::myProblem = &heq;
    Dim::myHeading = "Heat equation solver";
    inputContext = "heatequation";
    srcF = nullptr;
  }

  //! \brief The destructor zero out the integrand pointer (deleted by parent).
  virtual ~SIMHeatEquation()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();

    for (Material* mat : mVec)
      delete mat;

    delete srcF;
  }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override
  {
    if (!strcasecmp(elem->Value(),"thermoelasticity")) {
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        if (!strcasecmp(child->Value(),"isotropic")) {
          int code = this->parseMaterialSet(child,mVec.size());
          IFEM::cout <<"\tMaterial code "<< code <<":";
          mVec.push_back(new LinIsotropic());
          mVec.back()->parse(child);
          heq.setMaterial(mVec.back());
          wdc.setMaterial(mVec.back());
        }
      return true;
    }
    else if (strcasecmp(elem->Value(),inputContext.c_str()))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"anasol")) {
        IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!Dim::mySol)
          Dim::mySol = new AnaSol(child);

        // Define the analytical boundary traction field
        int code = 0;
        if (utl::getAttribute(child,"code",code))
          if (code > 0 && Dim::mySol->getScalarSecSol())
          {
            this->setPropertyType(code,Property::NEUMANN);
            Dim::myVectors[code] = Dim::mySol->getScalarSecSol();
          }
      }

      else if (!strcasecmp(child->Value(),"heatflux") ||
               !strcasecmp(child->Value(),"storedenergy")) {
        BoundaryFlux flux;
        utl::getAttribute(child,"set",flux.set);
        utl::getAttribute(child,"file",flux.file);
        utl::getAttribute(child,"stride",flux.timeIncr);
        if (flux.set.empty())
          utl::getAttribute(child,"code",flux.code);
        else {
          size_t oldcode = strcasecmp(child->Value(),"heatflux") ? senergy.size() :
                                                                   fluxes.size();
          flux.code = this->getUniquePropertyCode(flux.set,(oldcode+1)*1000);
        }
        strcasecmp(child->Value(),"heatflux") ? senergy.push_back(flux) :
                                                fluxes.push_back(flux);
      }

      else if (!strcasecmp(child->Value(),"environmentproperties")) {
        double T = 273.5, alpha = 1.0;
        utl::getAttribute(child,"T",T);
        utl::getAttribute(child,"alpha",alpha);
        wdc.setEnvTemperature(T);
        wdc.setEnvConductivity(alpha);
      }

      else if (!strcasecmp(child->Value(),"source")) {
        std::string type;
        utl::getAttribute(child, "type", type, true);
        if (type == "expression" && child->FirstChild() && !srcF) {
          IFEM::cout <<"\tSource function (expression)";
          srcF = utl::parseRealFunc(child->FirstChild()->Value(),type);
          IFEM::cout << std::endl;
          heq.setSource(srcF);
        }
      }

      else
        this->Dim::parse(child);

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "HeatEquation"; }

  //! \brief Initializes the temperature solution vectors.
  void initSol()
  {
    this->initSystem(Dim::opt.solver);

    size_t n, nSols = this->getNoSolutions();
    std::string str = "temperature1";
    this->initSolution(this->getNoDOFs(),nSols);
    for (n = 0; n < nSols; n++, str[11]++)
      this->registerField(str,solution[n]);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Dummy method.
  bool init(const TimeStep&) { return true; }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&)
  {
    this->pushSolution(); // Update solution vectors between time steps
    heq.advanceStep();
    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMHeatEquation::solveStep");

    if (Dim::msgLevel >= 0)
      IFEM::cout <<"\n  step = "<< tp.step
                 <<"  time = "<< tp.time.t << std::endl;

    Vector dummy;
    this->updateDirichlet(tp.time.t,&dummy);

    this->setMode(SIM::DYNAMIC);
    this->setQuadratureRule(Dim::opt.nGauss[0]);
    if (!this->assembleSystem(tp.time,solution))
      return false;

    if (!this->solveSystem(solution.front(),Dim::msgLevel-1,"temperature "))
      return false;

    if (Dim::msgLevel == 1)
    {
      size_t iMax[1];
      double dMax[1];
      double normL2 = this->solutionNorms(solution.front(),dMax,iMax,1);
      IFEM::cout <<"  Temperature summary: L2-norm         : "<< normL2
                 <<"\n                       Max temperature : "<< dMax[0]
                 << std::endl;
    }

    return true;
  }

  //! \brief Evaluates and prints out solution norms.
  void printFinalNorms(const TimeStep& tp)
  {
    Vectors gNorm;
    this->setMode(SIM::RECOVERY);
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,solution,gNorm))
      return;
    else if (gNorm.empty())
      return;

    IFEM::cout <<"L2 norm |t^h| = a(t^h,t^h)^0.5      : "<< gNorm[0](1);
    IFEM::cout <<"\nH1 norm |t^h| = a(t^h,t^h)^0.5      : "<< gNorm[0](2);
    if (this->haveAnaSol() && gNorm[0].size() >= 7)
      IFEM::cout <<"\nL2 norm |t|   = (t,t)^0.5           : "<< gNorm[0](4)
                 <<"\nH1 norm |t|   = a(t,t)^0.5          : "<< gNorm[0](6)
                 <<"\nL2 norm |e|   = (e,e)^0,5, e=t-t^h  : "<< gNorm[0](5)
                 <<"\nH1 norm |e|   = a(e,e)^0.5, e=t-t^h : "<< gNorm[0](7)
                 <<"\nExact relative error (%)            : "
                 << gNorm[0](7)/gNorm[0](6)*100.0;
    IFEM::cout << std::endl;
  }

  //! \brief Compute and save boundary heat flux or stored energy in a volume
  //! \param[in] bf Description of integration domain
  //! \param[in] tp Time stepping information
  //! \param[in] flux True to calculate a heat flux, false for stored energy
  bool saveIntegral(const BoundaryFlux& bf, const TimeStep& tp, bool flux)
  {
    if (bf.code == 0 || bf.timeIncr < 1 || bf.set.empty()) return true;
    if (tp.step < 1 || (tp.step-1)%bf.timeIncr > 0) return true;

    Vector integral;

    if (flux)
      integral = SIM::getBoundaryForce(solution,this,bf.code,tp.time);
    else {
      HeatEquationStoredEnergy<Integrand> energy(heq);
      energy.initBuffer(this->getNoElms());
      SIM::integrate(solution,this,bf.code,tp.time,&energy);
      energy.assemble(integral);
    }

    if (integral.empty())
      return false;

    std::ostream* os = &std::cout;
    std::ofstream of;

    if (Dim::myPid == 0) {
      if (bf.file.empty())
        std::cout << std::endl;
      else {
        of.open(bf.file.c_str(),
                tp.step == 1 ? std::ios::out : std::ios::app);
        os = &of;
      }

      char line[256];
      if (tp.step == 1) {
        *os << (flux ? "# Heat flux over surface" : "# Stored energy in volume")
            <<" with code "<< bf.code << std::endl;
        sprintf(line,"#%9s %11s\n", "time", flux ? "Flux" : "Energy");
        *os << line;
      }
      sprintf(line,"%10.6f %11.6g\n", tp.time.t, integral[0]);
      *os << line;
    }

    return true;
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    PROFILE1("SIMHeatEquation::saveStep");

    bool ok = true;
    for (const BoundaryFlux& f : fluxes)
      ok &= this->saveIntegral(f,tp,true);
    for (const BoundaryFlux& f : senergy)
      ok &= this->saveIntegral(f,tp,false);

    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    ok &= this->savePoints(solution.front(),tp.time.t,tp.step);
    utl::zero_print_tol = old;

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0 || !ok)
      return ok;

    int iDump = 1 + tp.step/Dim::opt.saveInc;

    // Write solution fields
    if (this->writeGlvS1(solution.front(),iDump,nBlock,
                         tp.time.t,"temperature",89) < 0)
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
  }

  using SIMsolution::getSolution;
  //! \brief Returns a reference to current solution vector.
  Vector& getSolution() { return solution.front(); }

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter, const std::string& prefix="")
  {
    exporter.registerField("theta","temperature",DataExporter::SIM,
                           DataExporter::PRIMARY, prefix);
    exporter.setFieldValue("theta", this, &this->getSolution(0));
  }

  //! \brief Set context to read from input file
  void setContext(int ctx)
  {
    std::stringstream str;
    str << "heatequation-" << ctx;
    inputContext = str.str();
  }

  //! \brief Sets the function of the initial temperature field.
  void setInitialTemperature(const RealFunc* f) { heq.setInitialTemperature(f); }
  //! \brief Returns the function of the initial temperature field.
  const RealFunc* getInitialTemperature() const { return heq.getInitialTemperature(); }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SerializeMap& data) const override
  {
    return this->saveSolution(data,this->getName());
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const SerializeMap& data) override
  {
    if (!this->restoreSolution(data,this->getName()))
      return false;

    heq.advanceStep();
    return true;
  }

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to ensure that threading groups are
  //! established for the patch faces subjected to boundary flux integration.
  bool preprocessB() override
  {
    for (const Property& p : Dim::myProps)
      if (std::find_if(fluxes.begin(),fluxes.end(),hasCode(p.pindx)) != fluxes.end())
        this->generateThreadGroups(p,SIMadmin::msgLevel < 2);

    return true;
  }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  bool initMaterial(size_t propInd) override
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    heq.setMaterial(mVec[propInd]);
    wdc.setMaterial(mVec[propInd]);
    return true;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override
  {
    typename Dim::SclFuncMap::const_iterator tit = Dim::myScalars.find(propInd);
    if (tit == Dim::myScalars.end())
      return false;

    heq.setFlux(tit->second);
    wdc.setFlux(tit->second);
    return true;
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the Robin property codes.
  void preprocessA() override
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    for (const Property& p : Dim::myProps)
      if (p.pcode == Property::ROBIN)
        if (Dim::myInts.find(p.pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p.pindx,&wdc));
  }

  Integrand                   heq;  //!< Main integrand

private:
  typename Integrand::WeakDirichlet wdc; //!< Weak dirichlet integrand
  std::vector<Material*>      mVec; //!< Material data
  RealFunc*                   srcF; //!< Source function

  std::string inputContext; //!< Input context

  std::vector<BoundaryFlux> fluxes;  //!< Heat fluxes to calculate
  std::vector<BoundaryFlux> senergy; //!< Stored energies to calculate
};


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
struct SolverConfigurator< SIMHeatEquation<Dim,Integrand> > {
  //! \brief Setup a heat equation simulator.
  //! \param sim The simulator to set up
  //! \param[in] props Setup properties
  //! \param[in] infile The input file to parse
  int setup(SIMHeatEquation<Dim,Integrand>& sim,
            const typename SIMHeatEquation<Dim,Integrand>::SetupProps& props, char* infile)
  {
    if (props.shareGrid)
      // Let the turbulence solver use the same grid as the velocity solver
      sim.clonePatches(props.share->getFEModel(),props.share->getGlob2LocMap());

    // Read the input file
    ASMstruct::resetNumbering();
    if (!sim.readModel(infile))
      return 2;

    // Preprocess the model and establish FE data structures
    if (!sim.preprocess())
      return 3;

    // Initialize the linear equation system solver
    sim.initSol();

    if (props.shareGrid)
      sim.setVTF(props.share->getVTF());

    return 0;
  }
};

#endif
