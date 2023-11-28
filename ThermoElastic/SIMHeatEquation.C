// $Id$
//==============================================================================
//!
//! \file SIMHeatEquation.C
//!
//! \date Aug 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for the Heat equation.
//!
//==============================================================================

#include "SIMHeatEquation.h"

#include "HeatEquation.h"
#include "HeatQuantities.h"

#include "AnaSol.h"
#include "ASMstruct.h"
#include "DataExporter.h"
#include "ForceIntegrator.h"
#include "Functions.h"
#include "IFEM.h"
#include "LinIsotropic.h"
#include "Profiler.h"
#include "Property.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TimeStep.h"
#include "Utilities.h"

#include <tinyxml2.h>
#include <fstream>


template<class Dim, class Integrand>
SIMHeatEquation<Dim,Integrand>::SIMHeatEquation (int order) :
  Dim(1), heq(Dim::dimension,order), wdc(Dim::dimension)
{
  Dim::myProblem = &heq;
  Dim::myHeading = "Heat equation solver";
  inputContext = "heatequation";
  srcF = nullptr;
}


template<class Dim, class Integrand>
SIMHeatEquation<Dim,Integrand>::~SIMHeatEquation ()
{
  Dim::myProblem = nullptr;
  Dim::myInts.clear();

  for (Material* mat : mVec)
    delete mat;

  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
  if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);

  delete srcF;
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"thermoelasticity")) {
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
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

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
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


template<class Dim, class Integrand>
void SIMHeatEquation<Dim,Integrand>::initSol ()
{
  this->initSystem(Dim::opt.solver);

  size_t n, nSols = this->getNoSolutions();
  std::string str = "temperature1";
  this->initSolution(this->getNoDOFs(),nSols);
  for (n = 0; n < nSols; n++, str[11]++)
    this->registerField(str,solution[n]);
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::
saveModel (char* fileName, int& geoBlk, int& nBlock)
{
  if (Dim::opt.format < 0) return true;

  nBlock = 0;
  return this->writeGlvG(geoBlk,fileName);
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::advanceStep (TimeStep&)
{
  this->pushSolution(); // Update solution vectors between time steps
  heq.advanceStep();
  return true;
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::solveStep (const TimeStep& tp)
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


template<class Dim, class Integrand>
void SIMHeatEquation<Dim,Integrand>::printFinalNorms (const TimeStep& tp)
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


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::
saveIntegral (const BoundaryFlux& bf, const TimeStep& tp, bool flux)
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


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::saveStep (const TimeStep& tp, int& nBlock)
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


template<class Dim, class Integrand>
void SIMHeatEquation<Dim,Integrand>::
registerFields (DataExporter& exporter, const std::string& prefix)
{
  exporter.registerField("theta","temperature",DataExporter::SIM,
                         DataExporter::PRIMARY, prefix);
  exporter.setFieldValue("theta", this, &this->getSolution(0));
}


template<class Dim, class Integrand>
void SIMHeatEquation<Dim,Integrand>::setContext (int ctx)
{
  std::stringstream str;
  str << "heatequation-" << ctx;
  inputContext = str.str();
}


template<class Dim, class Integrand>
void SIMHeatEquation<Dim,Integrand>::
setInitialTemperature (const RealFunc* f)
{
  heq.setInitialTemperature(f);
}


template<class Dim, class Integrand>
const RealFunc*
SIMHeatEquation<Dim,Integrand>::getInitialTemperature () const
{
  return heq.getInitialTemperature();
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::serialize (SerializeMap& data) const
{
  return this->saveSolution(data,this->getName());
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::deSerialize (const SerializeMap& data)
{
  if (!this->restoreSolution(data,this->getName()))
    return false;

  heq.advanceStep();
  return true;
}


template<class Dim, class Integrand>
void SIMHeatEquation<Dim,Integrand>::preprocessA ()
{
  Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

  for (Property& p : Dim::myProps)
    if (p.pcode == Property::ROBIN) {
      if (Dim::myInts.find(p.pindx) == Dim::myInts.end())
        Dim::myInts.insert(std::make_pair(p.pindx,&wdc));
    } else if (p.pcode == Property::DIRICHLET_ANASOL) {
      if (!Dim::mySol->getScalarSol())
        p.pcode = Property::UNDEFINED;
      else if (aCode[0] == abs(p.pindx))
        p.pcode = Property::DIRICHLET_INHOM;
      else if (aCode[0] == 0)
      {
        aCode[0] = abs(p.pindx);
        Dim::myScalars[aCode[0]] = Dim::mySol->getScalarSol();
        p.pcode = Property::DIRICHLET_INHOM;
      }
      else
        p.pcode = Property::UNDEFINED;
    } else if (p.pcode == Property::NEUMANN_ANASOL) {
      if (!Dim::mySol->getScalarSecSol())
        p.pcode = Property::UNDEFINED;
      else if (aCode[1] == p.pindx)
        p.pcode = Property::NEUMANN;
      else if (aCode[1] == 0)
      {
        aCode[1] = p.pindx;
        Dim::myVectors[aCode[1]] = Dim::mySol->getScalarSecSol();
        p.pcode = Property::NEUMANN;
      }
      else
        p.pcode = Property::UNDEFINED;
    }
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::preprocessB ()
{
  for (const Property& p : Dim::myProps)
    if (std::find_if(fluxes.begin(),fluxes.end(),hasCode(p.pindx)) != fluxes.end())
      this->generateThreadGroups(p,SIMadmin::msgLevel < 2);

  return true;
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size())
    propInd = mVec.size()-1;

  heq.setMaterial(mVec[propInd]);
  wdc.setMaterial(mVec[propInd]);
  return true;
}


template<class Dim, class Integrand>
bool SIMHeatEquation<Dim,Integrand>::initNeumann (size_t propInd)
{
  typename Dim::SclFuncMap::const_iterator tit = Dim::myScalars.find(propInd);
  if (tit == Dim::myScalars.end())
    return false;

  heq.setFlux(tit->second);
  wdc.setFlux(tit->second);
  return true;
}


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
int SolverConfigurator<SIMHeatEquation<Dim,Integrand>>::
setup (SIMHeatEquation<Dim,Integrand>& sim,
       const typename SIMHeatEquation<Dim,Integrand>::SetupProps& props,
       char* infile)
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


template class SIMHeatEquation<SIM1D,HeatEquation>;
template struct SolverConfigurator<SIMHeatEquation<SIM1D,HeatEquation>>;
template class SIMHeatEquation<SIM2D,HeatEquation>;
template struct SolverConfigurator<SIMHeatEquation<SIM2D,HeatEquation>>;
template class SIMHeatEquation<SIM3D,HeatEquation>;
template struct SolverConfigurator<SIMHeatEquation<SIM3D,HeatEquation>>;
