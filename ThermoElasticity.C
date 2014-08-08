// $Id$
//==============================================================================
//!
//! \file ThermoElasticity.C
//!
//! \date Aug 07 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for linear thermo-elasticity problems.
//!
//==============================================================================

#include "ThermoElasticity.h"
#include "FiniteElement.h"
#include "MaterialBase.h"
#include "ElmMats.h"
#include "Tensor.h"
#include "Functions.h"
#include "Utilities.h"
#include "tinyxml.h"


ThermoElasticity::ThermoElasticity (unsigned short int n, bool axS)
  : LinearElasticity(n,axS), myTemp0(NULL)
{
  this->registerVector("temperature1",&myTemp);
}


bool ThermoElasticity::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"initialtemperature"))
  {
    std::string type;
    utl::getAttribute(elem,"type",type,true);
    const TiXmlNode* tval = elem->FirstChild();
    if (tval)
    {
      std::cout <<"\tInitial temperature";
      myTemp0 = utl::parseRealFunc(tval->Value(),type);
      std::cout << std::endl;
    }
  }

  return true;
}


bool ThermoElasticity::initElement (const std::vector<int>& MNPC,
                                    LocalIntegral& elmInt)
{
  elmInt.vec.resize(2);

  // Extract displacement vector for this element
  int ierr = 0;
  if (!primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),elmInt.vec.front());

  // Extract temperature vector for this element
  if (!myTemp.empty() && ierr == 0)
    ierr = utl::gather(MNPC,1,myTemp,elmInt.vec.back());

  if (ierr == 0) return true;

  std::cerr <<" *** ThermoElasticity::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


bool ThermoElasticity::formInitStrainForces (ElmMats& elMat, const Vector& N,
                                             const Matrix& B, const Matrix& C,
                                             const Vec3& X, double detJW) const
{
  if (!eS || !myTemp0 || elMat.vec.size() < 2)
    return true; // No temperature field

  SymmTensor eps(nsd,axiSymmetry), sigma(nsd,axiSymmetry);

  // Strains due to thermal expansion
  double T0 = (*myTemp0)(X);
  double Th = elMat.vec.back().dot(N);
  double alpha = material->getThermalExpansion(Th);
  eps = alpha*(Th-T0)*detJW;

  // Stresses due to thermal expansion
  Vector sigma0;
  if (!C.multiply(eps,sigma0))
    return false;

  // Integrate external forces due to thermal expansion
  sigma = sigma0;
  return B.multiply(sigma,elMat.b[eS-1],true,true); // ES += B^T*sigma0
}


bool ThermoElasticity::evalSol (Vector& s,
                                const FiniteElement& fe, const Vec3& X,
                                const std::vector<int>& MNPC) const
{
  // Extract element displacement and temperatures
  Vector eV, eT;
  int ierr = 0;
  if (!primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),eV);

  // Extract temperature vector for this element
  if (!myTemp.empty() && ierr == 0)
    ierr = utl::gather(MNPC,1,myTemp,eT);

  if (ierr > 0)
  {
    std::cerr <<" *** ThermoElasticity::evalSol: Detected "<< ierr
              <<" node numbers out of range."<< std::endl;
      return false;
  }

  // Evaluate the deformation gradient, dUdX, and/or the strain tensor, eps
  Matrix Bmat;
  Tensor dUdX(nDF);
  SymmTensor eps(nsd,axiSymmetry);
  if (!this->kinematics(eV,fe.N,fe.dNdX,X.x,Bmat,dUdX,eps))
    return false;

  // Add strains due to temperature expansion
  if (myTemp0 && !eT.empty())
  {
    double T0 = (*myTemp0)(X);
    double Th = eT.dot(fe.N);
    eps += material->getThermalExpansion(Th)*(T0-Th);
  }

  // Calculate the stress tensor through the constitutive relation
  Matrix Cmat;
  SymmTensor sigma(nsd, axiSymmetry || material->isPlaneStrain()); double U;
  if (!material->evaluate(Cmat,sigma,U,fe,X,dUdX,eps))
    return false;

  // Congruence transformation to local coordinate system at current point
  if (locSys) sigma.transform(locSys->getTmat(X));

  s = sigma;
  s.push_back(sigma.vonMises());

  // Find the maximum values for each quantity. This block must be performed
  // serially on multi-threaded runs too, due to the update of the maxVal array
  // which is a member of the Elasticity class. Therefore the critical pragma.
#pragma omp critical
  for (size_t j = 0; j < s.size() && j < maxVal.size(); j++)
    if (fabs(s[j]) > fabs(maxVal[j].second))
      maxVal[j] = std::make_pair(X,s[j]);

  return true;
}
