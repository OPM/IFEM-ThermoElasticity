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
#include "MaterialBase.h"
#include "ElmMats.h"
#include "Tensor.h"
#include "Function.h"
#include "Utilities.h"


ThermoElasticity::ThermoElasticity (unsigned short int n, bool axS)
  : LinearElasticity(n,axS)
{
  this->registerVector("temperature1",&myTempVec);
}


bool ThermoElasticity::initElement (const std::vector<int>& MNPC,
                                    const FiniteElement&, const Vec3&, size_t,
                                    LocalIntegral& elmInt)
{
  elmInt.vec.resize(2);

  // Extract displacement vector for this element
  int ierr = 0;
  if (!primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),elmInt.vec.front());

  // Extract temperature vector for this element
  if (!myTempVec.empty() && ierr == 0)
    ierr = utl::gather(MNPC,1,myTempVec,elmInt.vec.back());

  if (ierr > 0)
  {
    std::cerr <<" *** ThermoElasticity::initElement: Detected "
              << ierr <<" node numbers out of range."<< std::endl;
    return false;
  }

#if INT_DEBUG > 2
  std::cout <<"Element displacement vector"<< elmInt.vec.front()
            <<"Element temperature vector"<< elmInt.vec.back();
#endif
  return true;
}


double ThermoElasticity::getThermalStrain (const Vector& eT, const Vector& N,
                                           const Vec3& X) const
{
  double T0 = myTemp0 ? (*myTemp0)(X) : 0.0;
  double Th = myTemp ? (*myTemp)(X) : eT.dot(N);
  return material->getThermalExpansion(Th)*(Th-T0);
}


bool ThermoElasticity::formInitStrainForces (ElmMats& elMat, const Vector& N,
                                             const Matrix& B, const Matrix& C,
                                             const Vec3& X, double detJW) const
{
  if (!eS || elMat.vec.size() < 2)
    return true; // No temperature field

  // Strains due to thermal expansion
  SymmTensor eps(nsd,axiSymmetry);
  eps = this->getThermalStrain(elMat.vec.back(),N,X)*detJW;

  // Stresses due to thermal expansion
  RealArray sigma0;
  if (!C.multiply(eps,sigma0))
    return false;

  // Integrate external forces due to thermal expansion
  SymmTensor sigma(nsd,axiSymmetry); sigma = sigma0;
  return B.multiply(sigma,elMat.b[eS-1],true,true); // ES += B^T*sigma0
}


bool ThermoElasticity::evalSol (Vector& s,
                                const FiniteElement& fe, const Vec3& X,
                                const std::vector<int>& MNPC) const
{
  // Extract element displacement and temperatures
  Vectors eV(2);
  int ierr = 0;
  if (!primsol.empty() && !primsol.front().empty())
    ierr = utl::gather(MNPC,npv,primsol.front(),eV.front());

  // Extract temperature vector for this element
  if (!myTempVec.empty() && ierr == 0)
    ierr = utl::gather(MNPC,1,myTempVec,eV.back());

  if (ierr > 0)
  {
    std::cerr <<" *** ThermoElasticity::evalSol: Detected "<< ierr
              <<" node numbers out of range."<< std::endl;
    return false;
  }

  return this->evalSol2(s,eV,fe,X);
}
