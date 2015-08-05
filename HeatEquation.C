// $Id$
//==============================================================================
//!
//! \file HeatEquation.C
//!
//! \date Aug 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for the heat equation.
//!
//==============================================================================

#include "HeatEquation.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "MaterialBase.h"
#include "Vec3Oper.h"
#include "HeatQuantities.h"
#include "WeakOperators.h"


HeatEquation::HeatEquation (unsigned short int n, int order) :
  nsd(n), bdf(order), mat(NULL), flux(NULL), init(NULL)
{
  primsol.resize(order+1);
}


bool HeatEquation::evalInt (LocalIntegral& elmInt,
                            const FiniteElement& fe,
                            const TimeDomain& time,
                            const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double theta=0;
  double rhocp=1.0, kappa=1.0;
  for (int t=0;t<bdf.getOrder();++t) {
    double val = fe.N.dot(elMat.vec[t+1]);
    if (t == 0) {
      rhocp = mat?mat->getMassDensity(X)*mat->getHeatCapacity(val):1.0;
      kappa = mat?mat->getThermalConductivity(val):1.0;
    }

    theta += -bdf[1+t]/time.dt*val;
  }

  WeakOperators::Laplacian(elMat.A[0], fe, kappa);
  WeakOperators::Mass(elMat.A[0], fe, rhocp*bdf[0]/time.dt);

  WeakOperators::Source(elMat.b.front(), fe, rhocp*theta);

  return true;
}


bool HeatEquation::evalBou (LocalIntegral& elmInt,
                            const FiniteElement& fe,
                            const Vec3& X, const Vec3& normal) const
{
  if (!flux) {
    std::cerr <<" *** HeatEquation::evalBou: No fluxes."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  if (elMat.b.empty()) {
    std::cerr <<" *** HeatEquation::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Evaluate the Neumann value
  double T = (*flux)(X);
  double val = fe.N.dot(elMat.vec[0]);
  double kappa = mat?mat->getThermalConductivity(val):1.0;

  // Integrate the Neumann value
  WeakOperators::Source(elMat.b.front(), fe, kappa*T);

  return true;
}


const char* HeatEquation::getField1Name (size_t, const char* prefix) const
{
  if (!prefix)
    return "theta";

  static std::string name;
  name = prefix + std::string(" theta");

  return name.c_str();
}


const char* HeatEquation::getField2Name (size_t i,
                                         const char* prefix) const
{
  if (i > 2)
    return 0;

  static const char* s[3] = { "theta_x", "theta_y", "theta_z" };
  if (!prefix)
    return s[i];

  static std::string name;
  name = prefix + std::string(" ") + s[i];

  return name.c_str();
}


LocalIntegral* HeatEquation::WeakDirichlet::getLocalIntegral(size_t nen,
                                                             size_t,
                                                             bool) const
{
  ElmMats* result = new ElmMats;

  result->withLHS = true;
  result->resize(1,1);
  result->A[0].resize(nen,nen);
  result->b[0].resize(nen);

  return result;
}


bool HeatEquation::WeakDirichlet::evalBou(LocalIntegral& elmInt,
                                          const FiniteElement& fe,
                                          const Vec3& X,
                                          const Vec3& normal) const
{
  if (!flux) {
    std::cerr <<" *** HeatEquation::evalBou: No flux function."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the Neumann value
  double q = (*flux)(X);
  double val = fe.N.dot(elMat.vec[0]);
  double kappa = mat?mat->getThermalConductivity(val):1.0;

  WeakOperators::Mass(elMat.A.front(), fe, -envCond);
  WeakOperators::Source(elMat.b.front(), fe, q);

  for (size_t i=1;i<=fe.N.size();++i) {
    for (size_t j=1;j<=fe.N.size();++j) {
      double dT=0;
      for (size_t k=1;k<=nsd;++k)
        dT += fe.dNdX(j,k)*normal[k-1];

      elMat.A.front()(i,j) += (kappa*dT+envT*envCond)*fe.N(i)*fe.detJxW;
    }
  }

  return true;
}


ForceBase* HeatEquation::getForceIntegrand(const Vec3*, AnaSol*) const
{
  return new HeatEquationFlux<HeatEquation>(*const_cast<HeatEquation*>(this));
}
