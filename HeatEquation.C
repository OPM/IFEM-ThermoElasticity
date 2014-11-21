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

  // loop over test functions (i) and basis functions (j)
  for (size_t i = 1; i <= fe.N.size(); ++i) {
    for (size_t j = 1; j <= fe.N.size(); ++j) {
      double laplace = 0.0;
      for (size_t k = 1;k <= nsd; ++k)
        laplace += fe.dNdX(i,k)*fe.dNdX(j,k);

      elMat.A[0](i,j) += (kappa*laplace+rhocp*bdf[0]/time.dt*fe.N(i)*fe.N(j))*fe.detJxW;
    }
  }

  elMat.b.front().add(fe.N, rhocp*theta*fe.detJxW);

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
  elMat.b.front().add(fe.N,kappa*T*fe.detJxW);

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


ForceBase* HeatEquation::getForceIntegrand(const Vec3*, AnaSol*) const
{
  return new HeatEquationFlux(*const_cast<HeatEquation*>(this));
}


bool HeatEquationFlux::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                                const TimeDomain& time,
                                const Vec3& X, const Vec3& normal) const
{
  HeatEquation& problem = static_cast<HeatEquation&>(myProblem);
  ElmNorm& elmNorm = static_cast<ElmNorm&>(elmInt);

  const Material* mat = problem.getMaterial();

  double theta = fe.N.dot(elmNorm.vec[0]);
  double kappa=mat?mat->getThermalConductivity(theta):1.0;

  // Temperature gradient in integration point
  Vec3 gradT;
  for (size_t i = 1;i <= fe.N.size();i++)
    for (size_t l = 1;l <= problem.getNoSpaceDim();l++)
      gradT[l-1] += fe.dNdX(i,l)*elmNorm.vec[0](i);

  elmNorm[0] += kappa*gradT*normal*fe.detJxW;

  return true;
}


bool HeatEquationStoredEnergy::evalInt (LocalIntegral& elmInt,
                                        const FiniteElement& fe,
                                        const TimeDomain& time,
                                        const Vec3& X) const
{
  HeatEquation& problem = static_cast<HeatEquation&>(myProblem);
  ElmNorm& elmNorm = static_cast<ElmNorm&>(elmInt);

  const Material* mat = problem.getMaterial();

  double theta = fe.N.dot(elmNorm.vec[0]);
  double theta0 = problem.initialTemperature(X);
  double rhocp = mat?mat->getMassDensity(X)*mat->getHeatCapacity(theta):1.0;

  elmNorm[0] = rhocp*(theta-theta0)*fe.detJxW;

  return true;
}
