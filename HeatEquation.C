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
#include "HeatQuantities.h"
#include "WeakOperators.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "MaterialBase.h"
#include "Vec3Oper.h"
#include "AnaSol.h"


HeatEquation::HeatEquation (unsigned short int n, int order)
  : nsd(n), bdf(order), mat(NULL), flux(NULL), init(NULL)
{
  primsol.resize(order+1);
  sourceTerm = nullptr;
}


bool HeatEquation::evalInt (LocalIntegral& elmInt,
                            const FiniteElement& fe,
                            const TimeDomain& time,
                            const Vec3& X) const
{
  Matrix& A = static_cast<ElmMats&>(elmInt).A.front();
  Vector& b = static_cast<ElmMats&>(elmInt).b.front();

  double theta = 0.0;
  double rhocp = 1.0, kappa = 1.0;
  for (int t = 1; t <= bdf.getOrder(); t++) {
    double val = fe.N.dot(elmInt.vec[t]);
    if (t == 1 && mat) {
      rhocp = mat->getMassDensity(X)*mat->getHeatCapacity(val);
      kappa = mat->getThermalConductivity(val);
    }
    theta -= bdf[t]/time.dt*val;
  }

  WeakOperators::Laplacian(A,fe,kappa);
  WeakOperators::Mass(A,fe,rhocp*bdf[0]/time.dt);
  WeakOperators::Source(b,fe,rhocp*theta+this->getSource(X));

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

  // Evaluate the Neumann value
  double T = (*flux)(X);
  double val = fe.N.dot(elmInt.vec.front());
  double kappa = mat ? mat->getThermalConductivity(val) : 1.0;

  // Integrate the Neumann value
  WeakOperators::Source(static_cast<ElmMats&>(elmInt).b.front(),fe,kappa*T);

  return true;
}


std::string HeatEquation::getField1Name (size_t, const char* prefix) const
{
  if (!prefix)
    return "theta";

  return prefix + std::string(" theta");
}


std::string HeatEquation::getField2Name (size_t i, const char* prefix) const
{
  if (i > 2)
    return "";

  static const char* s[3] = { "theta_x", "theta_y", "theta_z" };
  if (!prefix)
    return s[i];

  return prefix + std::string(" ") + s[i];
}


LocalIntegral* HeatEquation::WeakDirichlet::getLocalIntegral (size_t nen,
                                                              size_t,
                                                              bool) const
{
  ElmMats* result = new ElmMats(true);
  result->resize(1,1);
  result->redim(nen);

  return result;
}


bool HeatEquation::WeakDirichlet::evalBou (LocalIntegral& elmInt,
                                           const FiniteElement& fe,
                                           const Vec3& X,
                                           const Vec3& normal) const
{
  if (!flux) {
    std::cerr <<" *** HeatEquation::evalBou: No flux function."<< std::endl;
    return false;
  }

  Matrix& A = static_cast<ElmMats&>(elmInt).A.front();
  Vector& b = static_cast<ElmMats&>(elmInt).b.front();

  // Evaluate the Neumann value
  double q = (*flux)(X);
  double val = fe.N.dot(elmInt.vec.front());
  double kappa = mat ? mat->getThermalConductivity(val) : 1.0;

  WeakOperators::Mass(A,fe,-envCond);
  WeakOperators::Source(b,fe,q);

  for (size_t i = 1; i <= fe.N.size(); i++) {
    for (size_t j = 1; j <= fe.N.size(); j++) {
      double dT = normal * fe.dNdX.getRow(j);
      A(i,j) += (kappa*dT+envT*envCond)*fe.N(i)*fe.detJxW;
    }
    b(i) += q*fe.N(i)*fe.detJxW;
  }

  return true;
}


ForceBase* HeatEquation::getForceIntegrand (const Vec3*, AnaSol*) const
{
  return new HeatEquationFlux<HeatEquation>(*const_cast<HeatEquation*>(this));
}


NormBase* HeatEquation::getNormIntegrand (AnaSol* asol) const
{
  return new HeatEquationNorm(*const_cast<HeatEquation*>(this),asol);
}


HeatEquationNorm::HeatEquationNorm (HeatEquation& p, AnaSol* a) : NormBase(p)
{
  nrcmp = myProblem.getNoFields(2);
  anasol = a;
}


bool HeatEquationNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X) const
{
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  HeatEquation& hep = static_cast<HeatEquation&>(myProblem);
  const Material* mat = hep.getMaterial();

  // Evaluate the FE temperature and thermal conductivity at current point
  double Uh = fe.N.dot(elmInt.vec.front());
  double kappa = mat ? mat->getThermalConductivity(Uh) : 1.0;

  // Evaluate the FE heat flux vector, gradU = dNdX^T * eV
  Vector gradUh;
  if (!fe.dNdX.multiply(elmInt.vec.front(),gradUh,true))
    return false;

  size_t ip = 0;
  // Integrate the L2 norm, (U^h, U^h)
  pnorm[ip++] += Uh*Uh*fe.detJxW;

  // Integrate the energy norm, a(U^h,U^h)
  pnorm[ip++] = 0.5*kappa*gradUh.dot(gradUh)*fe.detJxW;
  ip++; // Currently no external energy yet

  if (anasol && anasol->getScalarSol()) {
    double T = (*anasol->getScalarSol())(X);
    pnorm[ip++] += T*T*fe.detJxW; // L2 norm of analytical solution
    pnorm[ip++] += (T-Uh)*(T-Uh)*fe.detJxW; // L2 norm of error
  }

  if (anasol && anasol->getScalarSecSol()) {
    Vec3 dT = (*anasol->getScalarSecSol())(X);
    pnorm[ip++] += 0.5*kappa*dT*dT*fe.detJxW;
    pnorm[ip++] += 0.5*kappa*(dT-gradUh)*(dT-gradUh)*fe.detJxW;
  }


  // TODO: Add energy-norm of exact error based on analytical solution
  size_t i, j;
  for (i = 0; i < pnorm.psol.size(); i++)
    if (!pnorm.psol[i].empty())
    {
      // Evaluate projected heat flux field
      Vector gradUr(nrcmp);
      for (j = 0; j < nrcmp; j++)
        gradUr[j] = pnorm.psol[i].dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(U^r,U^r)
      pnorm[ip++] += 0.5*kappa*gradUr.dot(gradUr)*fe.detJxW;
      // Integrate the estimated error in energy norm a(U^r-U^h,U^r-U^h)
      Vector error = gradUr - gradUh;
      pnorm[ip++] += 0.5*kappa*error.dot(error)*fe.detJxW;
    }

  return true;
}


size_t HeatEquationNorm::getNoFields (int group) const
{
  size_t nf = 1;
  if (group < 1)
    for (size_t i = 0; i < prjsol.size(); i++)
      nf += prjsol.empty() ? 0 : 1;
  else
    nf = anasol ? 7 : 3;

  return nf;
}


std::string HeatEquationNorm::getName (size_t i, size_t j,
                                       const char* prefix) const
{
  if (i == 0 || j == 0 || j > 4)
    return this->NormBase::getName(i,j,prefix);

  static const char* s[] = {
    "(theta^h,theta^h)^0.5",
    "a(theta^h,theta^h)^0.5",
    "(q,theta^h)^0.5",
    "(theta, theta)^0.5",
    "a(theta,theta)^0.5",
    "a(e,e)^0.5, e=theta-theta^h",
    "a(theta^r,theta^r)^0.5",
    "a(e,e)^0.5, e=theta^r-theta^h",
    "a(e,e)^0.5, e=theta-theta^r",
    "effectivity index"
  };

  size_t k = i > 1 ? j+3 : j-1;

  if (!prefix)
    return s[k];

  return prefix + std::string(" ") + s[k];
}


bool HeatEquationNorm::hasElementContributions (size_t i, size_t j) const
{
  return i > 1 || j != 2;
}


double HeatEquation::getSource (const Vec3& X) const
{
  if (sourceTerm)
    return (*sourceTerm)(X);
  return 0.0;
}
