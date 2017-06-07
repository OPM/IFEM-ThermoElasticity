// $Id$
//==============================================================================
//!
//! \file HeatQuantities.h
//!
//! \date April 28 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for heat-derived quantities.
//!
//==============================================================================

#ifndef HEAT_QUANTITIES_H
#define HEAT_QUANTITIES_H

#include "ElmNorm.h"
#include "FiniteElement.h"
#include "IntegrandBase.h"
#include "Vec3Oper.h"


/*!
  \brief Class representing the integrand for computing boundary heat fluxes.
*/

template<class HE> class HeatEquationFlux : public ForceBase
{
public:
  //! \brief Constructor for global force resultant integration.
  //! \param[in] p The heat equation problem to evaluate fluxes for
  explicit HeatEquationFlux(HE& p) : ForceBase(p) {}

  //! \brief Empty destructor.
  virtual ~HeatEquationFlux() {}

  using ForceBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time,
                       const Vec3& X, const Vec3& normal) const
  {
    HE& problem = static_cast<HE&>(myProblem);
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

  //! \brief Returns the number of force components.
  virtual size_t getNoComps() const { return 1; }
};


/*!
  \brief Class representing the integrand for computing stored energy in a topology set.
*/

template<class HE> class HeatEquationStoredEnergy : public ForceBase
{
public:
  //! \brief Constructor for global force resultant integration.
  //! \param[in] p The heat equation problem to evaluate fluxes for
  explicit HeatEquationStoredEnergy(HE& p) : ForceBase(p) {}

  //! \brief Empty destructor.
  virtual ~HeatEquationStoredEnergy() {}

  using ForceBase::evalInt;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X) const
  {
    HE& problem = static_cast<HE&>(myProblem);
    ElmNorm& elmNorm = static_cast<ElmNorm&>(elmInt);

    const Material* mat = problem.getMaterial();

    double theta = fe.N.dot(elmNorm.vec[0]);
    double theta0 = problem.initialTemperature(X);
    double rhocp = mat?mat->getMassDensity(X)*mat->getHeatCapacity(theta):1.0;

    elmNorm[0] = rhocp*(theta-theta0)*fe.detJxW;

    return true;
  }

  //! \brief Returns the number of force components.
  virtual size_t getNoComps() const { return 1; }

  using ForceBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen1 Number of nodes on element for basis 1
  //! \param[in] nen2 Number of nodes on element for basis 2
  //! \param[in] iEl Global element number (1-based)
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  //!
  //! \details This form is used for mixed formulations only.
  //! The default implementation just forwards to the single-basis version.
  //! Reimplement this method if your mixed formulation requires specialized
  //! local integral objects.
  virtual LocalIntegral* getLocalIntegral(size_t nen1, size_t nen2, size_t iEl,
                                          bool neumann = false) const
  {
    return this->getLocalIntegral(nen1,iEl,neumann);
  }

  using ForceBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  //! Reimplement this method for problems requiring the element center and/or
  //! the number of integration points during/before the integrand evaluations.
  virtual bool initElement(const std::vector<int>& MNPC,
			   const FiniteElement& fe,
			   const Vec3& X0, size_t nPt,
			   LocalIntegral& elmInt)
  {
    return myProblem.initElement(MNPC,elmInt);
   }

  //! \brief This is a volume integrand.
  virtual bool hasInteriorTerms() const { return true; }
};

#endif
