// $Id$
//==============================================================================
//!
//! \file HeatEquation.h
//!
//! \date Aug 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for the heat equation.
//!
//==============================================================================

#ifndef _ADVECTION_DIFFUSION_BDF_H
#define _ADVECTION_DIFFUSION_BDF_H

#include "IntegrandBase.h"
#include "BDF.h"

class Material;

/*!
  \brief Class representing the integrand of the heat equation.

  \details Time stepping is done using BDF1/BDF2
*/

class HeatEquation : public IntegrandBase
{
public:
  class WeakDirichlet : public IntegrandBase
  {
    public:
      //! \brief Default constructor.
      //! \param[in] n Number of spatial dimensions
      //! \param[in] form The solution formulation to use
      WeakDirichlet(unsigned short int n) :
        nsd(n), flux(nullptr), mat(nullptr), envT(273.5), envCond(1.0) {}

      //! \brief Empty destructor.
      virtual ~WeakDirichlet() {}

      //! \brief Returns that this integrand has no interior contributions.
      virtual bool hasInteriorTerms() const { return false; }
      //! \brief Returns a local integral contribution object for given element.
      //! \param[in] nen Number of nodes on element
      virtual LocalIntegral* getLocalIntegral(size_t nen, size_t, bool) const;
      //! \brief Evaluates the integrand at a boundary point.
      //! \param elmInt The local integral object to receive the contributions
      //! \param[in] fe Finite element data of current integration point
      //! \param[in] X Cartesian coordinates of current integration point
      //! \param[in] normal Boundary normal vector at current integration point
      virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                           const Vec3& X, const Vec3& normal) const;

      //! \brief Defines the material properties.
      virtual void setMaterial(Material* material) { mat = material; }
      //! \brief Defines the flux function.
      void setFlux(RealFunc* f) { flux = f; }
      //!< \brief Set temperature of environment
      void setEnvTemperature(double T) { envT = T; }
      //!< \brief Set conductivity of environment
      void setEnvConductivity(double alpha) { envCond = alpha; }
    protected:
      size_t nsd;       //!< Number of spatial dimensions
      RealFunc* flux;   //!< Flux function
      Material* mat;    //!< Material parameters
      double envT;      //!< Temperature of environment
      double envCond;   //!< Conductivity of environment
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] order Temporal order (1,2)
  HeatEquation(unsigned short int n = 3, int order = 1);

  //! \brief Empty destructor.
  virtual ~HeatEquation() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const TimeDomain& time, const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Advance time stepping scheme.
  virtual void advanceStep() { bdf.advanceStep(); }

  //! \brief Defines the material properties.
  virtual void setMaterial(Material* material) { mat = material; }

  //! \brief Obtain the current material.
  const Material* getMaterial() const { return mat; }

  //! \brief Defines the flux function.
  void setFlux(RealFunc* f) { flux = f; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual const char* getField1Name(size_t, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

  //! \brief Set the function for the initial temperature field
  //! \param f The function
  void setInitialTemperature(const RealFunc* f)  { init = f; }

  //! \brief Returns number of spatial dimensions.
  size_t getNoSpaceDim() const { return nsd; }

  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] X0 Reference point for torque computation
  //! \param[in] asol Pointer to analytical solution field (optional)
  virtual ForceBase* getForceIntegrand(const Vec3* X0, AnaSol* asol = 0) const;

  //! \brief Returns the initial temperature in a point.
  //! \param[in] X The coordinate of the point.
  //! \returns Initial temperature
  double initialTemperature(const Vec3& X) const { return init?(*init)(X):0.0; }

protected:
  size_t nsd;                //!< Number of spatial dimensions
  TimeIntegration::BDF bdf;  //!< BDF helper class
  Material* mat;             //!< Material parameters
  RealFunc* flux;            //!< Pointer to the flux field
  const RealFunc* init;      //!< Initial temperature function
};

/*!
  \brief Class representing the integrand for computing boundary heat fluxes.
*/

class HeatEquationFlux : public ForceBase
{
public:
  //! \brief Constructor for global force resultant integration.
  //! \param[in] p The heat equation problem to evaluate fluxes for
  HeatEquationFlux(HeatEquation& p)
    : ForceBase(p) {}

  //! \brief Empty destructor.
  virtual ~HeatEquationFlux() {}

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns the number of force components.
  virtual size_t getNoComps() const { return 1; }
};

class HeatEquationStoredEnergy : public ForceBase
{
public:
  //! \brief Constructor for global force resultant integration.
  //! \param[in] p The heat equation problem to evaluate fluxes for
  HeatEquationStoredEnergy(HeatEquation& p)
    : ForceBase(p) {}

  //! \brief Empty destructor.
  virtual ~HeatEquationStoredEnergy() {}

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X) const;

  //! \brief Returns the number of force components.
  virtual size_t getNoComps() const { return 1; }

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
  { return myProblem.initElement(MNPC, elmInt); }

  //! \brief This is a volume integrand
  virtual bool hasInteriorTerms() const { return true; }
};

#endif
