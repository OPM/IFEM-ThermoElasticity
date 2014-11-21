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

  size_t getNoSpaceDim() const { return nsd; }

  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] X0 Reference point for torque computation
  //! \param[in] asol Pointer to analytical solution field (optional)
  virtual ForceBase* getForceIntegrand(const Vec3* X0, AnaSol* asol = 0) const;

protected:
  size_t nsd;                //!< Number of spatial dimensions
  TimeIntegration::BDF bdf;  //!< BDF helper class
  Material* mat;             //!< Material parameters
  RealFunc* flux;            //!< Pointer to the flux field
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

#endif
