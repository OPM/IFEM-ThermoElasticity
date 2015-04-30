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
#include "LinIsotropic.h"

/*!
  \brief Class representing the integrand of the heat equation.

  \details Time stepping is done using BDF1/BDF2
*/

class HeatEquation : public IntegrandBase
{
public:
  typedef LinIsotropic MaterialType; //!< Material used in this integrand

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

#endif
