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

#ifndef _HEAT_EQUATION_H
#define _HEAT_EQUATION_H

#include "IntegrandBase.h"
#include "LinIsotropic.h"
#include "BDF.h"


/*!
  \brief Class representing the integrand of the heat equation.
  \details Time stepping is done using BDF1/BDF2
*/

class HeatEquation : public IntegrandBase
{
public:
  typedef LinIsotropic MaterialType; //!< Material used in this integrand

  //! \brief Class representing the weak Dirichlet integrand.
  class WeakDirichlet : public IntegrandBase
  {
  public:
    //! \brief Default constructor.
    //! \param[in] n Number of spatial dimensions
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
    void setMaterial(Material* material) { mat = material; }
    //! \brief Defines the flux function.
    void setFlux(RealFunc* f) { flux = f; }
    //! \brief Sets temperature of environment.
    void setEnvTemperature(double T) { envT = T; }
    //! \brief Sets conductivity of environment.
    void setEnvConductivity(double alpha) { envCond = alpha; }

  private:
    size_t nsd;     //!< Number of spatial dimensions
    RealFunc* flux; //!< Flux function
    Material* mat;  //!< Material parameters
    double envT;    //!< Temperature of environment
    double envCond; //!< Conductivity of environment
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
  void advanceStep() { bdf.advanceStep(); }

  //! \brief Defines the material properties.
  void setMaterial(Material* material) { mat = material; }

  //! \brief Defines the source term
  void setSource(RealFunc* src) { sourceTerm = src; }

  //! \brief Evaluates the source term (if any) at a specified point.
  //! \param[in] X Cartesian coordinate of current integration point
  double getSource(const Vec3& X) const;

  //! \brief Obtain the current material.
  const Material* getMaterial() const { return mat; }

  //! \brief Defines the flux function.
  void setFlux(RealFunc* f) { flux = f; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual std::string getField1Name(size_t, const char* prefix = NULL) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix = NULL) const;

  //! \brief Sets the function of the initial temperature field.
  void setInitialTemperature(const RealFunc* f)  { init = f; }
  //! \brief Returns the function of the initial temperature field.
  const RealFunc* getInitialTemperature() const { return init; }

  //! \brief Returns number of spatial dimensions.
  size_t getNoSpaceDim() const { return nsd; }

  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  virtual ForceBase* getForceIntegrand(const Vec3*, AnaSol* = NULL) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] asol Pointer to analytical solution field (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = NULL) const;

  //! \brief Returns the initial temperature in a point.
  //! \param[in] X The coordinate of the point.
  //! \returns Initial temperature
  double initialTemperature(const Vec3& X) const { return init?(*init)(X):0.0; }

private:
  size_t nsd;               //!< Number of spatial dimensions
  TimeIntegration::BDF bdf; //!< BDF helper class
  Material* mat;            //!< Material parameters
  RealFunc* flux;           //!< Pointer to the heat flux field
  const RealFunc* init;     //!< Initial temperature function
  RealFunc* sourceTerm;     //!< Pointer to source term
};


/*!
  \brief Class representing the integrand of heat equation energy norms.
*/

class HeatEquationNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The heat equation problem to evaluate norms for
  //! \param[in] a The analytical aolution (optional)
  HeatEquationNorm(HeatEquation& p, AnaSol* a = NULL);
  //! \brief Empty destructor.
  virtual ~HeatEquationNorm() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix) const;

  //! \brief Returns whether a norm quantity stores element contributions.
  virtual bool hasElementContributions(size_t i, size_t j) const;

private:
  AnaSol* anasol; //!< Analytical solution
};

#endif
