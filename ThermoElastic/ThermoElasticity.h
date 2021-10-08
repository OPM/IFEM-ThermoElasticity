// $Id$
//==============================================================================
//!
//! \file ThermoElasticity.h
//!
//! \date Aug 07 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for linear thermo-elasticity problems.
//!
//==============================================================================

#ifndef _THERMO_ELASTICITY_H
#define _THERMO_ELASTICITY_H

#include "LinearElasticity.h"


/*!
  \brief Class representing the integrand of a linear thermo-elasticity problem.
  \details Most methods of this class are inherited form the parent class.
  Only the \a getThermalStrain, \a formInitStrainForces and \a evalSol methods
  are reimplemented, to account for strains due to a varying temperature field.
*/

class ThermoElasticity : public LinearElasticity
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, an axisymmetric 3D formulation is assumed
  ThermoElasticity(unsigned short int n, bool axS = false);
  //! \brief Empty destructor.
  virtual ~ThermoElasticity() {}

  using LinearElasticity::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC, const FiniteElement&,
                   const Vec3&, size_t, LocalIntegral& elmInt) override;

  using LinearElasticity::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  bool evalSol(Vector& s, const FiniteElement& fe,
               const Vec3& X, const std::vector<int>& MNPC) const override;

protected:
  //! \brief Evaluates the thermal strain at current integration point.
  //! \param[in] eT Element temperature vector
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current integration point
  double getThermalStrain(const Vector& eT, const Vector& N,
                          const Vec3& X) const override;

  //! \brief Calculates integration point initial strain force contributions.
  //! \param elMat Element matrices for current element
  //! \param[in] N Basis function values at current point
  //! \param[in] B Strain-displacement matrix
  //! \param[in] C Constitutive matrix
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  bool formInitStrainForces(ElmMats& elMat, const Vector& N,
                            const Matrix& B, const Matrix& C,
                            const Vec3& X, double detJW) const override;

private:
  Vector myTempVec; //!< Current temperature at nodal points
};

#endif
