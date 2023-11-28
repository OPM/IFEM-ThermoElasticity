// $Id$
//==============================================================================
//!
//! \file SIMHeatEquation.h
//!
//! \date Aug 19 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for the Heat equation.
//!
//==============================================================================

#ifndef _SIM_HEAT_EQUATION_H_
#define _SIM_HEAT_EQUATION_H_

#include "SIMsolution.h"
#include "SIMconfigure.h"


class DataExporter;
class Material;
class RealFunc;
class SIMoutput;
class TimeStep;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Driver class for a heat equation simulator.
  \details The class incapsulates data and methods for solving a
  heat equation problem using NURBS-based finite elements.
*/

template<class Dim, class Integrand>
class SIMHeatEquation : public Dim, public SIMsolution
{
  //! \brief Struct containing parameters for boundary heat flux calculation.
  struct BoundaryFlux
  {
    std::string file; //!< Name of output file for boundary flux calculation
    std::string set;  //!< Name of topology set for boundary flux calculation
    int code;         //!< Code identifying the boundary for flux calculation
    int timeIncr;     //!< Time level increment for boundary flux calculation
    //! \brief Default constructor.
    BoundaryFlux(int c = 0) : code(c), timeIncr(1) {}
    //! \brief Constructor naming the topology set for force calculation.
    explicit BoundaryFlux(const std::string& s) : set(s), code(0), timeIncr(1) {}
  };

  //! \brief Helper class for searching among BoundaryForce objects.
  class hasCode
  {
    int myCode; //!< The property code to compare with
  public:
    //! \brief Constructor initializing the property code to search for.
    explicit hasCode(int code) : myCode(abs(code)) {}
    //! \brief Returns \e true if the BoundaryForce \a b has the code \a myCode.
    bool operator()(const BoundaryFlux& b) { return abs(b.code) == myCode; }
  };

public:
  //! \brief Struct containing setup parameters.
  struct SetupProps
  {
    bool shareGrid = false;     //!< \e true to share grid with another simulator
    SIMoutput* share = nullptr; //!< Simulator to share grid with
  };

  //! \brief Default constructor.
  //! \param[in] order Order of temporal integration (1 or 2)
  explicit SIMHeatEquation(int order);

  //! \brief The destructor zero out the integrand pointer (deleted by parent).
  virtual ~SIMHeatEquation();

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "HeatEquation"; }

  //! \brief Initializes the temperature solution vectors.
  void initSol();

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock);

  //! \brief Dummy method.
  bool init(const TimeStep&) { return true; }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&);

  //! \brief Computes the solution for the current time step.
  bool solveStep(const TimeStep& tp);

  //! \brief Evaluates and prints out solution norms.
  void printFinalNorms(const TimeStep& tp);

  //! \brief Compute and save boundary heat flux or stored energy in a volume
  //! \param[in] bf Description of integration domain
  //! \param[in] tp Time stepping information
  //! \param[in] flux True to calculate a heat flux, false for stored energy
  bool saveIntegral(const BoundaryFlux& bf, const TimeStep& tp, bool flux);

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock);

  using SIMsolution::getSolution;
  //! \brief Returns a reference to current solution vector.
  Vector& getSolution() { return solution.front(); }

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter, const std::string& prefix="");

  //! \brief Set context to read from input file
  void setContext(int ctx);

  //! \brief Sets the function of the initial temperature field.
  void setInitialTemperature(const RealFunc* f);
  //! \brief Returns the function of the initial temperature field.
  const RealFunc* getInitialTemperature() const;

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(SerializeMap& data) const override;

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const SerializeMap& data) override;

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the Robin property codes.
  void preprocessA() override;

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to ensure that threading groups are
  //! established for the patch faces subjected to boundary flux integration.
  bool preprocessB() override;

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  bool initMaterial(size_t propInd) override;

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override;

  Integrand                   heq;  //!< Main integrand

private:
  typename Integrand::WeakDirichlet wdc; //!< Weak dirichlet integrand
  std::vector<Material*>      mVec; //!< Material data
  RealFunc*                   srcF; //!< Source function

  std::string inputContext; //!< Input context

  std::vector<BoundaryFlux> fluxes;  //!< Heat fluxes to calculate
  std::vector<BoundaryFlux> senergy; //!< Stored energies to calculate
  int aCode[2] = {0}; //!< Analytical BC code (used by destructor)
};


//! \brief Partial specialization for configurator
template<class Dim, class Integrand>
struct SolverConfigurator< SIMHeatEquation<Dim,Integrand> > {
  //! \brief Setup a heat equation simulator.
  //! \param sim The simulator to set up
  //! \param[in] props Setup properties
  //! \param[in] infile The input file to parse
  int setup(SIMHeatEquation<Dim,Integrand>& sim,
            const typename SIMHeatEquation<Dim,Integrand>::SetupProps& props, char* infile);
};

#endif
