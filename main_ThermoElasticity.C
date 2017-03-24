// $Id$
//==============================================================================
//!
//! \file main_ThermoElasticity.C
//!
//! \date 05 August 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Main program for an isogeometric thermo-elastic solver.
//!
//==============================================================================

#include "IFEM.h"
#include "AppCommon.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMSolver.h"
#include "SIMHeatEquation.h"
#include "SIMThermalCoupling.h"
#include "SIMThermoElasticity.h"
#include "HDF5Writer.h"
#include "HeatEquation.h"
#include "XMLWriter.h"
#include "TimeIntUtils.h"
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


//! \brief Setup and launch the simulation.
//! \param[in] infile The input file to process
//! \param[in] tit The time integration method to use. Either BE or BDF2
template<class Dim>
int runSimulator(char* infile, TimeIntegration::Method tIt)
{
  typedef SIMHeatEquation<Dim,HeatEquation>               HeatSolver;
  typedef SIMThermoElasticity<Dim>                        ElasticitySolver;
  typedef SIMThermalCoupling<HeatSolver,ElasticitySolver> CoupledSolver;

  HeatSolver               tempModel(TimeIntegration::Order(tIt));
  ElasticitySolver         solidModel;
  CoupledSolver            model(tempModel,solidModel);
  SIMSolver<CoupledSolver> solver(model);

  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n=========================\n";

  if (ConfigureSIM(tempModel, infile) ||
      ConfigureSIM(solidModel, infile) || !solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  const SIMoptions& sopt = tempModel.opt;

  model.setupDependencies();
  if (sopt.restartFile.empty())
    model.init(solver.getTimePrm());
  else if (solver.restart(sopt.restartFile,sopt.restartStep) < 0)
    return 2;

  DataExporter* exporter = nullptr;
  if (tempModel.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(model,solver,sopt.hdf5,false,
                                     sopt.saveInc,sopt.restartInc);

  int res = solver.solveProblem(infile,exporter);

  delete exporter;
  return res;
}


/*!
  \brief Main program for an isogeometric thermo-elastic solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -2D : Use two-parametric simulation driver (plane stress)
  \arg -2Dpstrain : Use two-parametric simulation driver (plane strain)
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  bool twoD = false;
  char* infile = nullptr;
  TimeIntegration::Method tIt = TimeIntegration::BDF2;
  Elasticity::wantPrincipalStress = true;

  IFEM::Init(argc,argv,"Thermo-Elasticity solver");

  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strncmp(argv[i],"-2Dpstra",8))
      twoD = SIMElasticity<SIM2D>::planeStrain = true;
    else if (!strncmp(argv[i],"-2D",3))
      twoD = true;
    else if (!strncmp(argv[i],"-msg",4) && i < argc-1)
      SIMadmin::msgLevel = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-be"))
      tIt = TimeIntegration::BE;
    else if (!strcmp(argv[i],"-bdf2"))
      tIt = TimeIntegration::BDF2;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D[pstrain]] [-nGauss <n>]\n"
	      <<"       [-hdf5] [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]]\n";
    return 0;
  }

  std::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout) << std::endl;
  utl::profiler->stop("Initialization");

  if (twoD)
    return runSimulator<SIM2D>(infile,tIt);
  else
    return runSimulator<SIM3D>(infile,tIt);
}
