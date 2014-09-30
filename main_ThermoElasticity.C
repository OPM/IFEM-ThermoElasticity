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
#include "SIMThermoElasticity.h"
#include "SIMElasticityWrap.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "TimeIntUtils.h"


//! \brief Setup and launch simulation
//! \param[in] infile The input file to process
//! \param[in] restartfile File to restart from. NULL for no restart
//! \param[in] tit The time integration method to use. Either BE or BDF2
  template<class Dim>
int runSimulator(char* infile, char* restartfile, TimeIntegration::Method tIt)
{
  SIMHeatEquation<Dim> tempModel(tIt==TimeIntegration::BE?1:2);
  SIMElasticityWrap<Dim> solidModel;
  SIMThermoElasticity< SIMHeatEquation<Dim>, SIMElasticityWrap<Dim> > model(tempModel, solidModel);
  SIMSolver< SIMThermoElasticity< SIMHeatEquation<Dim>, SIMElasticityWrap<Dim> > > solver(model);

  if (ConfigureSIM(tempModel, infile)  ||
      ConfigureSIM(solidModel, infile) || !solver.read(infile))
    return 1;

  if (restartfile)
    SIM::handleRestart(model, solver, restartfile, tempModel.getDumpInterval(),
                       TimeIntegration::Steps(tIt));

  model.setupDependencies();
  model.init(solver.getTimePrm());

  DataExporter* exporter=NULL;
  if (tempModel.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(model, solver, tempModel.opt.hdf5,
                                     restartfile && tempModel.opt.hdf5 == restartfile,
                                     tempModel.getDumpInterval(),
                                     TimeIntegration::Steps(tIt));

  if (!solver.solveProblem(infile, exporter))
    return 5;

  delete exporter;
  return 0;
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
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SIMoptions dummy;
  int i;
  bool twoD = false;
  char* infile = NULL;
  char* restartfile = NULL;
  TimeIntegration::Method tIt = TimeIntegration::BDF2;

  int myPid = IFEM::Init(argc, argv);

  for (i = 1; i < argc; i++)
    if (dummy.parseOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      twoD = true;
    else if (!strncmp(argv[i],"-msg",4) && i < argc-1)
      SIMinput::msgLevel = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-be"))
      tIt = TimeIntegration::BE;
    else if (!strcmp(argv[i],"-bdf2"))
      tIt = TimeIntegration::BDF2;
    else if (!strcmp(argv[i],"-restart") && i < argc-1)
      restartfile = strtok(argv[++i],".");
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>]\n"
	      <<"       [-hdf5] [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]]\n";
    return 0;
  }

  if (myPid == 0)
  {
    const SIMoptions& opts = IFEM::getOptions();
    std::cout <<"\n >>> IFEM Thermo-Elasticity solver <<<"
	      <<"\n ====================================\n"
	      <<"\n Executing command:\n";
    for (i = 0; i < argc; i++) std::cout <<" "<< argv[i];
    std::cout <<"\n\nInput file: "<< infile
	      <<"\nEquation solver: "<< opts.solver
	      <<"\nNumber of Gauss points: "<< opts.nGauss[0];
    if (opts.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (opts.format ? "BINARY":"ASCII")
                <<"\nNumber of visualization points: "<< opts.nViz[0]
                <<" "<< opts.nViz[1];
      if (!twoD) std::cout <<" "<< opts.nViz[2];
    }

    if (opts.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (opts.discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    else if (opts.discretization == ASM::LRSpline)
      std::cout <<"\nLR-spline basis functions are used";
    std::cout << std::endl;
  }
  utl::profiler->stop("Initialization");

  if (twoD)
    return runSimulator<SIM2D>(infile, restartfile, tIt);
  else
    return runSimulator<SIM3D>(infile, restartfile, tIt);

  return 1; // Should not be here
}
