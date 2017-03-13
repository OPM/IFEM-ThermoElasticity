# IFEM ThermoElasticity

## Introduction

This module contains ThermoElasticity libraries and applications built
using the IFEM library.

### Getting all dependencies

1. Install IFEM from https://github.com/OPM/IFEM

### Getting the code

This is done by first navigating to the folder in which you want the application and typing

    git clone https://github.com/OPM/IFEM-Elasticity Elasticity
    git clone https://github.com/OPM/IFEM-ThermoElasticity ThermoElasticity

### Compiling the code

To compile, first navigate to the root catalogue `ThermoElasticity`.

1. `cd ThermoElasticity`
2. `mkdir Debug`
3. `cd Debug`
5. `cmake -DCMAKE_BUILD_TYPE=Debug ..`
6. `make `

this will compile the library and the thermoelasticity application. 
The binaries can be found in the 'bin' subfolder.
Change all instances of `Debug` with `Release` to drop debug-symbols, 
but get faster running code.

### Testing the code

IFEM is using cmake test system. To compile run all regression- and unit-tests, navigate to your build 
folder (i.e. `ThermoElasticity/Debug`) and type

    make check
