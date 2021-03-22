#!/bin/sh

# Test the program by mocking the qchem binary
EXTGAUSS_QCHEM_RUNDIR=$(realpath test-data)
EXTGAUSS_QCHEM_EXE=$(realpath test-data/qchem)
./qchem_g16 --rem $(realpath test-data/params.rem) 1 test-data/Gaussian.EIn Gaussian.EOu Gaussian.Em null null
