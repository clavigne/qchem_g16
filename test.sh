#!/bin/sh

# Test the program by mocking the qchem binary
mkdir -p fake-qchem-dir
cp test-data/* fake-qchem-dir

export EXTGAUSS_QCHEM_RUNDIR=$(realpath fake-qchem-dir)
export EXTGAUSS_QCHEM_EXE=$(realpath test-data/qchem)
./qchem_g16 --rem $(realpath test-data/params.rem) 1 test-data/Gaussian.EIn Gaussian.EOu Gaussian.Em null null
