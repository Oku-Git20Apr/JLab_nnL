#! /bin/bash

cp nnl_HRSAnalysis.cc ../src/HRSAnalysis.cc
cp nnl_HRSDetectorConstruction.cc ../src/HRSDetectorConstruction.cc
cp nnl_HRSPrimaryGeneratorAction.cc ../src/HRSPrimaryGeneratorAction.cc
cd ../build/
make -j
cd -
