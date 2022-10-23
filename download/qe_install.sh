#!/bin/bash

echo "Create directory"
mkdir -p QE
cd QE/
echo "COMPLETE ..."
echo " "

echo "Download QE file "
curl -O https://gitlab.com/QEF/q-e/-/archive/qe-7.0/q-e-qe-7.0.tar.gz
echo "COMPLETE ..."
echo " "

echo "Unzip QE file"
tar zxf q-e-qe-7.0.tar.gz
mv q-e-qe-7.0 qe-7.0
echo "COMPLETE ..."
echo " "

# Navigate to build directory
mkdir -p ./qe-7.0/_build
cd ./qe-7.0/_build

# Start build
# This is a serial build because I found errors when doing parallel
echo "Build QE programs"

>>qe_build.log 2>&1

echo "#########################################################################" >> qe_build.log 2>&1
echo "                                BUILD                                    " >> qe_build.log 2>&1
echo "#########################################################################" >> qe_build.log 2>&1
cmake \
-DQE_ENABLE_MPI=ON \
-DQE_ENABLE_TEST=ON \
-DQE_ENABLE_SCALAPACK=ON \
-DQE_FFTW_VENDOR=Intel_DFTI \
-DCMAKE_C_COMPILER=mpiicc \
-DCMAKE_Fortran_COMPILER=mpiifort \
-DCMAKE_INSTALL_PREFIX=../_install \
../ >> qe_build.log 2>&1

echo "#########################################################################" >> qe_build.log
echo "COMPLETE ..."
echo " "

# INSTALL 
echo "Installing QE programs (Grab some coffee, this will take a moment)"
echo "#########################################################################" >> qe_build.log
echo "                                INSTALL                                  " >> qe_build.log
echo "#########################################################################" >> qe_build.log
make install >> qe_build.log 2>&1
echo "#########################################################################" >> qe_build.log
echo "COMPLETE ..."
echo " "

echo " "
echo "!! INSTALLATION COMPLETE !!"
echo " "

# Print Usage guide
echo "- - - - - - - - - - - - - - - - - - - - - -"
echo "USAGE GUIDE"
echo " "
cd ../_install/bin
TMPBIN=$(pwd)
echo "add to ~/.bashrc:"
echo " "
echo "   export PATH="$TMPBIN':$PATH'
echo " "
echo "- - - - - - - - - - - - - - - - - - - - - -"
echo " "








