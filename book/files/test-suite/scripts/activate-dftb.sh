#!/bin/bash

installdir=/home/krojas/tutorial_files/apps/dftbplus-22.1-OpenMP/_install

export PATH=${installdir}/bin:$PATH
export LD_LIBRARY_PATH=${installdir}/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=${installdir}/lib/python3.8/site-packages/pythonapi-0.1-py3.8.egg:$PYTHONPATH
export DFTB_LIB=${installdir}/lib64