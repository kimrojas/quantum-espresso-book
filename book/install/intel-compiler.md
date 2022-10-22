# Intel Compiler and MKL

## A. In Supercomputer

In the usual working environment, you will be using a supercomputer to run the simulations. Here, the administator usually provides a ready to use Intel compiler and MKL library. In this case, you can just use the one provided.

Example: In my supercomputer, I can activate the provided Intel package using

```bash
module load intel intelmkl
```

*Please consult your administrator about this.*


## B. In computer without Intel Package

In the case that you are installing in a local computer or a supercomputer without Intel Compiler & MKL, please follow the instructions here.

There are two main files that you need to install. 
1. Intel oneAPI Base toolkit ([link](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html))
   - Python (with conda)
   - C compiler
   - Math Kernel Library (MKL)
2. Intel oneAPI HPC toolkit  ([link](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html))
   - Fortran compiler
   - MPI library

**Just go to the website, download the installer and install.**

Once done, you can activate the Intel compiler by using

```bash
source ~/intel/oneapi/setvars.sh
# This is the default installation directory
```