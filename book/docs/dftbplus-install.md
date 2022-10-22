# DFTB+ installation

There are two ways the install DFTB+:

1. via conda-forge (easier to setup but less optimized)
2. via compilation

In this documentation, the compilation method will be discussed in detail. 

## Compilation method

The `DFTB+` executable can be compiled with 4 parallelism types:

1. Serial (non-parallelized)
2. OpenMP (parallelize via threads, shared memory)
3. MPI (parallelize via mpi, distributed memory)
4. OpenMP-MPI hybrid (parallelize via thread and mpi, not discussed)

For simplicity and stability **OpenMP** parallelism will be used. This is also the default parallelism of `DFTB+`

```{admonition} What to use: Serial, OpenMP or MPI
:class: tip
**Serial** is just good for testing. 

**OpenMP** seems good for relatively small systems that uses a single node.

**MPI** *may* be better for larger systems that requires more than a single node. For small systems, MPI takes longer than OpenMP. 

Between **OpenMP** and **MPI**, OpenMP parallization is the default way of using DFTB+. Its usage is also heavily documented. MPI, on the otherhand, is not. The non-MPI version supports more excited state methods, while the MPI version has better parallelism for many tasks.
```

### Load compilers and python environment

```bash
# Load compilers
module load cmake/3.18.3
module load intel/2020.2.254
module load intelmpi/2020.2.254
module load python/3.8

# Load python environment
source activate tutorial
```

### Download and use the DFTB+ installer
We will use the DFTB+ version 22.1 (Other releases are available [here](https://github.com/dftbplus/dftbplus/releases/))

We first create a `tutorial_files` directory in $HOME to keep our files clean.

#### 1. Prepare the installation directory

```bash
# Download source files
tutordir=~/tutorial_files
mkdir -p $tutordir/apps && cd $tutordir/apps
```

#### 2. Install using the automation script

Download the installer. There are two installers available: (a) OpenMP version and (b) MPI version. Feel free to choose any, but here we use OpenMP version.

```bash
wget https://raw.githubusercontent.com/kimrojas/gofee-book/master/book/files/install-dftb-openmp.py
chmod +x install-dftb-openmp.py
./install-dftb-openmp.py
```

```{admonition} Want the MPI version?
:class: tip
In case you want the MPI version, it is available here: 
[install-dftb-mpi.py](https://raw.githubusercontent.com/kimrojas/gofee-book/master/book/files/install-dftb-mpi.py)

It works and has been tested very slightly. Hard to use.
```

### Declare the path

The `install-dftb-mpi.py` installer script will return an advice on how to properly declare the PATHing of DFTB+.
In my case, it looks like this:

```
SHOWING IMPORTANT PATHS
BASE DIRECTORY = /home/krojas/tutorial_files/apps/dftbplus-22.1-OpenMP/_install
- - - - -
add to PATH            : /home/krojas/tutorial_files/apps/dftbplus-22.1-OpenMP/_install/bin
add to LD_LIBRARY_PATH : /home/krojas/tutorial_files/apps/dftbplus-22.1-OpenMP/_install/lib64
add to PYTHONPATH      : /home/krojas/tutorial_files/apps/dftbplus-22.1-OpenMP/_install/lib/python3.8/site-packages/pythonapi-0.1-py3.8.egg
set variable DFTB_LIB  :  /home/krojas/tutorial_files/apps/dftbplus-22.1-OpenMP/_install/lib64
- - - - -
SLAKOS FILES
https://dftb.org/parameters/download/all-sk-files
```

There are a few ways to declare this: (a) Update startup script (bashrc or bash_profile), (b) custom module environment, or (c) exporter script. For simplicity, let's create an exporter script.

Create an `activate_dftb.sh` file:

```bash
#!/bin/bash
installdir=/home/krojas/tutorial_files/apps/dftbplus-22.1-OpenMP/_install

export PATH=${installdir}/bin:$PATH
export LD_LIBRARY_PATH=${installdir}/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=${installdir}/lib/python3.8/site-packages/pythonapi-0.1-py3.8.egg:$PYTHONPATH
export DFTB_LIB=${installdir}/lib64
```

From now on, we just use this script to activate dftb+,

```bash
source ~/tutorial_files/apps/activate_dftb.sh
```


**INSTALLATION COMPLETE**

<hr>

## Conda method (a pre-compiled version)

A pre-compiled DFTB+ executable is available in conda-forge via conda installation. 

It can be simply downloaded by

```bash
conda install -n tutorial -c conda-forge 'dftbplus=*=mpi_openmpi_*' dftbplus-tools dftbplus-python
```

**INSTALLATION COMPLETE**