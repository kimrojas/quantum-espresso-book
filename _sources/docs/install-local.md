# Compile on local computer

There are a few general advantages of using Quantum Espresso from a local machine such as doing basic tutorials and testing scripts for input file creations. There are two ways to install quantum espresso in your local machine:

1. via CONDA (serial or single-core executable)
2. via compilation (can be serial or parallelized executable)


## Conda method

With an activated conda:

```bash
conda create -n qe -c conda-forge  pyton=3.8 qe
```

---

## Compilation method

There are many options for the compilers, math library and MPI but for simplicity (in Intel-based system) we will use [oneAPI Toolkits](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.5fpuqp). You will need the OneAPI [base](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) and [HPC](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html) Toolkits. 

For Ryzen-based systems, you might need to use their ryzen-optimized compilers and math library (AOCC and AOCL). Installation on Ryzen-based systems will not be tackled here but you can refer to this [article](https://developer.amd.com/spack/quantum-espresso/) for more information. 

The procedure are as follows: 

### Download the source files

The source file releases can be found in [QE Releases](https://gitlab.com/QEF/q-e/-/tags). In this specific tutorial, we use the [QE 7.0 Release](https://gitlab.com/QEF/q-e/-/releases/qe-7.0)

```bash
# Download
wget https://gitlab.com/QEF/q-e/-/archive/qe-7.0/q-e-qe-7.0.tar.gz
# Extract
tar zxvf q-e-qe-7.0.tar.gz
```

### Activate compiler environment

We need to activate the required depdendency modules.

:::{note}
Replace the `<oneapi-intel-path>` with the path to your `setvars.sh` script.
:::

```bash
source <oneapi-intel-path>/setvars.sh
```

### Build and install

Follow the following commands:

:::{note}
Replace `make -j8` with `make -jN` where `N` is your allowed cpu to be used. For example, if my computer have 4 CPUs only, then `N`=4.
:::

```bash
# Initialize build directory
mkdir q-e-qe-7.0/_build
cd q-e-qe-7.0/_build

# Build
cmake \
    -DQE_ENABLE_MPI=ON \
    -DQE_ENABLE_TEST=ON \
    -DQE_ENABLE_SCALAPACK=ON \
    -DQE_FFTW_VENDOR=Intel_DFTI \
    -DCMAKE_C_COMPILER=mpiicc \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_INSTALL_PREFIX=../_install \
    -DQE_ENABLE_LIBXC=ON \
    ../

# Compile
make -j8

# Install
make install
```

### Add to PATH

Add the compiled executables (fancy way to say "program") to the PATH so it can be discovered by the system

```{note}
Replace `<full-QE-directory-path>` with the full path of your quantum espresso directory.
```

```bash
echo 'export PATH=<full-QE-directory-path>/q-e-qe-7.0/_install/bin:$PATH' >> ~/.bashrc

source ~/.bashrc
```