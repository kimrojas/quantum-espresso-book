# Installation on Smith cluster

## Using available Quantum Espresso

A compiled and optimized quantum espresso is already made available. To use it:

### 1. Add the custom modules to `~/.bashrc`
```bash
echo "module use --append /home/krojas/share/modulefiles" >> ~/.bashrc

# Refresh environment
source ~/.bashrc
```

### 2. Activate the environment
```bash
module load qe/7.0
```

---

## Compiling from source files

### 1. Download the source files

The source file releases can be found in [QE Releases](https://gitlab.com/QEF/q-e/-/tags). In this specific tutorial, we use the [QE 7.0 Release](https://gitlab.com/QEF/q-e/-/releases/qe-7.0)

```bash
# Download
wget https://gitlab.com/QEF/q-e/-/archive/qe-7.0/q-e-qe-7.0.tar.gz
# Extract
tar zxvf q-e-qe-7.0.tar.gz
```

### 2. Activate compiler environment

We need to activate the required depdendency modules.

```bash
module load cmake/3.18.3
module load intel/2020.2.254
module load intelmpi/2020.2.254
module load python/3.8
module load libxc/5.2.2
module load git/2.17
```

### 3. Build and install

Follow the following commands:

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

### 4. Add to PATH

Add the compiled executables (fancy way to say "program") to the PATH so it can be discovered by the system

```{note}
Replace `<full-QE-directory-path>` with the full path of your quantum espresso directory.
```

```bash
echo 'export PATH=<full-QE-directory-path>/q-e-qe-7.0/_install/bin:$PATH' >> ~/.bashrc

source ~/.bashrc
```




