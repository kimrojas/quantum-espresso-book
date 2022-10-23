# Build Quantum Espresso

## A. Prepare environment

If you followed the previous steps closely, then both **intel compilers and python environment** is activated already. If not, activate it

```{admonition} Activate required environment
:class: tip

To activate the Intel package, simply run:
`source ~/intel/oneapi/setvars.sh`

To activate the python environment, simply run:
`source activate qe` or `conda activate qe`
```

## B. Check the requirements

We do a final check if the requirements are available.

```bash
which mpirun
which mpiifort
which mpiicc
which cmake
which git
```

## C. Download the installer script

For your convenience, I have prepared an installer script that will handle all the installation procedures. It is currently set to install [**Quantum Espresso v7.0**](https://gitlab.com/QEF/q-e/-/releases/qe-7.0)

```bash
wget https://raw.githubusercontent.com/kimrojas/quantum-espresso-book/master/download/qe_install.sh
```


## A. CONDA installation

1. Check if conda is available

   ```bash
   conda --version
   ```

**if conda was not found, proceed to step 2 (SKIP IF CONDA IS AVAILABLE)**

2. Install Miniconda [WEBLINK](https://docs.conda.io/en/latest/miniconda.html)
   1. Visit the website
   2. Download the installer for your system
   3. Install (choose `initialize` when asked)

3. Recheck if conda was installed successfully

   ```bash
   conda --version
   ```

## B. Create the environment

1. Download the environment information file (`qe_environment.yml`)

   ```bash
   wget https://raw.githubusercontent.com/kimrojas/quantum-espresso-book/master/download/qe_environment.yml
   ```

2. Create the environment using the file

   ```bash
   conda env create -f qe_environment.yml
   ```

   ```{admonition} QE Environment Activation?
   :class: tip

   To activate the environment, simply run:

   `source activate qe`
   or
   `conda activate qe`

   ```