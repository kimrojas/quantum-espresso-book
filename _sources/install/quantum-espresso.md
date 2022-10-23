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

Download the installer script
```bash
wget https://raw.githubusercontent.com/kimrojas/quantum-espresso-book/master/download/qe_install.sh
```

## D. Run the installer script

To run the script:

```bash
bash qe_install.sh
```

This will take some time. When the script finished installing, it will show you the path of the executables. You need to add the `export ...` line to your `~/.bashrc` file. 

In my case:

```
- - - - - - - - - - - - - - - - - - - - - -
USAGE GUIDE

add to ~/.bashrc:

   export PATH=/home/krojas/QE/qe-7.0/_install/bin:$PATH

- - - - - - - - - - - - - - - - - - - - - -
```