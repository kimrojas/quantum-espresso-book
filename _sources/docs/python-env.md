# Python environment

This part can be considered as the most difficult part due to latest version of dependencies breaking compatability.

## Creating the python environment

You can create the python environment using conda's environment solver and yaml file. Conda creation simply uses more updated packages (here I set some version to preserve compatability). The yaml file is an exact way of installing all packages (more stable way of copying my working environment).


### VIA Conda creation

Feel free to change the envname.

```bash
# Load compilers
module load python/3.8

envname='tutorial'

conda create -c conda-forge -n $envname python=3.8.3 ase=3.22.1 cymem cython mpi4py numpy=1.22.1 pytest scikit-learn scipy dscribe matplotlib
```

```{note}
GOFEE calls some modules of ASE that are already depreciated with the recent version (3.22.1), hence, you will encounter depreciation warnings such as `("VisibleDeprecationWarning: Use get_global_number_of_atoms() instead")`. This particular `warning` will be printed out in the standard output. If this bothers your workflow, you can set to ignore this warning, simply add the following lines to the start of your python file.

```python
import warning
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
```

<hr>

### VIA exact `YAML` specifications

Use:

```bash
conda env create -f environment.yml
```

**The `environment.yml` YAML file**

```yaml
name: tutorial
channels:
  - conda-forge
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - ase=3.22.1=pyhd8ed1ab_1
  - attrs=21.4.0=pyhd8ed1ab_0
  - brotli=1.0.9=h166bdaf_7
  - brotli-bin=1.0.9=h166bdaf_7
  - ca-certificates=2022.6.15=ha878542_0
  - certifi=2022.6.15=py38h578d9bd_0
  - click=8.1.3=py38h578d9bd_0
  - cycler=0.11.0=pyhd8ed1ab_0
  - cymem=2.0.6=py38hfa26641_3
  - cython=0.29.30=py38hfa26641_0
  - dbus=1.13.6=hfdff14a_1
  - dscribe=1.2.1=py38h1fd1430_1
  - expat=2.4.8=h27087fc_0
  - flask=2.1.3=pyhd8ed1ab_0
  - fontconfig=2.14.0=h8e229c2_0
  - fonttools=4.34.4=py38h0a891b7_0
  - freetype=2.10.4=h0708190_1
  - gettext=0.19.8.1=hf34092f_1004
  - giflib=5.2.1=h36c2ea0_2
  - glib=2.66.3=h58526e2_0
  - gst-plugins-base=1.14.5=h0935bb2_2
  - gstreamer=1.14.5=h36ae1b5_2
  - icu=64.2=he1b5a44_1
  - importlib-metadata=4.11.4=py38h578d9bd_0
  - iniconfig=1.1.1=pyh9f0ad1d_0
  - itsdangerous=2.1.2=pyhd8ed1ab_0
  - jinja2=3.1.2=pyhd8ed1ab_1
  - joblib=1.1.0=pyhd8ed1ab_0
  - jpeg=9e=h166bdaf_2
  - kiwisolver=1.4.3=py38h43d8883_0
  - lcms2=2.12=hddcbb42_0
  - ld_impl_linux-64=2.36.1=hea4e1c9_2
  - lerc=3.0=h9c3ff4c_0
  - libblas=3.9.0=15_linux64_openblas
  - libbrotlicommon=1.0.9=h166bdaf_7
  - libbrotlidec=1.0.9=h166bdaf_7
  - libbrotlienc=1.0.9=h166bdaf_7
  - libcblas=3.9.0=15_linux64_openblas
  - libclang=9.0.1=default_hb4e5071_5
  - libdeflate=1.12=h166bdaf_0
  - libffi=3.2.1=he1b5a44_1007
  - libgcc-ng=12.1.0=h8d9b700_16
  - libgfortran-ng=12.1.0=h69a702a_16
  - libgfortran5=12.1.0=hdcd56e2_16
  - libglib=2.66.3=hbe7bbb4_0
  - libgomp=12.1.0=h8d9b700_16
  - libiconv=1.16=h516909a_0
  - liblapack=3.9.0=15_linux64_openblas
  - libllvm11=11.1.0=hf817b99_3
  - libllvm9=9.0.1=default_hc23dcda_7
  - libopenblas=0.3.20=pthreads_h78a6416_0
  - libpng=1.6.37=h753d276_3
  - libstdcxx-ng=12.1.0=ha89aaad_16
  - libtiff=4.4.0=hc85c160_1
  - libuuid=2.32.1=h7f98852_1000
  - libwebp=1.2.2=h3452ae3_0
  - libwebp-base=1.2.2=h7f98852_1
  - libxcb=1.13=h7f98852_1004
  - libxkbcommon=0.10.0=he1b5a44_0
  - libxml2=2.9.10=hee79883_0
  - libzlib=1.2.12=h166bdaf_2
  - llvmlite=0.38.1=py38h38d86a4_0
  - lz4-c=1.9.3=h9c3ff4c_1
  - markupsafe=2.1.1=py38h0a891b7_1
  - matplotlib=3.5.2=py38h578d9bd_0
  - matplotlib-base=3.5.2=py38h826bfd8_0
  - mpi=1.0=mpich
  - mpi4py=3.1.3=py38h97ac3a3_1
  - mpich=4.0.2=h846660c_100
  - munkres=1.1.4=pyh9f0ad1d_0
  - ncurses=6.3=h27087fc_1
  - nspr=4.32=h9c3ff4c_1
  - nss=3.78=h2350873_0
  - numba=0.55.2=py38hdc3674a_0
  - numpy=1.22.1=py38h6ae9a64_0
  - openjpeg=2.4.0=hb52868f_1
  - openssl=1.1.1q=h166bdaf_0
  - packaging=21.3=pyhd8ed1ab_0
  - pcre=8.45=h9c3ff4c_0
  - pillow=9.2.0=py38h0ee0e06_0
  - pip=22.1.2=pyhd8ed1ab_0
  - pluggy=1.0.0=py38h578d9bd_3
  - pthread-stubs=0.4=h36c2ea0_1001
  - py=1.11.0=pyh6c4a22f_0
  - pybind11=2.9.2=py38h43d8883_1
  - pybind11-global=2.9.2=py38h43d8883_1
  - pyparsing=3.0.9=pyhd8ed1ab_0
  - pyqt=5.12.3=py38ha8c2ead_3
  - pytest=7.1.2=py38h578d9bd_0
  - python=3.8.3=cpython_he5300dc_0
  - python-dateutil=2.8.2=pyhd8ed1ab_0
  - python_abi=3.8=2_cp38
  - qt=5.12.5=hd8c4c69_1
  - readline=8.1.2=h0f457ee_0
  - scikit-learn=1.1.1=py38hf80bbf7_0
  - scipy=1.8.1=py38h1ee437e_0
  - setuptools=63.2.0=py38h578d9bd_0
  - six=1.16.0=pyh6c4a22f_0
  - sparse=0.13.0=pyhd8ed1ab_0
  - sqlite=3.39.1=h4ff8645_0
  - threadpoolctl=3.1.0=pyh8a188c0_0
  - tk=8.6.12=h27826a3_0
  - tomli=2.0.1=pyhd8ed1ab_0
  - tornado=6.2=py38h0a891b7_0
  - unicodedata2=14.0.0=py38h0a891b7_1
  - werkzeug=2.1.2=pyhd8ed1ab_1
  - wheel=0.37.1=pyhd8ed1ab_0
  - xorg-libxau=1.0.9=h7f98852_0
  - xorg-libxdmcp=1.1.3=h7f98852_0
  - xz=5.2.5=h516909a_1
  - zipp=3.8.0=pyhd8ed1ab_0
  - zlib=1.2.12=h166bdaf_2
  - zstd=1.5.2=h8a70e8d_2
  - pip:
    - pyqt5-sip==4.19.18
    - pyqtchart==5.12
    - pyqtwebengine==5.12.1
```

<hr>


## General information on the requirements

### GOFEE

1. python (3.6.3 or newer)
2. ASE (3.17 or newer)
3. cython (0.28 or newer)
4. cymem (1.31.2 or newer)
5. mpi4py (3.0 or newer)

### DFTB+
1. numpy (any except 1.23.1)
2. ASE

```{warning}
NUMPY version 1.23.1 breaks the ASE and DFTB+ compatability which is a core process in the calculations.
```

### Custom requirements

The following requirements are tested (locally and on smith). The versions are hard set for future reproducability.

1. python=3.8
2. ase=3.22.1
3. cymem=2.0.6
4. cython=0.29.30
5. mpi4py=3.1.3
6. numpy=1.22.1
7. pytest=7.1.2
8. scikit-learn=1.1.1
9. scipy=1.8.1
10. dscribe=1.2.1
11. matplotlib=3.5.2


