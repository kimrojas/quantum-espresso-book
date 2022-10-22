# Python Environment

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