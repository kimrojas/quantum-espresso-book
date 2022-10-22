#!/usr/bin/env python3
import os
import shutil
import sys
import subprocess as sp

'''
DFTB+ VERSION 21.2 (OpenMP)
SMITH INSTALLATION using INTEL COMPILERS
REQUIREMENTS:
      1. Pre-load environment modules = cmake, intel, intelmpi, python
      2. Activate the required python environment


'''

# SETTINGS
version = '22.1'
variant = 'OpenMP'
if '--no-test' in sys.argv:
    test=0
else:
    test=1


# LINKS
tarlink = f'https://github.com/dftbplus/dftbplus/releases/download/{version}/dftbplus-{version}.tar.xz'
tarfile = f'dftbplus-{version}.tar.xz'
outdir = f'dftbplus-{version}'
srcdir = f'{outdir}-{variant}'

def process(line):
    print('\033[96m'+line+'\033[0m')

# Download and extract Source code
process('DOWNLOADING FILE')
sp.run(f'wget -q -O {tarfile} {tarlink}', shell=True)
print("    "+tarfile+" downloaded")
process("EXTRACTING FILE")
sp.run(f'tar xf {tarfile}', shell=True)
print("    "+tarfile+" extracted to "+srcdir)
os.remove(tarfile)


# Prepare source directory
process("PREPARING SOURCE FILES")
os.rename(outdir, srcdir)
os.chdir("./" + srcdir)
sp.run('./utils/get_opt_externals', shell=True)

# Prepare input parameters and other essentials
# **OpenMP version parameters**
builddir = "_build"
srcdir_full = os.getcwd()
installdir = "_install"
installdir_full = os.path.join(srcdir_full,installdir)
COMPILER_OPT = 'FC=mpiifort  CC=mpiicc'
python_opt = '-DENABLE_DYNAMIC_LOADING=1 -DWITH_PYTHON=1 -DBUILD_SHARED_LIBS=1 -DWITH_API=1'
# python_opt = '-DWITH_PYTHON=1 -DWITH_API=1'
ase_opt = '-DWITH_SOCKETS=1'

# python_opt = ''
# ase_opt = ''


CMAKE_OPT = f"-DCMAKE_INSTALL_PREFIX={installdir_full} {python_opt} {ase_opt} -DTEST_OMP_THREADS=2"

# Rebuild build directory
process("REBUILDING BUILD DIRECTORY")
loc = os.getcwd()
builddir_full = os.path.join(loc,builddir)
shutil.rmtree(builddir_full, ignore_errors=True)
os.makedirs(builddir)
print("    "+"Rebuilding successful:   "+builddir_full)

# CMAKE configuration protocol
process("CMAKE CONFIGURATION PROTOCOL")
command = f"{COMPILER_OPT} cmake {CMAKE_OPT} -B {builddir} ./"
logfile = open('config.log', 'w')
proc = sp.Popen(command, stdout=sp.PIPE, stderr=sp.STDOUT, universal_newlines=True, shell=True)
for line in proc.stdout:
    sys.stdout.write("    "+line)
    logfile.write(line)
proc.wait()
logfile.close()


# CMAKE build protocol
process("CMAKE BUILD PROTOCOL")
command = f"cmake --build {builddir} -- -j"
logfile = open('build.log', 'w')
proc = sp.Popen(command, stdout=sp.PIPE, stderr=sp.STDOUT, universal_newlines=True, shell=True)
for line in proc.stdout:
    sys.stdout.write("    "+line)
    logfile.write(line)
proc.wait()
logfile.close()

# CMAKE test protocol
process("CMAKE TEST PROTOCOL")
if test:
    os.chdir(builddir)
    command = f"ctest -j8"
    logfile = open('test.log', 'w')
    proc = sp.Popen(command, stdout=sp.PIPE, stderr=sp.STDOUT, universal_newlines=True, shell=True)
    for line in proc.stdout:
        sys.stdout.write("    "+line)
        logfile.write(line)
    proc.wait()
    logfile.close()
    os.chdir('../')
else:
    print('Skipping tests ...')


# CMAKE install protocol
process("CMAKE INSTALL PROTOCOL")
os.environ['DESTDIR']=''
pythonpath = f'{installdir_full}/lib/python3.8/site-packages'
command = f"PYTHONPATH={pythonpath} cmake --install {builddir}"
logfile = open('install.log', 'w')
proc = sp.Popen(command, stdout=sp.PIPE, stderr=sp.STDOUT, universal_newlines=True, shell=True)
for line in proc.stdout:
    sys.stdout.write("    "+line)
    logfile.write(line)
proc.wait()
logfile.close()

# PATHing guide
process("SHOWING IMPORTANT PATHS")
logfile = open("buildlog.path",'w')
paths = []
paths.append(f"BASE DIRECTORY = {installdir_full}")
paths.append(f"- - - - -")
paths.append(f"add to PATH            : {installdir_full}/bin")
paths.append(f"add to LD_LIBRARY_PATH : {installdir_full}/lib64")
paths.append(f"add to PYTHONPATH      : {installdir_full}/lib/python3.8/site-packages/pythonapi-0.1-py3.8.egg")

paths.append(f"set variable DFTB_LIB  :  {installdir_full}/lib64")

paths.append(f"- - - - -")

paths.append(f"SLAKOS FILES")
paths.append(f"https://dftb.org/parameters/download/all-sk-files")


logfile.write('SHOWING IMPORTANT PATHS'+'\n')
for p in paths:
    logfile.write(p+'\n')
    print(p)

logfile.close()


