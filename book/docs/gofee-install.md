# GOFEE Installation

## Installation procedure

### Download the file

```bash
cd ~/tutorial_files/apps
wget http://grendel-www.cscaa.dk/mkb/_downloads/694ca887012ad412602f2a45ee7b2ea2/gofee_stable.tar.gz
tar zxvf gofee_stable.tar.gz
```

### Build GOFEE package

```bash
cd gofee-master
./build_code
```

### Create an activation file

Create a gofee activator `activate_gofee.sh` with the following contents:

```bash
#!/bin/bash
installdir=/home/krojas/tutorial_files/apps/gofee-master
export PYTHONPATH=${installdir}:$PYTHONPATH
```
