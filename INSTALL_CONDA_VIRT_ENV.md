# Install Python packages and netCDF fortran interfaces using Conda virtual environment

Here, we recommend using conda virtual environment to install the python packages and netCDF fortran interfaces in linux system.

## Step 1 Install Conda (see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) for more details)
First, download conda to the directory you will install conda using ```wget``` (please replace ```2020.11``` with your preferred version):
```sh
wget https://repo.continuum.io/archive/Anaconda3-2020.11-Linux-x86_64.sh
```
Then, install conda:
```sh
bash Anaconda-2020.11-Linux-x86_64.sh
```
Once the conda is installed locally, set the package channel to [```conda-forge```](https://conda-forge.org/):
```sh
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## Step 2 Create a new virtual environment (see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for more details)
To create an environment with Python 3.7 (please replace ```2020.11``` with your preferred name), use:
```sh
conda create --name dart-pflotran python=3.7
```
When conda asks you to proceed, type ```y```:
```sh
proceed ([y]/n)?
```
Once the virtual environment is created, to activate this environment, use:
```sh
conda activate dart-pflotran
```
To deactivate an active environment, use:
```sh
conda deactivate
```

## Step 3 Install netCDF fortran interface and NCO utility
To install netCDF fortran interface, use (make sure your virtual environment is activated if you are using it):
```sh
conda install netcdf-fortran nco
```
For pinklady user, please install the following specific version of netcdf-fortran (to be compatible with gfortran on pinklady):
```sh
conda install netcdf-fortran=4.5.2 nco
```
Then, revise or add ```LD_LIBRARY_PATH``` to your virtual environment (so that when you run ```conda activate dart-pflotran```, the environment variables are set to the values you wrote into the file; and when you run ```conda deactivate```, those variables are erased.): 
```sh
cd $CONDA_PREFIX # go to the path of the current virtual environment
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
```
Edit ```./etc/conda/activate.d/env_vars.sh``` as follows:
```sh
#!/bin/sh

export OLD_LD_LIBRARY=${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH}
```
Edit ```./etc/conda/deactivate.d/env_vars.sh``` as follows:
```sh
#!/bin/sh

unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${OLD_LD_LIBRARY}
```
Once that is done, deactivate and activate the virtual environment again so that what have been written in ```env_vars.sh``` would be loaded.


## Step 4 Install Python packages
To install the required Python packages, use (make sure your virtual environment is activated if you are using it):
```sh
conda install numpy scipy f90nml h5py pandas netcdf4 jupyter matplotlib
```