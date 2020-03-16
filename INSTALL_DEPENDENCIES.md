# Install the Required Dependencies

DART-PFLOTRAN interface requires multiple preinstalled packages. Here, I'll detail of how to install some of them.

## Install Anaconda in Unix system through command-line
The installation of anaconda on Windows and Mac is straightforward. One can refer it from [the official website](https://docs.anaconda.com/anaconda/install/). For unix systems (e.g., supercomputer), one might prefer the installation through command-line. The procedures are as follows:
- Download the installable package by using ```curl``` in a specified location  (make sure to replace the version with the one you prefer from the [archive](https://repo.continuum.io/archive)):
  ```sh
  curl -O https://repo.continuum.io/archive/Anaconda3-2019.10-Linux-x86_64.sh 
  ```
- Install the package:
  ```sh
  bash Anaconda3-2019.10-Linux-x86_64.sh
  ```

## Create conda virtual environment and install the required packages
A more detailed document of conda environment is available at the [doc](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands). Here, I'll illustrate the procedures of installing the packages used in DART-PFLOTRAN: 
- Create a virtual environment (named as ```geosci```) with ```python 3.7```:
    ```sh
    conda create --name geosci python=3.7
    ```
- Config the conda package channel to [conda-forge](https://conda-forge.org/):
    ```sh
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    ```
- Install the required dependencies:
    ```sh
    conda install numpy scipy pandas h5py f90nml matplotlib seaborn
    conda install netcdf4 netcdf-fortran
    conda install nco
    ```

## Configure the environment variables for using netcdf-fortran
Note that one can also install netcdf-fortran through apt-get in unix systems. However, if sudo is not accessible, the netcdf and its fortran/python interface can be installed through conda. To allow the usage of netcdf fortran, we need to do some configuations as follows:
- Check the locations of netcdf ```includedir```/```libdir``` and whether netcdf-fortran interace is installed by typing:
    ```sh
    nc-config --all
    ```
    , which gives us the following information if the netcdf-fortran is installed successfully:
    ```sh
    --has-fortran   -> yes
     --fc            -> /home/conda/feedstock_root/build_artifacts/netcdf-fortran_1571147053770/_build_env/bin/x86_64-conda_cos6-linux-gnu-gfortran
     --fflags        -> /home/jian449/software/anaconda3/envs/geosci/include
     --flibs         -> -L/home/jian449/software/anaconda3/envs/geosci/lib
     --has-f90       -> TRUE
     --has-f03       -> FALSE

    --prefix        -> /home/jian449/software/anaconda3/envs/geosci
    --includedir    -> /home/jian449/software/anaconda3/envs/geosci/include
    --libdir        -> /home/jian449/software/anaconda3/envs/geosci/lib
    --version       -> netCDF 4.7.1
    ```
- To make sure the linked library of netcdf-fortran is usable by a fortran compiler (e.g., ```gfortran```), one can either add the ```libdir``` to ```LD_LIBRARY_PATH``` in your bash profile or add it in the specific conda environment. Obviously, the latter is preferable. To do that, follow the following procedures --
    Enter your environment directory and create some subdirectories and files:
    ```sh
    cd $CONDA_PREFIX
    mkdir -p ./etc/conda/activate.d
    mkdir -p ./etc/conda/deactivate.d
    touch ./etc/conda/activate.d/env_vars.sh
    touch ./etc/conda/deactivate.d/env_vars.sh
    ```
    Edit ```./etc/conda/activate.d/env_vars.sh``` as follows:
    ```sh
    #!/bin/sh

    export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=/home/jian449/software/anaconda3/envs/geosci/lib:${LD_LIBRARY_PATH}
    ```
    Edit ```./etc/conda/deactivate.d/env_vars.sh``` as follows:
    ```sh
    #!/bin/sh

    export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}
    unset OLD_LD_LIBRARY_PATH
    ```

- To use it, do the following (replace ```gfortran``` with your own fortran compiler):
    ```sh
    gfortran [source_fortran_file] -o [target_file] -I/home/jian449/software/anaconda3/envs/geosci/include -lnetcdff
    ```