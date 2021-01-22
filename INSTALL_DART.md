# Install DART for running the DART-PFLOTRAN interface

A detailed installation instruction for Steps 1 and 2 are available at [here](https://www.image.ucar.edu/DAReS/DART/DART2_Starting.html).

## Step 1: Download DART Manhattan version from the offical website

Complete the [online form](https://www2.cisl.ucar.edu/software/dart/download) and gain the instructions for downloading DART. **Note that** be sure to get the Manhattan version.



## Step 2: Configure Fortran and NetCDF in DART

Once the DART package is downloaded, in your DART directory {DART}, open the make makefile template in ```{DART}/manhattan/build_templates/mkmf.template```. Then, modify the following information:

  ```sh
  MPIFC = # the mpif90 utility
  MPILD = # the mpif90 utility
  FC = # the Fortran compiler
  LC = # the name of the loader; typically, the same as the Fortran compiler
  ```

Next, modify the location of your NetCDF installation containing ```netcdf.mod``` and ```typesizes.mod``` in the ```include``` subdirectory (if you installed netcdf-fortran using [conda virtual environment](./INSTALL_CONDA_VIRT_ENV), the ```NETCDF``` location is the path for the virtual environment, e.g., ```[parentdir_of_anaconda3]/anaconda3/envs/dart-pflotran```):

  ```shell
  NETCDF = # the location of your netCDF installation
  ```



## Step 3: Allow the permission to execute some utility functions
Enable the permission of using the following files by entering:
```shell script
chmod +x {DART}/manhattan/build_templates/mkmf
chmod +x {DART}/manhattan/assimilation_code/modules/utilities/fixsystem
```



## Step 4: Configure DART for PFLOTRAN-DART interface usage

Replace the following files in DART:
- first, replace ```{DART}/manhattan/observations/obs_converters/utilities/obs_utilities_mod.f90``` with ```{DART-PFLOTRAN}/smoother/obs_utilities_mod.f90```;
- second, replace ```{DART}/manhattan/assimilation_code/modules/utilities/types_mod.f90``` with ```{DART-PFLOTRAN}/smoother/types_mod.f90```.

The new files contains the following changes in the original DART Manhattan version to allow the usage of the DART-PFLOTRAN interface:
- Modifying ```{DART}/manhattan/observations/obs_converters/utilities/obs_utilities_mod.f90```, to allow the subroutine ```set_obs_def_location()``` to read cartesian 3D grids by replacing:

  ```fortran
  type(obs_def_type)    :: obs_def
  
  call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
  call set_obs_def_type_of_obs(obs_def, okind)
  ```

  with

  ```fortran
  type(obs_def_type)    :: obs_def
  real(r8)              :: loc_info(4)
  
  loc_info(1) = lat
  loc_info(2) = lon
  loc_info(3) = vval
  loc_info(4) = real(vkind)
  
  call set_obs_def_location(obs_def, set_location(loc_info))
  call set_obs_def_type_of_obs(obs_def, okind)
  ```

- Modifying ```{DART}/manhattan/assimilation_code/modules/utilities/types_mod.f90```, to increase the maximum allowed observation type name length by replacing:

  ```fortran
  integer, parameter :: obstypelength    = 31
  ```

  with

  ```fortran
  integer, parameter :: obstypelength    = 64
  ```
