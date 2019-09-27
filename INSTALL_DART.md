# Install DART for running the PFLOTRAN-DART interface

A detailed installation instruction for Steps 1 and 2 are available at [here](https://www.image.ucar.edu/DAReS/DART/DART2_Starting.html).

## Step 1: Download DART-Manhattan from the offical website

Complete the [online form](https://www2.cisl.ucar.edu/software/dart/download) and gain the instructions for downloading DART. **Note that** be sure to get the Manhattan version.



## Step 2: Configure Fortran and NetCDF in DART

Once the DART package is downloaded, in your DART directory {DART}, open the make makefile template in ```{DART}/Manhattan/build_templates/mkmf.template```. Then, modify the following information:

- the locations of your Fortran compiler:

  ```sh
  FC = # the Fortran compiler
  LC = # the name of the loader; typically, the same as the Fortran compiler
  ```

- the location of your NetCDF installation containing ```netcdf.mod``` and ```typesizes.mod```:

  ```shell
  NETCDF = # the location of your netCDF installation
  ```



## Step 3: Configure DART for PFLOTRAN-DART interface usage

Fix the following bugs in DART Manhattan version to allow the usage of the PFLOTRAN-DART interface:

- Modify ```{DART}/Manhattan/observations/obs_converters/utilities/obs_utilities_mod.f90```, to allow the subroutine ```set_obs_def_location()``` to read cartesian 3D grids by replacing:

  ```fortran
  type(obs_def_type)    :: obs_def
  
  call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
  call set_obs_def_type_of_obs(obs_def, okind)
  ```

  with

  ```fortran
  type(obs_def_type)    :: obs_def
  real(r8)              :: loc_info(4)
  
  loc_info(1) = lon
  loc_info(2) = lat
  loc_info(3) = vval
  loc_info(4) = real(vkind)
  
  call set_obs_def_location(obs_def, set_location(loc_info))
  call set_obs_def_type_of_obs(obs_def, okind)
  ```

- Modify ```{DART}/Manhattan/assimilation_code/modules/utilities/types_mod.f90```, to increase the maximum allowed observation type name length by replacing:

  ```fortran
  integer, parameter :: obstypelength    = 31
  ```

  with

  ```fortran
  integer, parameter :: obstypelength    = 64
  ```