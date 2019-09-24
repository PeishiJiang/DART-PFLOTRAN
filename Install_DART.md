# Install DART for running the PFLOTRAN-DART interface

A detailed installation instruction for Steps 1 and 2 are available at [here](https://www.image.ucar.edu/DAReS/DART/DART2_Starting.html).

## Step 1: Download DART-Manhattan from the offical website

Complete the [online form](https://www2.cisl.ucar.edu/software/dart/download) and gain the instructions for downloading DART. **Note that** be sure to get the Manhattan version.

##Step 2: Configure Fortran and NetCDF in DART

Once the DART package is downloaded, in your DART directory {DART}, open the make makefile template in ```{DART}/Manhattan/build_templates/mkmf.template```. Then, modify the following information:

- the locations of your Fortran compiler:

  ```
  FC = # the Fortran compiler
  LC = # the name of the loader; typically, the same as the Fortran compiler
  ```

- the location of your NetCDF installation containing ```netcdf.mod``` and ```typesizes.mod```:

  ```
  NETCDF = # the location of your netCDF installation
  ```

## Step 3: Configure DART for PFLOTRAN-DART interface usage

Fix the following bugs in DART Manhattan version to allow the usage of the PFLOTRAN-DART interface:

- change...
- Change...