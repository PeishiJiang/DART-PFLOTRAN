# PFLOTRAN-DART

This is an interface package for integrating [PFLOTRAN](www.pflotran.org) and [DART](https://www.image.ucar.edu/DAReS/DART/). This objective is to allow user to conduct data assimilation on PFLOTRAN by using DART utilities.

## Prerequisites

It requires the installation of PFLOTRAN, DART, and some Python pacakges.

### Install PFLOTRAN

Please refer to PFLOTRAN's official [instruction](https://www.pflotran.org/documentation/user_guide/how_to/installation/linux.html#linux-install). For [NERSC Cori](https://nersc.gov/) users, a detailed instruction is available at [here](https://github.com/pnnl-sbrsfa/how-to-guide/blob/master/Compile-PFLOTRAN-on-Cori.md). 

### Install DART

Please refer [here](./Install_DART.md) for DART installation.

### Install Python packages

The required Python pacakge includes:

```

```

### Move the PFLOTRAN-DART repository into DART

Put the PFLOTRAN-DART repository in DART-compliant models repository by:

```
mv {PFLOTRAN-DART} {DART}/Manhattan/models/pflotran
```

## File structure

The main structure of this {PFLOTRAN-DART} repository is shown below:

```
.
+-- README.md         # The README file for a general introduction
+-- Install_DART.md   # The procedures of installing and configuring DART
+-- model_mod.F90     # The interface for linking PFLOTRAN and DART
|
+-- obs_kind/         # The repository containing the DEFAULT_obs_kind_mod.F90 file
+-- utils/            # The utility repository
+-- work/             # The repository containing shell scripts and compiling files
+-- applications/     # The application repository containing for running DART-PFLOTRAN
```

- ```model_mod.F90```: This file provides the Fortran interfaces for a minimal implementation of shaping PFLOTRAN as a DART-compliant model. A detailed introduction of the introduced Fortran interfaces can be found [here](https://www.image.ucar.edu/DAReS/DART/Manhattan/models/template/model_mod.html).

- ```work```: The folder provides a set of scripts for integrating PFLOTRAN and DART. It includes (1) shell scripts for running PFLOTRAN with DART (i.e., ```run_filter.csh``` and ```advance_model.csh```); (2) the template for input namelists file (i.e., ```input.nml.template```); (3) the shell script for converting NetCDF observation data to [DART format](https://www.image.ucar.edu/DAReS/DART/DART2_Observations.html#obs_seq_overview) (i.e., ```dart_seq_convert.csh```); (4) the shell script for [check](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/model_mod_check/model_mod_check.html) ```model_mod.F90```  (i.e., ```check_model_mod.csh```); and (5) other mkmf files and path names files required by the previous shell scripts. 

  ```
  ./work/
  +-- run_filter.csh
  +-- advance_model.csh
  +-- check_model_mod.csh
  +-- dart_seq_convert.csh
  +-- input.nml.template
  +-- mkmf_filter
  +-- mkmf_preprocess
  +-- mkmf_model_mod_check
  +-- mkmf_convert_nc
  +-- path_names_filter
  +-- path_names_preprocess
  +-- path_names_model_mod_check
  +-- path_names_convert_nc
  ```

- ```obs_kind```: The folder contains the ```DEFAULT_obs_kind_mod.F90``` defining a list of DART variable generic quantity, including PFLOTRAN's variables.

  ```
  ./obs_kind/
  +-- DEFAULT_obs_kind_mod.F90
  +-- DEFAULT_obs_kind_mod_copy.F90
  ```

- ```utils```: This folder provides a set of utility scripts for DART format observation conversion, preparing DART's prior data in NetCDF format, modifying ```DEFAULT_obs_kind_mod.F90``` by adding new PFLOTRAN's variables, and preparing the input namelist file.

  ```
  ./utils/
  +-- convert_nc_template.F90
  +-- csv2nc.py
  +-- list2dartqty.py
  +-- prepare_convert_nc.py
  +-- prepare_convertnc_nml.py
  +-- prepare_input_nml.py
  +-- prepare_pflotran_inpara.py
  +-- prepare_prior_nc.py
  ```

- ```applications``` This folder is where the applications should be put. A default ```template``` folder is provided, where the folder structure of each application should follow. An example of 1D thermal model is also provided in ```1dthermal``` folder. For running PFLOTRAN-DART, we suggest users to utilize the jupyter notebooks under ```workflow_ipynb``` to run through the work flow.

  ```
  ./applications/
  +-- workflow_ipynb/
  |   +-- DART_PFLOTRAN_integrat.ipynb
  +-- template/
  |   +-- dart_inout/
  |   +-- obs_type/
  |   +-- pflotran_input/
  |   +-- pflotran_output/
  |   +-- work/
  +-- 1dthermal/
  ```

  

## Configuration

## Usage



