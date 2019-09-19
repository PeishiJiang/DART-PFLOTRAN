"""Generate the PFLOTRAN input file (i.e., PFLOTRAN.in) before and during the assimilation."""

# Author: Peishi Jiang

import os
import sys
import copy
import h5py
import numpy as np

# Let's define the observation and state vector model first
# Codes are modified from Kewei's 1D thermal model
# Credit: https://github.com/pnnl-sbrsfa/DA-HEF-paper/blob/master/src/util.py
#


class TH1D:
###############################################################################
#
#                Define the model parameter class TH1D
#
# Details:
#
# TH1D.state_vector:Nd X Ne, where Nd is the number of varibles and Ne is the number
#               of realizations
# TH1D.state_vector_range: Nd X 2, the first column is the lower bound for each model
#               parameter and the second column is the upper bound for each paramter
# TH1D.z: 1-D array containing the coordinates of all grid centers
#
###############################################################################
    def __init__(self,da_para,nreaz,hz,spinup_length,**kwargs):

        self.da_para = da_para
        self.hz = hz
        self.dz = 0.01
        self.z =(np.arange(-hz,0,self.dz)+np.arange(-hz+self.dz,self.dz,self.dz))/2
        self.spinup = spinup_length

        if 'permeability' in da_para:
            logperm_mean = kwargs['logperm_mean']
            logperm_sd = kwargs['logperm_sd']
            logperm_low_bound = kwargs['logperm_low_bound']
            logperm_up_bound = kwargs['logperm_up_bound']
            init_logperm = np.random.normal(logperm_mean,logperm_sd,nreaz)
            logperm_range = np.array([logperm_low_bound,logperm_up_bound])
            init_logperm[init_logperm<logperm_range[0]] = logperm_range[0]
            init_logperm[init_logperm>logperm_range[1]] = logperm_range[1]
            init_perm = 10**(init_logperm)
        if 'thermal conductivity' in da_para:
            th_cond_mean = kwargs['th_cond_mean']
            th_cond_sd = kwargs['th_cond_sd']
            th_cond_low_bound = kwargs['th_cond_low_bound']
            th_cond_up_bound = kwargs['th_cond_up_bound']
            init_th_cond = np.random.normal(th_cond_mean,th_cond_sd,nreaz)
            th_cond_range = np.array([th_cond_low_bound,th_cond_up_bound])
            init_th_cond[init_th_cond<th_cond_range[0]] = th_cond_range[0]
            init_th_cond[init_th_cond>th_cond_range[1]] = th_cond_range[1]
        if 'porosity' in da_para:
            poro_mean = kwargs['poro_mean']
            poro_sd = kwargs['poro_sd']
            poro_low_bound = kwargs['poro_low_bound']
            poro_up_bound = kwargs['poro_up_bound']
            init_poro = np.random.normal(poro_mean,poro_sd,nreaz)
            poro_range = np.array([poro_low_bound,poro_up_bound])
            init_poro[init_poro<poro_range[0]] = poro_range[0]
            init_poro[init_poro>poro_range[1]] = poro_range[1]
        if 'flux' in da_para:
            flux_mean = kwargs['flux_mean']
            flux_sd = kwargs['flux_sd']
            flux_low_bound = kwargs['flux_low_bound']
            flux_up_bound = kwargs['flux_up_bound']
            flux_range = np.array([flux_low_bound,flux_up_bound])
            init_flux = np.random.normal(flux_mean,flux_sd,nreaz)
            init_flux[init_flux<flux_range[0]] = flux_range[0]
            init_flux[init_flux>flux_range[1]] = flux_range[1]

        if len(da_para) == 1:
            self.state_vector = np.zeros((1,nreaz))
            self.state_vector_range = np.zeros((1,2))
            if 'permeability' in da_para:
                self.state_vector = np.array([np.log10(init_perm)])
                self.state_vector_range = np.array([logperm_range])
            elif 'flux' in da_para:
                self.state_vector = np.array([init_flux])
                self.state_vector_range = np.array([flux_range])
            else:
                raise Exception("Please choose 'permeability' or 'flux'")
        elif len(da_para) == 2:
            self.state_vector = np.zeros((2,nreaz))
            if 'permeability' in da_para:
                self.state_vector = np.array([np.log10(init_perm)])
                self.state_vector_range = np.array([logperm_range])
            elif 'flux' in da_para:
                self.state_vector = np.array([init_flux])
                self.state_vector_range = np.array([flux_range])
            else:
                raise Exception("Please choose 'permeability' or 'flux'")
            if 'thermal conductivity' in da_para:
                self.state_vector = np.concatenate((self.state_vector,np.array([init_th_cond])))
                self.state_vector_range = np.concatenate((self.state_vector_range,np.array([th_cond_range])))
            elif 'porosity' in da_para:
                self.state_vector = np.concatenate((self.state_vector,np.array([init_poro])))
                self.state_vector_range = np.concatenate((self.state_vector_range,np.array([poro_range])))
            else:
                raise Exception("Please choose 'thermal conductivity' or 'porosity'")
        elif len(da_para) == 3:
            self.state_vector = np.zeros((3,nreaz))
            if 'permeability' in da_para:
                self.state_vector = np.array([np.log10(init_perm)])
                self.state_vector_range = np.array([logperm_range])
            elif 'flux' in da_para:
                self.state_vector = np.array([init_flux])
                self.state_vector_range = np.array([flux_range])
            else:
                raise Exception("Please choose 'permeability' or 'flux'")
            self.state_vector = np.concatenate((self.state_vector,np.array([init_th_cond])))
            self.state_vector = np.concatenate((self.state_vector,np.array([init_poro])))
            self.state_vector_range = np.concatenate((self.state_vector_range,np.array([th_cond_range])))
            self.state_vector_range = np.concatenate((self.state_vector_range,np.array([poro_range])))
        else:
            raise Exception("Maximum number of parameters is 3")


class Obs:
################################################################################
 #    obs: obsevation object which contains the observation time and dataself.
 #          obs.time: type: numpy.array, observation time series
 #          obs.ntime: total number of observation time points
 #          obs.value: type: numpy.array, observation data, e.g. Temperature
 #          obs.err_sd_ratio: ratio of the standard deviation of error to the observation value
 #          obs.coord: coordinates of the observation points
 #          obs.nobs: total number of observation points
################################################################################
    def __init__(self,obs_start_time,obs_end_time,obs_timestep,obs_coord,obs_error_type,obs_error,obs_data):
        self.start_time = obs_start_time
        self.end_time   = obs_end_time
        self.timestep   = obs_timestep
        self.time       = obs_data[:,0]
        self.ntime      = self.time.size
        self.value      = obs_data[:,1:]
        self.coord      = obs_coord
        self.error_type = obs_error_type
        self.error      = obs_error
        self.nobs       = obs_coord.size


def generate_dbase(nreaz,mod,filename):
################################################################################
#
#   GenerateDbase: generate h5 file Dbase.h5 after each assimlation.
#   Dbase is a keyword in PFLOTRAN that makes the scalar value realization dependent.
#   In the TH1D model, the Dbase.h5 contains two datasets, the first is Permeability
#   and the second is ThermalConductivity. This function will be called in each iteration
#   to update the parameters.
#
#   Details:
#        nreaz: number of realizations
#        mod: model-specific object, for TH1D model, it contains the permeability
#            and thermal conductivity and associated hard limits
#
################################################################################
    # filename = "./pflotran_results/Dbase.h5"
    if os.path.isfile(filename):
        os.remove(filename)

    h5file = h5py.File(filename,'w')
    variables = []
    if 'permeability' in mod.da_para:
        variables.append("PERMEABILITY")
    elif 'flux' in mod.da_para:
        # variables.append("Flux_top")
        variables.append("FLOW_FLUX")
    else:
        raise Exception("Please choose 'permeability' or 'flux'")
    if 'thermal conductivity' in mod.da_para:
        # variables.append('ThermalConductivity')
        variables.append('THERMAL_CONDUCTIVITY')
    if 'porosity' in mod.da_para:
        variables.append('POROSITY')

    values = copy.deepcopy(mod.state_vector)
    if 'permeability' in mod.da_para:
        values[0] = 10**(values[0])
#     if 'flux' in mod.da_para:
#         values[0] = values[0]

    for i in range(len(variables)):
        if h5file.get(variables[i]):
            del h5file[variables[i]]
        h5dset = h5file.create_dataset(variables[i],data=values[i])
#    mod.state_vector[0] = np.log(mod.state_vector[0])
    h5file.close()


########################################################
# Parameter settings
########################################################
# Read them from user input
pflotran_in_file   = sys.argv[1]
pflotran_para_file = sys.argv[2]
obs_timestep       = float(sys.argv[3])              # unit:s, the time interval that temperatures are collected
obs_error          = float(sys.argv[4])              # If the error type is 'absolute', the error means the accuracy of temperature measurement with unit degree C. If the error type is 'relative', the error means the percentage of temperature measurement.
nreaz              = int(sys.argv[5])                # number of ensemble members
spinup             = bool(sys.argv[6])

#configure ES-MDA
niter = 2 # number of iterations
alpha = np.array([2.0, 2.0]) # iteration coefficient for each iteration
da_para = ['flux','thermal conductivity','porosity'] # parameters to be estimated
# da_para = ['permeability'] # parameters to be estimated
da_timestep = -999 # unit: sec, data assimilation time step

# logperm_mean = -11.0          # unit: log(m2), mean of prior log(permeability)
# logperm_sd = 1.0          # unit: log(m2), S.D. of prior log(permeability)
# logperm_up_bound = -9.0     # unit: log(m2), upper bound of log(permeability)
# logperm_low_bound = -13.0   # unit: log(m2), lower bound of log(permeability)

flux_mean = 0.0          # unit: m/day, mean flux
flux_sd = 0.5          # unit: m/day, standard deviation of flux
flux_up_bound = 5.0     # unit: m/day, maximum downwelling flux
flux_low_bound = -5.0   # unit: m/day, maximum upwelling flux

th_cond_mean = 1.5          # unit: W/(mK)-1, mean of prior thermal conductivity
th_cond_sd = 0.5          # unit: W/(mK)-1, S.D. of prior thermal conductivity
th_cond_up_bound = 2.5          # unit: W/(mK)-1, upper bound of thermal conductivity
th_cond_low_bound = 0.5          # unit: W/(mK)-1, lower bound of thermal conductivity

poro_mean = 0.3          # mean porosity
poro_sd = 0.1          # standard deviation of porosity
poro_up_bound = 0.7          # upper bound of porosity
poro_low_bound = 0.01          # lower bound of porosity

# Configure observation
# obs_length = 3.4756944444444446     # unit:day, length of observation data used for flux estimation
obs_length = 1.     # unit:day, length of observation data used for flux estimation
#TODO should it be read from the observation file?
# from temperature.csv
therm_loc = [-0.01, -0.05, -0.65] # unit:m, location of thermistor, negative means below the riverbed
obs_error_type = 'absolute'    # 'absolute' and 'relative'. 'absolute' means the absolute measurement error in degree C, 'relative' means a perentage of the observation value
# Configure model domain and PFLOTRAN running environment
hz = 0.64          # unit: m, height of the 1-D column
spinup_length = .5 #unit: day, spinup time


#----------------------------------------------------------
kwargs1 = {}
if 'permeability' in da_para:
    kwargs1.update({'logperm_mean':logperm_mean})
    kwargs1.update({'logperm_sd':logperm_sd})
    kwargs1.update({'logperm_low_bound':logperm_low_bound})
    kwargs1.update({'logperm_up_bound':logperm_up_bound})

if 'flux' in da_para:
    kwargs1.update({'flux_mean':flux_mean})
    kwargs1.update({'flux_sd':flux_sd})
    kwargs1.update({'flux_low_bound':flux_low_bound})
    kwargs1.update({'flux_up_bound':flux_up_bound})

if 'thermal conductivity' in da_para:
    kwargs1.update({'th_cond_mean':th_cond_mean})
    kwargs1.update({'th_cond_sd':th_cond_sd})
    kwargs1.update({'th_cond_low_bound':th_cond_low_bound})
    kwargs1.update({'th_cond_up_bound':th_cond_up_bound})

if 'porosity' in da_para:
    kwargs1.update({'poro_mean':poro_mean})
    kwargs1.update({'poro_sd':poro_sd})
    kwargs1.update({'poro_low_bound':poro_low_bound})
    kwargs1.update({'poro_up_bound':poro_up_bound})

spinup_length_sec = spinup_length*86400
obs_length_sec = obs_length*86400
obs_start_time = spinup_length_sec
obs_end_time = obs_start_time+obs_length_sec


########################################################
# Obtain the parameters and the state vector information
########################################################
mod = TH1D(da_para,nreaz,hz,spinup_length,**kwargs1)


########################################################
# Obtain the time information
########################################################
obs_coord = np.array(therm_loc[1:-1])-np.array(therm_loc[0])
obs_data = np.loadtxt('../pflotran_input/obs_data.dat',skiprows=1)
obs_start_idx = int(obs_start_time/obs_timestep)
obs_end_idx = int(obs_end_time/obs_timestep)
obs_data = obs_data[obs_start_idx:obs_end_idx+1,:]

obs = Obs(obs_start_time,obs_end_time,obs_timestep,obs_coord,obs_error_type,obs_error,obs_data)


########################################################
# Write the PFLOTRAN input card file
########################################################
# Read the input card template
with open('../pflotran_input/1dthermal_template.in','r') as f:
    pflotranin = f.readlines()

# Write the new input card
with open(pflotran_in_file,'w') as f:
    for i,s in enumerate(pflotranin):
        if 'NXYZ' in s:
            pflotranin[i] = "  NXYZ 1 1 {}".format(int(mod.hz*100)) + "\n"
            pflotranin[i+2] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
        if 'REGION all' in s and 'COORDINATES' in pflotranin[i+1]:
            pflotranin[i+2] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
        if 'REGION bottom' in s and "FACE" in pflotranin[i+1]:
            pflotranin[i+3] = "    0.d0 0.d0 {}".format(-mod.hz) + "d0" + "\n"
            pflotranin[i+4] = "    1.d0 1.d0 {}".format(-mod.hz) + "d0" + "\n"
        if 'SNAPSHOT_FILE' in s:
            if spinup:
                pflotranin[i+1] = "   PERIODIC TIME {}".format(259200) + " s" +"\n"
            else:
                pflotranin[i+1] = "    PERIODIC TIME {}".format(obs.timestep) + " s" +"\n"
        if "FLOW_CONDITION flow_top" in s and "TYPE" in pflotranin[i+1]:
            pflotranin[i+2] = "    FLUX NEUMANN" + "\n"
            pflotranin[i+5] = "\n"
            pflotranin[i+6] = "  FLUX DBASE_VALUE FLOW_FLUX m/day" +"\n"
        if 'FLOW_CONDITION initial' in s and "TYPE" in pflotranin[i+1]:
            pflotranin[i+6] = "  TEMPERATURE " + str(np.mean(obs.value[0,:])) + "d0" + "\n"
        if 'THERMAL_CONDUCTIVITY_WET' in s:
            if 'thermal conductivity' in mod.da_para:
                pflotranin[i] = "  THERMAL_CONDUCTIVITY_WET DBASE_VALUE THERMAL_CONDUCTIVITY" + "\n"
        if 'POROSITY' in s:
            if 'porosity' in mod.da_para:
                pflotranin[i] = "  POROSITY DBASE_VALUE POROSITY" + "\n"
        if "FINAL_TIME" in s:
            pflotranin[i] = "  FINAL_TIME {} sec".format(obs.time[-1])+"\n"
        if "MODE TH" in s:
            if not spinup:
             pflotranin[i+1] = "      OPTIONS"+"\n"
             pflotranin[i+2] = "        REVERT_PARAMETERS_ON_RESTART"+"\n"
             pflotranin[i+3] = "      /"+"\n"
        if "FILENAME 1dthermal" in s:
            if not spinup:
#                        print(pflotranin[i+2])
                pflotranin[i-1] = "  RESTART"+"\n"
                pflotranin[i] = "    FILENAME 1dthermal-restart.chk "+" \n"
                pflotranin[i+1] = "    REALIZATION_DEPENDENT"+"\n"
                if obs.time[0] < 2e4:
                  pflotranin[i+2] = "#    RESET_TO_TIME_ZERO /" +"\n"
                else:
                  pflotranin[i+2] = "#    RESET_TO_TIME_ZERO /" +"\n"
        if "FINAL_TIME" in s:
            if spinup:
                pflotranin[i] = "  FINAL_TIME {} sec".format(mod.spinup*86400)+"\n"
            else:
                pflotranin[i] = "  FINAL_TIME {} sec".format(obs.time[-1])+"\n"

    # Write them into the file
    # for line in pflotranin:
        # f.write(line)
    f.writelines(pflotranin)

print("Finished generating the input card for PFLOTRAN...")


########################################################
# Generate the DBase for PFLOTRAN
########################################################
generate_dbase(nreaz, mod, pflotran_para_file)

print("Finished generating the DBASE for PFLOTRAN...")
