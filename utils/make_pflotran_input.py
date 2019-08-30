"""Generate the PFLOTRAN input file (i.e., PFLOTRAN.in) before and during the assimilation."""

# Author: Peishi Jiang

import numpy as np

# Let's define the observation and state vector model first
# Credit: https://github.com/pnnl-sbrsfa/DA-HEF-paper/blob/master/src/util.py
#

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
    def __init__(self,obs_start_time,obs_end_time,obs_timestep,obs_coord,obs_error_type,obs_error,obs_data,obs_flow_dir):
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
        self.flow_dir   = obs_flow_dir



# Obtain the time information
obs =

# Obtain the parameters and the state vector information
mod =

with open('./pflotran_results/1dthermal.in','r+') as f:
    pflotranin = f.readlines()
    if with_head:
        for i,s in enumerate(pflotranin):
            if "PERM_ISO" in s:
                pflotranin[i] = "    PERM_ISO DBASE_VALUE Permeability" + "\n"
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
                    pflotranin[i+1] = "   PERIODIC TIME {}".format(obs.timestep) + " s" +"\n"
            if 'FLOW_CONDITION flow_top' in s and "TYPE" in pflotranin[i+1]:
                pflotranin[i+5] = "  DATUM FILE ../pflotran_inputs/head_top.dat" + "\n"
            if 'FLOW_CONDITION flow_bottom' in s and "TYPE" in pflotranin[i+1]:
                pflotranin[i+5] = "  DATUM FILE ../pflotran_inputs/head_bottom.dat" + "\n"
            if 'FLOW_CONDITION initial' in s and "TYPE" in pflotranin[i+1]:
                pflotranin[i+6] = "  TEMPERATURE " + str(np.mean(obs.value[0,:])) + "d0" + "\n"
            if 'THERMAL_CONDUCTIVITY_WET' in s:
                if 'thermal conductivity' in mod.da_para:
                    pflotranin[i] = "  THERMAL_CONDUCTIVITY_WET DBASE_VALUE ThermalConductivity"+ "\n"
            if 'POROSITY' in s:
                if 'porosity' in mod.da_para:
                    pflotranin[i] = "  POROSITY DBASE_VALUE Porosity"
            if "FILENAME 1dthermal" in s:
                if not spinup:
                    pflotranin[i-1] = "  RESTART"+"\n"
                    pflotranin[i] = "    FILENAME 1dthermal-restart.chk "+"0.d0 \n"
                    pflotranin[i+1] = "    REALIZATION_DEPENDENT"+"\n"
                    pflotranin[i+2] = "    RESET_TO_TIME_ZERO" +"\n"
            if 'FINAL_TIME' in s:
                if spinup:
                    pflotranin[i] = "  FINAL_TIME " + str(86400) + "  sec" + "\n"
                else:
                    pflotranin[i] = "  FINAL_TIME " + str(obs.time[-1]) + "  sec" + "\n"
