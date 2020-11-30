import numpy as np
from netCDF4 import num2date, date2num, Dataset

nc_unstruct = "unstructured_parastate.nc"

root_nc_unstruct = Dataset(nc_unstruct, 'w')

xloc_para = np.array([1.5, 2.5], dtype=float)
yloc_para = np.array([9], dtype=float)
zloc_para = np.array([24], dtype=float)
para_time = np.array([30000], dtype=float)/86400.

xloc_state = np.array([1,2,3,4], dtype=float)
yloc_state = np.array([6,8,10,13], dtype=float)
zloc_state = np.array([23,35], dtype=float)
state_time = np.array([10000,20000,50000,70000], dtype=float)/86400.

ntime_state, ntime_para = len(state_time), len(para_time)
nxloc_state, nyloc_state, nzloc_state = len(xloc_state), len(yloc_state), len(zloc_state)
nloc_state = nxloc_state * nyloc_state * nzloc_state
nxloc_para, nyloc_para, nzloc_para = len(xloc_para), len(yloc_para), len(zloc_para)
nloc_para = nxloc_para * nyloc_para * nzloc_para

temperature = np.random.random([ntime_state, nloc_state]) * 8.
flowrate    = np.random.random([ntime_para, nloc_para]) * 4.

###############################
# Read attributes and values from nc_struct file
###############################
# Create the dimensions
state_loc_d    = root_nc_unstruct.createDimension('state_location', nxloc_state*nyloc_state*nzloc_state)
para_loc_d    = root_nc_unstruct.createDimension('para_location', nxloc_para*nyloc_para*nzloc_para)
state_time_d   = root_nc_unstruct.createDimension('state_time', ntime_state)
para_time_d   = root_nc_unstruct.createDimension('para_time', ntime_para)
time_d   = root_nc_unstruct.createDimension('time', None)
member_d = root_nc_unstruct.createDimension('member', 1)

# Create the variables for the dimensions
modeltime_v = root_nc_unstruct.createVariable('time', 'f8', ('time', ))
state_time_v   = root_nc_unstruct.createVariable('state_time', 'f8', ('state_time', ))
para_time_v   = root_nc_unstruct.createVariable('para_time', 'f8', ('para_time', ))
member = root_nc_unstruct.createVariable('member', 'i8', ('member', ))
state_xloc_v   = root_nc_unstruct.createVariable('state_x_location', 'f8', ('state_location', ))
state_yloc_v   = root_nc_unstruct.createVariable('state_y_location', 'f8', ('state_location', ))
state_zloc_v   = root_nc_unstruct.createVariable('state_z_location', 'f8', ('state_location', ))
para_xloc_v   = root_nc_unstruct.createVariable('para_x_location', 'f8', ('para_location', ))
para_yloc_v   = root_nc_unstruct.createVariable('para_y_location', 'f8', ('para_location', ))
para_zloc_v   = root_nc_unstruct.createVariable('para_z_location', 'f8', ('para_location', ))

state_time_v.units    = 'days'
state_time_v[:]       = state_time
state_time_v.calendar = 'none'
para_time_v.units    = 'days'
para_time_v[:]       = para_time
para_time_v.calendar = 'none'
modeltime_v.units = 'days'
modeltime_v[:]    = 2.3
modeltime_v.calendar = 'none'
member[:]     = 1
modeltime_v.type = 'dimension_value'
member.type, state_time_v.type, para_time_v.type = 'dimension_value', 'dimension_value', 'dimension_value'
state_xloc_v.units , state_xloc_v.type = 'm',               'dimension_value'
state_yloc_v.units , state_yloc_v.type = 'm',               'dimension_value'
state_zloc_v.units , state_zloc_v.type = 'm',               'dimension_value'
para_xloc_v.units , para_xloc_v.type = 'm',               'dimension_value'
para_yloc_v.units , para_yloc_v.type = 'm',               'dimension_value'
para_zloc_v.units , para_zloc_v.type = 'm',               'dimension_value'

state_xloc_grid, state_yloc_grid, state_zloc_grid,  = np.meshgrid(xloc_state[:], yloc_state[:], zloc_state[:])
state_xloc_v[:], state_yloc_v[:], state_zloc_v[:] = state_xloc_grid.flatten(), state_yloc_grid.flatten(), state_zloc_grid.flatten() 
para_xloc_grid, para_yloc_grid, para_zloc_grid,  = np.meshgrid(xloc_para[:], yloc_para[:], zloc_para[:])
para_xloc_v[:], para_yloc_v[:], para_zloc_v[:] = para_xloc_grid.flatten(), para_yloc_grid.flatten(), para_zloc_grid.flatten() 

# Create the physical variables
temperature2 = root_nc_unstruct.createVariable(
                    'TEMPERATURE', 'f8', ('state_time', 'state_location'))
temperature2.type = 'observation_value'
temperature2.unit = '[C]'
temperature2[:] = temperature

flow_flux2   = root_nc_unstruct.createVariable(
                    'FLOW_FLUX', 'f8', ('para_time', 'para_location'))
flow_flux2.type = 'observation_value'
flow_flux2.unit = ''
flow_flux2[:] = flowrate

root_nc_unstruct.close()