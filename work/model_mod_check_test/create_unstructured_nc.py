import numpy as np
from netCDF4 import num2date, date2num, Dataset

nc_unstruct = "unstructured.nc"

root_nc_unstruct = Dataset(nc_unstruct, 'w')

xloc = np.array([1,2,3,4], dtype=float)
yloc = np.array([6,8,10,13], dtype=float)
zloc = np.array([23,35], dtype=float)
state_time = np.array([10000,20000,50000,70000], dtype=float)/86400.

nxloc, nyloc, nzloc, ntime = len(xloc), len(yloc), len(zloc), len(state_time)
nloc = nxloc * nyloc * nzloc

temperature = np.random.random([ntime, nloc]) * 8.
flowrate    = np.random.random([ntime, nloc]) * 4.

###############################
# Read attributes and values from nc_struct file
###############################
# Create the dimensions
loc_d    = root_nc_unstruct.createDimension('location', nxloc*nyloc*nzloc)
state_time_d   = root_nc_unstruct.createDimension('state_time', ntime)
time_d   = root_nc_unstruct.createDimension('time', None)
member_d = root_nc_unstruct.createDimension('member', 1)

# Create the variables for the dimensions
modeltime2 = root_nc_unstruct.createVariable('time', 'f8', ('time', ))
time2   = root_nc_unstruct.createVariable('state_time', 'f8', ('state_time', ))
member2 = root_nc_unstruct.createVariable('member', 'i8', ('member', ))
xloc2   = root_nc_unstruct.createVariable('x_location', 'f8', ('location', ))
yloc2   = root_nc_unstruct.createVariable('y_location', 'f8', ('location', ))
zloc2   = root_nc_unstruct.createVariable('z_location', 'f8', ('location', ))
time2.units    = 'days'
time2[:]       = state_time
time2.calendar = 'none'
modeltime2.units = 'days'
modeltime2[:]    = 2.3
modeltime2.calendar = 'none'
member2[:]     = 1
member2.type, time2.type = 'dimension_value', 'dimension_value'
yloc2.units , yloc2.type = 'm',               'dimension_value'
zloc2.units , zloc2.type = 'm',               'dimension_value'
xloc2.units , xloc2.type = 'm',               'dimension_value'
modeltime2.type = 'dimension_value'

# zloc_grid, yloc_grid, xloc_grid,  = np.meshgrid(zloc[:], yloc[:], xloc[:])
xloc_grid, yloc_grid, zloc_grid,  = np.meshgrid(xloc[:], yloc[:], zloc[:])
xloc2[:], yloc2[:], zloc2[:] = xloc_grid.flatten(), yloc_grid.flatten(), zloc_grid.flatten() 

# Create the physical variables
temperature2 = root_nc_unstruct.createVariable(
                    'TEMPERATURE', 'f8', ('state_time', 'location'))
temperature2.type = 'observation_value'
temperature2.unit = '[C]'
temperature2[:] = temperature

flow_flux2   = root_nc_unstruct.createVariable(
                    'FLOW_FLUX', 'f8', ('state_time', 'location'))
flow_flux2.type = 'observation_value'
flow_flux2.unit = ''
flow_flux2[:] = flowrate

root_nc_unstruct.close()