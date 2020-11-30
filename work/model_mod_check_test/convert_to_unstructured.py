import numpy as np
from netCDF4 import num2date, date2num, Dataset

nc_struct = "prior_ensemble_template.nc"
nc_unstruct = "prior_ensemble_template_unstructured.nc"

root_nc_struct   = Dataset(nc_struct, 'r')
root_nc_unstruct = Dataset(nc_unstruct, 'w')


###############################
# Read attributes and values from nc_struct file
###############################
# global attributes
model_time = root_nc_struct.model_time
time_unit  = root_nc_struct.time_unit
assim_window = root_nc_struct.assimilation_window

# dimensions
nxloc = root_nc_struct.dimensions['x_location'].size
nyloc = root_nc_struct.dimensions['y_location'].size
nzloc = root_nc_struct.dimensions['z_location'].size
ntime = root_nc_struct.dimensions['time'].size
member = root_nc_struct.dimensions['member'].size

# variables
time = root_nc_struct.variables['time']
member = root_nc_struct.variables['member']
xloc = root_nc_struct.variables['x_location']
yloc = root_nc_struct.variables['y_location']
zloc = root_nc_struct.variables['z_location']
temperature = root_nc_struct.variables['TEMPERATURE']
flow_flux   = root_nc_struct.variables['FLOW_FLUX']

###############################
# Read attributes and values from nc_struct file
###############################
# Create the dimensions
loc_d    = root_nc_unstruct.createDimension('location', nxloc*nyloc*nzloc)
state_time_d   = root_nc_unstruct.createDimension('state_time', None)
# time_d   = root_nc_unstruct.createDimension('time', 1)
time_d   = root_nc_unstruct.createDimension('time', None)
member_d = root_nc_unstruct.createDimension('member', 1)

# Create the variables for the dimensions
modeltime2 = root_nc_unstruct.createVariable('time', 'f8', ('time', ))
time2   = root_nc_unstruct.createVariable('state_time', 'f8', ('state_time', ))
member2 = root_nc_unstruct.createVariable('member', 'i8', ('member', ))
xloc2   = root_nc_unstruct.createVariable('x_location', 'f8', ('location', ))
yloc2   = root_nc_unstruct.createVariable('y_location', 'f8', ('location', ))
zloc2   = root_nc_unstruct.createVariable('z_location', 'f8', ('location', ))
time2.units    = time.units
time2[:]       = time[:]
time2.calendar = 'none'
modeltime2.units = time.units
modeltime2[:]    = time[:]
modeltime2.calendar = 'none'
member2[:]     = member[:]
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
                    temperature.name, 'f8', ('state_time', 'location'))
temperature2.type = temperature.type
temperature2.unit = temperature.unit
temperature2[:] = np.expand_dims(temperature[:][0,:,:,:].T.flatten(), axis=0)

flow_flux2   = root_nc_unstruct.createVariable(
                    flow_flux.name, 'f8', ('state_time', 'location'))
flow_flux2.type = flow_flux.type
flow_flux2.unit = flow_flux.unit
flow_flux2[:] = np.expand_dims(flow_flux[:][0,:,:,:].T.flatten(), axis=0)

root_nc_struct.close()
root_nc_unstruct.close()