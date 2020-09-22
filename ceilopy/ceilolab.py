import xarray as xr
import pandas as pd
import netCDF4
import numpy as np


def read_L1(path2file):
    
    def check_shape(var, shape):
        isshape = np.unique(var[:]).shape
        assert(isshape == shape), f'I assumed this variable had a dimension of {shape}. Actual shape is {isshape}'
        return
    
    def get_ncattr(var):
        atr = {}
        for nca in var.ncattrs():
            atr[nca] = var.getncattr(nca)
        return atr
    
    nc = netCDF4.Dataset(path2file)

    # get datetime index from the date_stamp variable ... is easier than the time variable
    dt = nc.variables['date_stamp'][:]
    dt = pd.to_datetime(dt[:,0])
    dt.name = 'datetime'

    # range coordinate
    # range has no attributes not even the unit :-|
    var = nc.variables['range']
    coord_range = pd.Index(var[:])
    coord_range.name = 'range'

    # variables with constant value
    # non has an attribute
    var_list = ['name', 'message_type', 'version',  'period', 'profile_scale', 'profile_resolution', 'profile_length']
    variables = []
    for vn in var_list:
        var = nc.variables[vn]
        check_shape(var,(1,))
        data = var[:][0,0]
        attrs = get_ncattr(var)
        variables.append(dict(name = vn, data = data, attrs = attrs))

    # tilt_angle
    # no attrs
    vn = 'tilt_angle'
    var = nc.variables[vn]
    data = pd.Series(var[:][:,0], index = dt)
    attrs =  get_ncattr(var)
    variables.append(dict(name = vn, data = data, attrs = attrs))

    # cloud status
    # no attr, so I think this gives info on, 0: no cb detected, 1: cb detected, 2: mulitple cb detected .... but  I am not sure
    vn = 'cloud_status'
    var = nc.variables[vn]
    check_shape(var, (3,))
    data = pd.Series(var[:][:,0], index = dt)
    attrs =  get_ncattr(var)
    variables.append(dict(name = vn, data = data, attrs = attrs))

    # cloud data
    # this looks like the cloud base for the 3 layers the ceilometer is able to detect
    # no attrs
    vn = 'cloud_data'
    var = nc.variables[vn]
    coord_cloud_layer = pd.Index([1,2,3])
    coord_cloud_layer.name = 'cloud_layer'
    data = pd.DataFrame(var[:], columns=coord_cloud_layer, index=dt)
    attrs =  get_ncattr(var)
    variables.append(dict(name = vn, data = data, attrs = attrs))

    # status_bits
    # no ncattrs ... hopfully the manual wil help with that
    vn = 'status_bits'
    var = nc.variables[vn]
    data = pd.Series(var[:][:,0], index = dt)
    attrs =  get_ncattr(var)
    variables.append(dict(name = vn, data = data, attrs = attrs))

    # rcs_910
    # no attrs ... we really need to know the units here ... I think
    vn = 'rcs_910'
    var = nc.variables[vn]
    data = pd.DataFrame(var[:], columns=coord_range, index=dt)
    attrs =  get_ncattr(var)
    variables.append(dict(name = vn, data = data, attrs = attrs))

    # put everything into a xarray dataset
    ds = xr.Dataset()

    for var in variables:
        ds[var['name']] = var['data']
        ds[var['name']].attrs = var['attrs']

    return ds