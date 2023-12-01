#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 12:12:26 2023

@author: hagen
"""
import pathlib as _pl
import numpy as _np
import xarray as _xr
import pandas as _pd
import magic as _magic

import ceilopy.ceilolab as cl

class CorruptFileError(Exception):
    """Exception raised when File is not whats expected.
    """
    def __init__(self, message):
        super().__init__(message)
        
class MissingSerialNumberError(Exception):
    """Exception raised when File does not contain Serial number.
    """
    def __init__(self, message):
        super().__init__(message)
        
class SerialNumberMissmatchError(Exception):
    """Exception raised when Files doe not have the same serial number.
    """
    def __init__(self, message):
        super().__init__(message)
        
class IncompleteProfileError(Exception):
    pass

#### Warnings
class IncompleteProfileWarning(Warning):
    pass

# def read_netcdf_level1(file, 
#             parent = None, 
#             ignore1 = ['name','message_type','version',
#                        'date_stamp',
#                        'period','tilt_angle',
#                       'cloud_status','cloud_data','status_bits','profile_scale',
#                       'profile_resolution','profile_length'],
#             time_shift_error = 'raise'):
#     if isinstance(file, (str, _pl.Path)):
#         file = [file]
        
#     assert(isinstance(file, (_pd.Series,list, _np.array))), f'File type not recognized: {type(file)}'
    
#     # ignore1 = ['name','message_type','version',
#     #            'date_stamp',
#     #            'period','tilt_angle',
#     #           'cloud_status','cloud_data','status_bits','profile_scale',
#     #           'profile_resolution','profile_length']
    
#     if not _np.all([_magic.from_file(fn.as_posix()) == 'Hierarchical Data Format (version 5) data' for fn in file]):
#         fnc = '\n\t'.join([fn.as_posix() for fn in file])
#         raise CorruptFileError(f'At least one of the following can not be identified as a netcdf file: \n\t {fnc}')
#     # parent.tp_file = file.copy()    
#     # cause problems:L1 = _xr.open_mfdataset(file, drop_variables=ignore1) # , concat_dim = 'timeDim' the usage of this kwarg seems to have changed and is not needed anymore
#     L1 = _xr.open_mfdataset(file, drop_variables=ignore1, concat_dim = 'timeDim', combine='nested')
#     # return L1
#     L1 = L1.assign_coords(time = _pd.to_datetime(L1.time.values, unit = 's'))
#     # parent.tp_ = 
#     for var in L1.variables:
#         if 'timeDim' in L1[var].dims:
#             L1[var] = L1[var].swap_dims({'timeDim':'time'})
        
#     if time_shift_error == 'raise':
#         tdiff = L1.time.values[1:] - L1.time.values[:-1]
#         tdiffmin = tdiff.min()/_pd.to_timedelta(1,'s')
        
#         excepted_mindiff = -60
#         assert(tdiffmin > excepted_mindiff), f'There is a jump in the time of {tdiffmin} s, which is lareger than the excepted value of {excepted_mindiff}. Probably due to time update? Programming required'
#         if 0:
#             # the maxdiff is problematic, as there are cases where large chuncs of data is missing
#             tdiffmax = tdiff.max()/_pd.to_timedelta(1,'s')
#             excepted_maxdiff = 120 # 
#             assert(tdiffmax < excepted_maxdiff), f'There is a jump in the data of more than {excepted_maxdiff} s ({tdiffmax} s). Can we except that?'
#         L1 = L1.sortby('time') # this is still needed if a negative diff was observed but allowed
#     elif time_shift_error == 'ignore':
#         pass
#     else:
#         assert(False), f'{time_shift_error} is not a valid value for time_shift_error'
#     return cl.CeilometerData(L1)

def read_netcdf_level1(file, 
            parent = None, 
            ignore = [],
            # name','message_type','version',
            # 'date_stamp',
            # 'period','tilt_angle',
            # 'cloud_status','cloud_data','status_bits','profile_scale',
            # 'profile_resolution','profile_length'],
            time_shift_error = 'raise'):
    if isinstance(file, (str, _pl.Path)):
        file = [file]
        
    assert(isinstance(file, (_pd.Series,list, _np.array))), f'File type not recognized: {type(file)}'
    
    # ignore1 = ['name','message_type','version',
    #            'date_stamp',
    #            'period','tilt_angle',
    #           'cloud_status','cloud_data','status_bits','profile_scale',
    #           'profile_resolution','profile_length']
    
    if not _np.all([_magic.from_file(fn.as_posix()) == 'Hierarchical Data Format (version 5) data' for fn in file]):
        fnc = '\n\t'.join([fn.as_posix() for fn in file])
        raise CorruptFileError(f'At least one of the following can not be identified as a netcdf file: \n\t {fnc}')
    # parent.tp_file = file.copy()    
    # cause problems:L1 = _xr.open_mfdataset(file, drop_variables=ignore1) # , concat_dim = 'timeDim' the usage of this kwarg seems to have changed and is not needed anymore
    L1 = _xr.open_mfdataset(file, drop_variables=ignore, concat_dim = 'timeDim', combine='nested')
    
    #### make time stamps
    L1 = L1.assign_coords(time = _pd.to_datetime(L1.time.values, unit = 's'))
    # parent.tp_ = 
    for var in L1.variables:
        if 'timeDim' in L1[var].dims:
            L1[var] = L1[var].swap_dims({'timeDim':'time'})
        
    if time_shift_error == 'raise':
        tdiff = L1.time.values[1:] - L1.time.values[:-1]
        tdiffmin = tdiff.min()/_pd.to_timedelta(1,'s')
        
        excepted_mindiff = -60
        assert(tdiffmin > excepted_mindiff), f'There is a jump in the time of {tdiffmin} s, which is lareger than the excepted value of {excepted_mindiff}. Probably due to time update? Programming required'
        if 0:
            # the maxdiff is problematic, as there are cases where large chuncs of data is missing
            tdiffmax = tdiff.max()/_pd.to_timedelta(1,'s')
            excepted_maxdiff = 120 # 
            assert(tdiffmax < excepted_maxdiff), f'There is a jump in the data of more than {excepted_maxdiff} s ({tdiffmax} s). Can we except that?'
        L1 = L1.sortby('time') # this is still needed if a negative diff was observed but allowed
    elif time_shift_error == 'ignore':
        pass
    else:
        assert(False), f'{time_shift_error} is not a valid value for time_shift_error'
    ds = L1
    return cl.CeilometerData(ds)

def read_netcdf_level2(p2f, ):
    ds = _xr.open_dataset(p2f)

    #### make time index
    ds = ds.assign_coords(time = _pd.to_datetime(ds.time.values, unit = 's'))
    
    #### change the original time index to the new one for all variables that need them
    for var in ds.variables:
        if 'timeDim' in ds[var].dims:
            ds[var] = ds[var].swap_dims({'timeDim':'time'})
    
    ds['cloud_status'] = ds.cloud_status.isel(cloud_statusDim = 0)
    
    dst = ds.cloud_data.rename({'cloud_dataDim': 'cloud_layer'})
    dst = dst.assign_coords(cloud_layer = _np.int8(dst.cloud_layer + 1))
    dst = dst.where(dst > 0)
    ds['cloud_data'] = dst 
    
    # For some reason the coordinate variables are sometimes missing
    # try:
    ds.attrs['site_latitude'] = float(ds.latitude.values[0])
    ds.attrs['site_longitude'] = float(ds.longitude.values[0])
    # except AttributeError:
    #     ds.attrs['site_latitude'] = _np.nan
    #     ds.attrs['site_longitude'] = _np.nan
    
    ds = ds.drop_duplicates('time')
    return cl.CeilometerData(ds)
#### read hist files



def read_hist_level2(p2f):
    def read_data_line(rein, labels, n = 5, n_levels = 1540):
        line = rein.readline().strip()
        if line == '':
            return False
        data = dict(zip(labels, line.split(', ')))
        
        timestamp = _pd.to_datetime(data['CREATEDATE'])
        timestampfromunixtime = _pd.to_datetime(int(data['UNIXTIME']), unit = 's')
        assert(timestamp == timestampfromunixtime), f'Creation time not equal to Unixtime, {timestamp} vs {timestampfromunixtime}.'
        # decode the backscatter profile
        bsp = data['BS_PROFILE']
        # if len(bsp)/n_levels != n:
        #     return data
        global tp_data
        tp_data = dict(data = data,
                       bsp = bsp,
                       timestamp = timestamp, 
                       line = line)
        if len(bsp)/n_levels != n:
            raise IncompleteProfileError(f'The length of the hex-string that discribes the profile is too short, {len(bsp)/n_levels} instead of {n}')
        
        # split in equal sections of n
        psps = [bsp[i:i+n] for i in range(0, len(bsp), n)]
        
        # hex to int
        psps = [int(bs,16) for bs in psps]
        
        # bit shift operation, not 100% sure how to interpret this ... after the first 19 bits you do a roll over?
        psps = [bs if bs < (1<<(n*4 - 1)) else bs - (1<<(n*4)) for bs in psps ]
        
        # make profile datafram
        level_thickness = 10
        levels = range(level_thickness,(n_levels + 1) * level_thickness, level_thickness) 
        
        # convert to xarray dataarray
        df = _pd.DataFrame(psps, columns=['backscatter_profile'], index = levels)
        df.index.name = 'range'
        
        # add time stamp
        ds = _xr.Dataset()
        da = _xr.DataArray(df.backscatter_profile)
        ds['profile_data'] = da
        ds['unixtime'] = int(data['UNIXTIME'])
        ds['period'] = int(data['PERIOD'])
        ds['ceilometer'] = data['CEILOMETER']
        ds = ds.expand_dims({'time':[timestamp,]})
        return ds

    with open(p2f, 'r') as rein:
        #fist line
        line = rein.readline()#.decode()
        assert(line.strip() == 'History file'), f'first line indicates this is not a history file, it sais {line.strip()}.'
        
        # labels
        line = rein.readline()#.decode()
        labels = [l.strip() for l in line.split(', ')]
        
        i = 0
        while 1:
            try:
                dst = read_data_line(rein, labels)
            except IncompleteProfileError:
                _warnings.warn("The his file incountered an incomplete profile.", IncompleteProfileWarning)
                continue
            if not dst:
                # print('this is the end')
                break
            if i == 0:
                ds = dst
            else:
                ds = _xr.concat([ds, dst], 'time') 
            i += 1
    
    ds.attrs['site_latitude'] = _np.nan
    ds.attrs['site_longitude'] = _np.nan
    ds = ds.drop_duplicates('time')
    
    return cl.CeilometerData(ds)

def read_hist_level3(file, parent = None):
    def read_file(fn):
        cols = ['CREATEDATE',' CEILOMETER',' CLOUD_STATUS',' CLOUD_1',' CLOUD_2',
                ' CLOUD_3'] # What columns to keep.
        his3 = _pd.read_csv(fn, skiprows=1, header=0, sep=',',
                           na_values='-9999', index_col=0, parse_dates=True,
                           # infer_datetime_format=True, #deprecated in pandas, If format is not given format will be infered by default now
                           usecols=cols)
        his3.index.rename('time', inplace=True)  
        his3.columns = [col.strip() for col in his3.columns]
        return his3
    
    if isinstance(file, (str, _pl.Path)):
        file = [file]
        
    assert(isinstance(file, (_pd.Series,list, _np.array))), f'File type not recognized: {type(file)}'
    df = _pd.concat([read_file(fn) for fn in file], sort = True)
    #### testpoint
    #parent.tp_dfcc = df.copy()
    # assert(df.index.duplicated().sum() == 0), 'There are duplicates in the hist file ... I would think this should not be happening. if it does un-comment the following line'
    
    #### make unified xarray dataset
    df = df[~df.index.duplicated(keep='first')] # Remove duplicates
    
    ds = _xr.Dataset()
    ds['cloud_status'] = df.CLOUD_STATUS
    dfhcl = df.loc[:,['CLOUD_1', 'CLOUD_2', 'CLOUD_3']]
    dfhcl.columns = [_np.int8(c[-1]) for c in dfhcl.columns]
    dfhcl.columns.name = 'cloud_layer'
    
    dfhcl[dfhcl <= 0] = _np.nan
    
    ds['cloud_data'] = dfhcl
    
    ds = ds.drop_duplicates('time')
    return cl.CeilometerData(ds)

