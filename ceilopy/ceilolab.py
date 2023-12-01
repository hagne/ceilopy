import xarray as _xr
import pandas as _pd
import numpy as _np
import pathlib as _pl
import traceback as _tb
import datetime as _dt
from email.mime.text import MIMEText as _MIMEText
import smtplib as _smtplib
import pathlib as __pl
import configparser as _cp

import matplotlib.colors as _mplcolors
import matplotlib.pyplot as _plt
import matplotlib.dates as _mdates
# import magic as _magic
# import warnings as _warnings


settings = """
[notify]
email_address = None
smtp = localhost

"""
def generate_config(p2sf):
    if not p2sf.parent.is_dir():
        p2sf.parent.mkdir()
    with open(p2sf, 'w') as raus:
        raus.write(settings)

def load_config():
    p2sf = __pl.Path.home().joinpath('.ceilopy/config.ini')

    if not p2sf.is_file():
        generate_config(p2sf)

    config = _cp.ConfigParser()
    config.read(p2sf)
    return config

# class CorruptFileError(Exception):
#     """Exception raised when File is not whats expected.
#     """
#     def __init__(self, message):
#         super().__init__(message)
        
# class MissingSerialNumberError(Exception):
#     """Exception raised when File does not contain Serial number.
#     """
#     def __init__(self, message):
#         super().__init__(message)
        
# class SerialNumberMissmatchError(Exception):
#     """Exception raised when Files doe not have the same serial number.
#     """
#     def __init__(self, message):
#         super().__init__(message)



# def read_L1(file, 
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
#     return L1

# def read_L1v2(file, 
#             parent = None, 
#             ignore = [],
#             # name','message_type','version',
#             # 'date_stamp',
#             # 'period','tilt_angle',
#             # 'cloud_status','cloud_data','status_bits','profile_scale',
#             # 'profile_resolution','profile_length'],
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
#     L1 = _xr.open_mfdataset(file, drop_variables=ignore, concat_dim = 'timeDim', combine='nested')
    
#     #### make time stamps
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
#     return L1

# def read_netcdf_level2_file(p2f, ):
#     ds = _xr.open_dataset(p2f)

#     #### make time index
#     ds = ds.assign_coords(time = _pd.to_datetime(ds.time.values, unit = 's'))
    
#     #### change the original time index to the new one for all variables that need them
#     for var in ds.variables:
#         if 'timeDim' in ds[var].dims:
#             ds[var] = ds[var].swap_dims({'timeDim':'time'})
    
#     ds['cloud_status'] = ds.cloud_status.isel(cloud_statusDim = 0)
    
#     dst = ds.cloud_data.rename({'cloud_dataDim': 'cloud_layer'})
#     dst = dst.assign_coords(cloud_layer = _np.int8(dst.cloud_layer + 1))
#     dst = dst.where(dst > 0)
#     ds['cloud_data'] = dst 
    
#     # For some reason the coordinate variables are sometimes missing
#     # try:
#     ds.attrs['site_latitude'] = float(ds.latitude.values[0])
#     ds.attrs['site_longitude'] = float(ds.longitude.values[0])
#     # except AttributeError:
#     #     ds.attrs['site_latitude'] = _np.nan
#     #     ds.attrs['site_longitude'] = _np.nan
    
#     ds = ds.drop_duplicates('time')
#     return ds
# #### read hist files

# class IncompleteProfileError(Exception):
#     pass

# class IncompleteProfileWarning(Warning):
#     pass

# def read_hist_level2_file(p2f):
#     def read_data_line(rein, labels, n = 5, n_levels = 1540):
#         line = rein.readline().strip()
#         if line == '':
#             return False
#         data = dict(zip(labels, line.split(', ')))
        
#         timestamp = _pd.to_datetime(data['CREATEDATE'])
#         timestampfromunixtime = _pd.to_datetime(int(data['UNIXTIME']), unit = 's')
#         assert(timestamp == timestampfromunixtime), f'Creation time not equal to Unixtime, {timestamp} vs {timestampfromunixtime}.'
#         # decode the backscatter profile
#         bsp = data['BS_PROFILE']
#         # if len(bsp)/n_levels != n:
#         #     return data
#         global tp_data
#         tp_data = dict(data = data,
#                        bsp = bsp,
#                        timestamp = timestamp, 
#                        line = line)
#         if len(bsp)/n_levels != n:
#             raise IncompleteProfileError(f'The length of the hex-string that discribes the profile is too short, {len(bsp)/n_levels} instead of {n}')
        
#         # split in equal sections of n
#         psps = [bsp[i:i+n] for i in range(0, len(bsp), n)]
        
#         # hex to int
#         psps = [int(bs,16) for bs in psps]
        
#         # bit shift operation, not 100% sure how to interpret this ... after the first 19 bits you do a roll over?
#         psps = [bs if bs < (1<<(n*4 - 1)) else bs - (1<<(n*4)) for bs in psps ]
        
#         # make profile datafram
#         level_thickness = 10
#         levels = range(level_thickness,(n_levels + 1) * level_thickness, level_thickness) 
        
#         # convert to xarray dataarray
#         df = _pd.DataFrame(psps, columns=['backscatter_profile'], index = levels)
#         df.index.name = 'range'
        
#         # add time stamp
#         ds = _xr.Dataset()
#         da = _xr.DataArray(df.backscatter_profile)
#         ds['profile_data'] = da
#         ds['unixtime'] = int(data['UNIXTIME'])
#         ds['period'] = int(data['PERIOD'])
#         ds['ceilometer'] = data['CEILOMETER']
#         ds = ds.expand_dims({'time':[timestamp,]})
#         return ds

#     with open(p2f, 'r') as rein:
#         #fist line
#         line = rein.readline()#.decode()
#         assert(line.strip() == 'History file'), f'first line indicates this is not a history file, it sais {line.strip()}.'
        
#         # labels
#         line = rein.readline()#.decode()
#         labels = [l.strip() for l in line.split(', ')]
        
#         i = 0
#         while 1:
#             try:
#                 dst = read_data_line(rein, labels)
#             except IncompleteProfileError:
#                 _warnings.warn("The his file incountered an incomplete profile.", IncompleteProfileWarning)
#                 continue
#             if not dst:
#                 # print('this is the end')
#                 break
#             if i == 0:
#                 ds = dst
#             else:
#                 ds = _xr.concat([ds, dst], 'time') 
#             i += 1
    
#     ds.attrs['site_latitude'] = _np.nan
#     ds.attrs['site_longitude'] = _np.nan
#     ds = ds.drop_duplicates('time')
    
#     return ds

# ##### Read Level3 hist files. #############################################
# def read_level3_hist(file, parent = None):
#     def read_file(fn):
#         cols = ['CREATEDATE',' CEILOMETER',' CLOUD_STATUS',' CLOUD_1',' CLOUD_2',
#                 ' CLOUD_3'] # What columns to keep.
#         his3 = _pd.read_csv(fn, skiprows=1, header=0, sep=',',
#                            na_values='-9999', index_col=0, parse_dates=True,
#                            # infer_datetime_format=True, #deprecated in pandas, If format is not given format will be infered by default now
#                            usecols=cols)
#         his3.index.rename('time', inplace=True)  
#         his3.columns = [col.strip() for col in his3.columns]
#         return his3
    
#     if isinstance(file, (str, _pl.Path)):
#         file = [file]
        
#     assert(isinstance(file, (_pd.Series,list, _np.array))), f'File type not recognized: {type(file)}'
#     df = _pd.concat([read_file(fn) for fn in file], sort = True)
#     #### testpoint
#     #parent.tp_dfcc = df.copy()
#     # assert(df.index.duplicated().sum() == 0), 'There are duplicates in the hist file ... I would think this should not be happening. if it does un-comment the following line'
    
#     #### make unified xarray dataset
#     df = df[~df.index.duplicated(keep='first')] # Remove duplicates
    
#     ds = _xr.Dataset()
#     ds['cloud_status'] = df.CLOUD_STATUS
#     dfhcl = df.loc[:,['CLOUD_1', 'CLOUD_2', 'CLOUD_3']]
#     dfhcl.columns = [_np.int8(c[-1]) for c in dfhcl.columns]
#     dfhcl.columns.name = 'cloud_layer'
    
#     dfhcl[dfhcl <= 0] = _np.nan
    
#     ds['cloud_data'] = dfhcl
    
#     ds = ds.drop_duplicates('time')
#     return ds

class CeilometerData(object):
    def __init__(self, dataset, site_name = None):
        self.dataset = dataset
        self.site_name = site_name
        
    def decorate_dataset(self,
                          version = 0,
                          site_name = 'None',
                          serial_no = 'None',
                          input_files = [],
                         ):
        """
        Decorates a dataset with units and other attributes.
    
        Parameters
        ----------
        ds : xr.Dataset
            Any dataset that comes out of the read function above.
        version : TYPE, optional
            DESCRIPTION. The default is 0.
        serial_no : TYPE, optional
            DESCRIPTION. The default is 'None'.
        input_files : TYPE, optional
            DESCRIPTION. The default is [].
        reindex : TYPE, optional
            DESCRIPTION. The default is False.
         : TYPE
            DESCRIPTION.
    
        Returns
        -------
        dsr : TYPE
            DESCRIPTION.
    
        """
        ds = self.dataset
        dsr = _xr.Dataset()
        dsr['backscatter_profile']=ds.profile_data.astype(_np.float32)
        
        dsr.backscatter_profile.attrs = {'long_name':'2-D ceilometer signal backscatter profile.',
                                         'units':'10e-9 m^-1 sr^-1',
                                         'comments':'Range-corrected-scattering'}
        
        # dsr['cloud_status'] = ds.cloud_status.isel(cloud_statusDim = 0).astype(_np.float32)
        dsr['cloud_status'] = ds.cloud_status.astype(_np.float32)
        
        dsr.cloud_status.attrs = {'long_name':'Cloud detection status.',
                                          'units':'1',
                                          'flag_values':'0,1,2,3,4',
                                          'flag_0':'No significant backscatter.',
                                          'flag_1':'One cloud layer detected.',
                                          'flag_2':'Two cloud layers detected.',
                                          'flag_3':'Three cloud layers detected.',
                                          'flag_4':'''Full obscuration/vertical visibility mode. First cloud_base will report vertical visibility and second cloud_base will report highest signal''',
                                          'comments':'''When cloud_status=4 there is an optically thick cloud that obscures the signal. Therefore it is not possible to discern additional cloud layers above it so the vertical visibility and highest signal are reported instead.'''}
        
        #### cloud data
        # dst = ds.cloud_data.rename({'cloud_dataDim': 'cloud_layer'})
        # dst = dst.assign_coords(cloud_layer = _np.int8(dst.cloud_layer + 1))
        # dst = dst.where(dst > 0) #for The L1 it looks like 0 is invalid
        dsr['cloud_data'] = ds.cloud_data.astype(_np.float32)
        
        dsr.cloud_data.attrs = {'long_name':'Cloud base heights.',
                                        'units':'m',
                                        'comments':'''A 2D array containing all three cloud bases at each timestep. -999 if no significant signal.''',
                                        'cloud_base_1':'''First cloud base height or vertical visibility if cloud_status=4''',
                                        'cloud_base_2':'''Second cloud base height or highest received signal if cloud_status=4''',
                                        'cloud_base_3':'Third cloud base height'}
        
        
        #### attributes of the cloud layer
        dsr.cloud_layer.attrs = dict(long_name = "Cloud layer index",
                                    units = "1" ,
                                    axis = "Z" ,
                                    positive = "up",)
        
        dsr.range.attrs = {'long_name': 'Distance from ground',
                           'units': 'm',
                           'standard_name' : 'distance',
                           'axis' : 'Z',
                           'positive': 'up'}
        
        #### global attributes
        dsr.attrs = {'title':'Ceilometer cloud product',
                     'version':version,
                     'institution':'NOAA/GML/GRAD',
                     'author':'hagen.telg@noaa.gov',
                     'source':'Vaisala CL51 ceilometer',
                     'serial_number': serial_no,
                     'input_files': ', '.join(input_files),#[fn.name for fn in poutg.path2raw]),
                     'Conventions':'CF-1.8',
                     'comments': ('The "time" coordinate was re-indexed to the full minute by back-filling the nearest '
                                  'valid data value within the following minute.The data values have not undergone any processing other than what the Vaisala '
                                  'software has applied. In addition, no QC has been applied other than a visual inspection for impossible values or obvious errors.')
                                                             }
        
        dsr.attrs['site_name'] = site_name
        dsr.attrs['site_latitude'] = f"{ds.attrs['site_latitude']:0.6f}"
        dsr.attrs['site_longitude'] = f"{ds.attrs['site_longitude']:0.6f}"
                    
        
        return CeilometerData(dsr)

    def remove_timeinconsitancies(self):
        """
        Fixes negative jumps in time, sorry McFly. This is only a quick-fix that 
        does not actually remove an apperent inconsistancy in the timestamps.
        The function takes the data from after the jump replaces all the data that 
        is older but ahead of it.
        Multiple jumps can be handled.
        
        Parameters
        ----------
        dsr : TYPE
            DESCRIPTION.
    
        Returns
        -------
        dsr : TYPE
            DESCRIPTION.
        log : TYPE
            DESCRIPTION.
    
        """
        dsr = self.dataset
        excepted_mindiff = -60 # such smoll differences will simply be fixed by ordering the index
        df = _pd.DataFrame({'time': dsr.time[:-1]})
        
        log = []
        #I am looping in case there are muliple time jumps
        while (mindiff := ((dsr.time.values[1:] - dsr.time.values[:-1]).min())/_pd.to_timedelta(1,'s')) < excepted_mindiff:
            logt = {'mindiff': mindiff}
            
            df['tdiff'] = dsr.time.values[1:] - dsr.time.values[:-1]
            df['tdiffs'] = df.tdiff/_pd.to_timedelta(1,'s')
            
            #index of last ocurring negative time jump
            idx = (df[::-1].tdiffs < excepted_mindiff).idxmax()
            idx_end = idx + 1
            tsaj = df.loc[idx_end].time #first timestamp after the jump back
            
            # find the last timestamp and its index that is smaller than this one
            idxstart = ((df.time - tsaj) < _pd.to_timedelta(0, 's'))[::-1].idxmax()
            idxval =  list(range(dsr.time.shape[0]))
    
            # select from datset
            dsr = dsr.isel(time = idxval[:idxstart] + idxval[idx_end:])
            log.append(logt)
        # if 0:
        dsr = dsr.sortby("time") # in case there is a small inconsitancy
        dsr = dsr.drop_duplicates('time') # just in case
        cdi = CeilometerData(dsr)
        cdi.log = log
        return cdi

    def reindex(self):
        """
        Reindex on 1min timeintervals and fill any gaps with nan ... every minute 
        of the day has a timestamp! Reindexing is backfilling and limited to 1 
        minute, so the closest data after the timestamp up to one minute is used as
        the value at the new index.
    
        Parameters
        ----------
        dsr : TYPE
            DESCRIPTION.
    
        Returns
        -------
        dsr : TYPE
            DESCRIPTION.
    
        """
        dsr = self.dataset
        time = dsr.time.to_pandas()
        start = time.iloc[0].date()
        end = start + _pd.to_timedelta(1, 'd')
        newtime = _pd.date_range(start = start, end = end, freq = _pd.to_timedelta(1, 'minute'), inclusive = 'left')
        
        assert(time.iloc[0].date() == time.iloc[-1].date()), f'Fist and last timestamp are not on the same date!! {time.iloc[0].date()} vs {time.iloc[-1].date()}'
        dsr['instrumen_reported_time'] = dsr.time.copy() #this ensures we still have the original timestamp for each value
        dsr.instrumen_reported_time.attrs = {'long_name': "Instrument reported time",
                                             'comments': 'Timestamp that was reported by the instrument. The "time" coordinate was re-indexed to the full minute (see global comments).'}
        dsr = dsr.reindex(time = newtime, method = 'bfill', tolerance = '1min')
        return CeilometerData(dsr)
    
    def plot_quicklooks(self):
        ds = self.dataset
        site = ds.attrs['site_name']
        f,aa = _plt.subplots(2,2, 
                            # sharex=True, 
                            width_ratios=[30,1],  gridspec_kw={'hspace': 0, 'wspace': 0.05})
        scale = 1.2
        f.set_figwidth(f.get_figwidth() * scale * 1.4)
        f.set_figheight(f.get_figheight() * scale)
        aa = aa.flatten()
        if 1:
            a = aa[0]
            pc = ds.backscatter_profile.plot.pcolormesh(x = 'time', ax = a, add_colorbar = False)
            pc.set_norm(_mplcolors.LogNorm(vmin = 5e1, vmax = 2e5))
            pc.set_cmap(_plt.cm.inferno_r)
            
            dt = _pd.to_datetime(ds.time.values[0])
            # a.set_title(f'{dt.year}-{dt.month:02d}-{dt.day:02d}')
            
            a.set_title(f'{dt.year}-{dt.month:02d}-{dt.day:02d} - {site} - L51 ceilometer backscatter (top), cloud-layers (bottom)')
            # a.text(0.1,0.9, f'Site: {site}',transform = a.transAxes)
            a.set_yscale('log')
            a.set_ylim(bottom = 2e2)
            a.xaxis.set_ticks([])
        
        a = aa[2]
        out = ds.cloud_data.plot.line(x = 'time', ls = '', marker = '_', ax = a)
        leg = a.get_legend()
        leg.set_title('layer idx.')
        a.set_xlim(aa[0].get_xlim())
        a.xaxis.set_major_locator(_mdates.HourLocator(interval = 2))
        a.xaxis.set_minor_locator(_mdates.HourLocator(interval = 1))
        a.xaxis.set_major_formatter(_mdates.DateFormatter('%H'))
        a.tick_params(axis = 'x', 
                      reset = True,
                      # labelrotation = 0, 
                     )
        a.set_xlabel('Houre of day [UTC]')
        
        # for a in aa[,3]]:
        aa[3].axis('off')
        
        cbar = f.colorbar(pc, cax = aa[1],
                   # location = 'center', 
                   fraction = 1)
        cbar.set_label('Range-corrected backscattering')
        
        return f,aa,pc

#### Dataset operations
# def decorate_dataset( ds,
#                       version = 0,
#                       site_name = 'None',
#                       serial_no = 'None',
#                       input_files = [],
#                      ):
#     """
#     Decorates a dataset with units and other attributes.

#     Parameters
#     ----------
#     ds : xr.Dataset
#         Any dataset that comes out of the read function above.
#     version : TYPE, optional
#         DESCRIPTION. The default is 0.
#     serial_no : TYPE, optional
#         DESCRIPTION. The default is 'None'.
#     input_files : TYPE, optional
#         DESCRIPTION. The default is [].
#     reindex : TYPE, optional
#         DESCRIPTION. The default is False.
#      : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     dsr : TYPE
#         DESCRIPTION.

#     """
#     dsr = _xr.Dataset()
#     dsr['backscatter_profile']=ds.profile_data.astype(_np.float32)
    
#     dsr.backscatter_profile.attrs = {'long_name':'2-D ceilometer signal backscatter profile.',
#                                      'units':'10e-9 m^-1 sr^-1',
#                                      'comments':'Range-corrected-scattering'}
    
#     # dsr['cloud_status'] = ds.cloud_status.isel(cloud_statusDim = 0).astype(_np.float32)
#     dsr['cloud_status'] = ds.cloud_status.astype(_np.float32)
    
#     dsr.cloud_status.attrs = {'long_name':'Cloud detection status.',
#                                       'units':'1',
#                                       'flag_values':'0,1,2,3,4',
#                                       'flag_0':'No significant backscatter.',
#                                       'flag_1':'One cloud layer detected.',
#                                       'flag_2':'Two cloud layers detected.',
#                                       'flag_3':'Three cloud layers detected.',
#                                       'flag_4':'''Full obscuration/vertical visibility mode. First cloud_base will report vertical visibility and second cloud_base will report highest signal''',
#                                       'comments':'''When cloud_status=4 there is an optically thick cloud that obscures the signal. Therefore it is not possible to discern additional cloud layers above it so the vertical visibility and highest signal are reported instead.'''}
    
#     #### cloud data
#     # dst = ds.cloud_data.rename({'cloud_dataDim': 'cloud_layer'})
#     # dst = dst.assign_coords(cloud_layer = _np.int8(dst.cloud_layer + 1))
#     # dst = dst.where(dst > 0) #for The L1 it looks like 0 is invalid
#     dsr['cloud_data'] = ds.cloud_data.astype(_np.float32)
    
#     dsr.cloud_data.attrs = {'long_name':'Cloud base heights.',
#                                     'units':'m',
#                                     'comments':'''A 2D array containing all three cloud bases at each timestep. -999 if no significant signal.''',
#                                     'cloud_base_1':'''First cloud base height or vertical visibility if cloud_status=4''',
#                                     'cloud_base_2':'''Second cloud base height or highest received signal if cloud_status=4''',
#                                     'cloud_base_3':'Third cloud base height'}
    
    
#     #### attributes of the cloud layer
#     dsr.cloud_layer.attrs = dict(long_name = "Cloud layer index",
#                                 units = "1" ,
#                                 axis = "Z" ,
#                                 positive = "up",)
    
#     dsr.range.attrs = {'long_name': 'distance from ground',
#                        'units': 'm',
#                        'standard_name' : 'distance',
#                        'axis' : 'Z',
#                        'positive': 'up'}
    
#     #### global attributes
#     dsr.attrs = {'title':'Ceilometer cloud product',
#                  'version':version,
#                  'institution':'NOAA/GML/GRAD',
#                  'author':'hagen.telg@noaa.gov',
#                  'source':'Vaisala CL51 ceilometer',
#                  'serial_number': serial_no,
#                  'input_files': ', '.join(input_files),#[fn.name for fn in poutg.path2raw]),
#                  'Conventions':'CF-1.8',
#                  'comments': ('The "time" coordinate was re-indexed to the full minute by back-filling the nearest '
#                               'valid data value within the following minute.The data values have not undergone any processing other than what the Vaisala '
#                               'software has applied. In addition, no QC has been applied other than a visual inspection for impossible values or obvious errors.')
#                                                          }
    
#     dsr.attrs['site_name'] = site_name
#     dsr.attrs['site_latitude'] = f"{ds.attrs['site_latitude']:0.6f}"
#     dsr.attrs['site_longitude'] = f"{ds.attrs['site_longitude']:0.6f}"
                
#     return dsr

# def remove_timeinconsitancies(dsr):
#     """
#     Fixes negative jumps in time, sorry McFly. This is only a quick-fix that 
#     does not actually remove an apperent inconsistancy in the timestamps.
#     The function takes the data from after the jump replaces all the data that 
#     is older but ahead of it.
#     Multiple jumps can be handled.
    
#     Parameters
#     ----------
#     dsr : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     dsr : TYPE
#         DESCRIPTION.
#     log : TYPE
#         DESCRIPTION.

#     """
#     excepted_mindiff = -60 # such smoll differences will simply be fixed by ordering the index
#     df = _pd.DataFrame({'time': dsr.time[:-1]})
    
#     log = []
#     #I am looping in case there are muliple time jumps
#     while (mindiff := ((dsr.time.values[1:] - dsr.time.values[:-1]).min())/_pd.to_timedelta(1,'s')) < excepted_mindiff:
#         logt = {'mindiff': mindiff}
        
#         df['tdiff'] = dsr.time.values[1:] - dsr.time.values[:-1]
#         df['tdiffs'] = df.tdiff/_pd.to_timedelta(1,'s')
        
#         #index of last ocurring negative time jump
#         idx = (df[::-1].tdiffs < excepted_mindiff).idxmax()
#         idx_end = idx + 1
#         tsaj = df.loc[idx_end].time #first timestamp after the jump back
        
#         # find the last timestamp and its index that is smaller than this one
#         idxstart = ((df.time - tsaj) < _pd.to_timedelta(0, 's'))[::-1].idxmax()
#         idxval =  list(range(dsr.time.shape[0]))

#         # select from datset
#         dsr = dsr.isel(time = idxval[:idxstart] + idxval[idx_end:])
#         log.append(logt)
#     # if 0:
#     dsr = dsr.sortby("time") # in case there is a small inconsitancy
#     dsr = dsr.drop_duplicates('time') # just in case
#     return dsr, log

# def reindex(dsr):
#     """
#     Reindex on 1min timeintervals and fill any gaps with nan ... every minute 
#     of the day has a timestamp! Reindexing is backfilling and limited to 1 
#     minute, so the closest data after the timestamp up to one minute is used as
#     the value at the new index.

#     Parameters
#     ----------
#     dsr : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     dsr : TYPE
#         DESCRIPTION.

#     """
#     time = dsr.time.to_pandas()
#     start = time.iloc[0].date()
#     end = start + _pd.to_timedelta(1, 'd')
#     newtime = _pd.date_range(start = start, end = end, freq = _pd.to_timedelta(1, 'minute'), inclusive = 'left')
    
#     assert(time.iloc[0].date() == time.iloc[-1].date()), f'Fist and last timestamp are not on the same date!! {time.iloc[0].date()} vs {time.iloc[-1].date()}'
#     dsr['instrumen_reported_time'] = dsr.time.copy() #this ensures we still have the original timestamp for each value
#     dsr.instrumen_reported_time.attrs = {'long_name': "Instrument reported time",
#                                          'comments': 'Timestamp that was reported by the instrument. The "time" coordinate was re-indexed to the full minute (see global comments).'}
#     dsr = dsr.reindex(time = newtime, method = 'bfill', tolerance = '1min')
#     return dsr
#######
#### products

# class Cl51CloudProdRetriever_v1p3():
#     def __init__(self, poutg, version,
#                  # check_serial = True,
#                  ):
#         """
#         There was a time when the product was a aggregate of different functions, hist and netcdf. 
#         In todays version this might not be needed anymore... Consider consoidating this into the processor class unless?!?

#         Parameters
#         ----------
#         poutg : TYPE
#             DESCRIPTION.
#         version : TYPE
#             DESCRIPTION.
#         # check_serial : TYPE, optional
#             DESCRIPTION. The default is True.
#          : TYPE
#             DESCRIPTION.

#         Returns
#         -------
#         None.

#         """
#         self.version = version
#         self.poutg = poutg
# #         self.p2fnout = poutg.path2fn_out.unique()[0]
#         self._product_dataset = None
#         # self.get_serial_numbers()
#         # if check_serial:
#         #     self.check_serial()
        
#     # def get_serial_numbers(self):
#     #     def get_serial(row):
#     #         # Extract serial numbers from files
#     #         key = row.file_type
#     #         file = row.path2raw
#     #         if key in ['L1', 'L2', 'bl']:
#     #             serial = file.name[-11:-3]  # Extract serial number from L1 filename.
#     #         # elif key == 'L3':
#     #         #     serial = files['L3'][-11:-3]
#     #         elif key in ['H2','H3','hist']:
#     #             h = _pd.read_csv(file, skiprows=1, header=0, sep=',')
#     #             serial = h[' CEILOMETER'][0].strip() # Extract serial number from H2 file.
#     #         else:
#     #             raise KeyError('File type unknown')
#     #         return serial
        
#     #     self.poutg['sn'] = self.poutg.apply(get_serial, axis = 1)
        
#     # def check_serial(self, error_handling = 'raise'):
#     #     """
#     #     Currently not used anyChecks if the serial numbers in all the files are the same. In early 
#     #     measurments the serial number was not stored ... use error_handling to
#     #     deal with occuring errors.

#     #     Parameters
#     #     ----------
#     #     error_handling : str, optional
#     #         How to deal with errors. The default is 'raise'.
#     #         raise: raises occuring errors
#     #         allow_empty: do not raise an error if serial number is not available

#     #     Raises
#     #     ------
#     #     KeyError
#     #         DESCRIPTION.

#     #     Returns
#     #     -------
#     #     serial : TYPE
#     #         DESCRIPTION.

#     #     """
#     #     sn_series = self.poutg['sn'].copy()
#     #     # self.poutg['sn'] = sn_series.copy()
#     #     valid = ['raise', 'allow_empty']
#     #     assert(error_handling in valid), f'error_handling got an unexpected value ({error_handling}. Choose from: {valid})'
#     #     if error_handling == 'allow_empty':
#     #         sn_series = sn_series[sn_series.apply(lambda x: len(x)) != 0]
#     #     if sn_series.unique().shape[0] != 1:
#     #         if len(sn_series[sn_series.apply(lambda x: len(x)) != 0]) != len(sn_series):
#     #             fnj = '\n\t'.join([fn.as_posix() for fn in self.poutg.path2raw])
#     #             raise MissingSerialNumberError(f'At least one of the following files is missing a serial number:\n\t{fnj}')
#     #         raise SerialNumberMissmatchError(f'Serial numbers ({sn_series.unique()}) do not match')
        
#     @property
#     def product_dataset(self):
#         if isinstance(self._product_dataset, type(None)):
#             poutg = self.poutg
            
#             try:
#                 dsl = []
#                 input_files = []
#                 for p in poutg.path2raw:
#                     dst = read_netcdf_level2_file(p)
#                     dsl.append(dst)
#                     input_files.append(p.name)
#                 ds = _xr.concat(dsl, 'time')
#             except OSError:
#                 dsl = []
#                 input_files = []
                
#                 #### load hist lev 2 file
#                 files = poutg.path2althist_l2.unique()
#                 assert(files.shape[0] ==1), 'I thought there is always just 1 hist file a day, even if there is a reboot ... unlike the netcdffiles'
#                 dsh2 = read_hist_level2_file(files[0])
#                 input_files.append(files[0].name)
                
#                 #### load hist lev 3 file
#                 files = poutg.path2althist_l3.unique()
#                 assert(files.shape[0] ==1), 'I thought there is always just 1 hist file a day, even if there is a reboot ... unlike the netcdffiles'
#                 dsh3 = read_level3_hist(files[0])
#                 input_files.append(files[0].name)    
                
#                 ds = _xr.merge([dsh2, dsh3])
#                 # raise

#             self.tp_dspd = ds.copy()
#             ds = decorate_dataset(ds, version=self.version, site_name = poutg.site[0], serial_no=poutg.serial_no[0], input_files=input_files)
#             ds,log = remove_timeinconsitancies(ds)
            
#             if len(log) > 0:
#                 print('--------')
#                 print(poutg.iloc[0].name)
#                 print(log)
                
#             ds = reindex(ds)
            
            
            
            
            
# #             ds = read_L1(poutg[poutg.file_type == 'bl'].path2raw, parent = self,
# #                          ignore1=['name', 'message_type', 'version', 'date_stamp', 'period', 'tilt_angle', 'status_bits', 'profile_scale', 'profile_resolution', 'profile_length'])

            
# #             dsr = _xr.Dataset()
# #             dsr['backscatter_profile']=ds.rcs_910.astype(_np.float32)
            
# #             dsr.backscatter_profile.attrs = {'long_name':'2-D ceilometer signal backscatter profile.',
# #                                              'units':'10e-9 m^-1 sr^-1',
# #                                              'comments':'Range-corrected-scattering'}
            
# #             dsr['cloud_status'] = ds.cloud_status.isel(cloud_statusDim = 0).astype(_np.float32)
            
# #             dsr.cloud_status.attrs = {'long_name':'Cloud detection status.',
# #                                               'units':'1',
# #                                               'flag_values':'0,1,2,3,4',
# #                                               'flag_0':'No significant backscatter.',
# #                                               'flag_1':'One cloud layer detected.',
# #                                               'flag_2':'Two cloud layers detected.',
# #                                               'flag_3':'Three cloud layers detected.',
# #                                               'flag_4':'''Full obscuration/vertical visibility mode. First cloud_base will report vertical visibility and second cloud_base will report highest signal''',
# #                                               'comments':'''When cloud_status=4 there is an optically thick cloud that obscures the signal. Therefore it is not possible to discern additional cloud layers above it so the vertical visibility and highest signal are reported instead.'''}
            
# #             dst = ds.cloud_data.rename({'cloud_dataDim': 'cloud_layer'})
# #             dst = dst.assign_coords(cloud_layer = _np.int8(dst.cloud_layer + 1))
# #             dst = dst.where(dst > 0) #for The L1 it looks like 0 is invalid
# #             dsr['cloud_data'] = dst.astype(_np.float32)
            
# #             dsr.cloud_data.attrs = {'long_name':'Cloud base heights.',
# #                                             'units':'m',
# #                                             'comments':'''A 2D array containing all three cloud bases at each timestep. -999 if no significant signal.''',
# #                                             'cloud_base_1':'''First cloud base height or vertical visibility if cloud_status=4''',
# #                                             'cloud_base_2':'''Second cloud base height or highest received signal if cloud_status=4''',
# #                                             'cloud_base_3':'Third cloud base height'}
            
# #             dsr.attrs = {'title':'Ceilometer cloud product',
# #                          'version':self.version,
# #                          'institution':'NOAA/GML/GRAD',
# #                          'author':'hagen.telg@noaa.gov',
# #                          'source':'Vaisala CL51 ceilometer',
# #                          'serial_number': poutg.sn.unique()[0],
# #                          'input_files': ', '.join([fn.name for fn in poutg.path2raw]),
# #                          'Conventions':'CF-1.8',
# #                          'comments':''' The "time" coordinate was re-indexed to the full minute by back-filling the nearest 
# # valid data value within the following minute.The data values have not undergone any processing other than what the Vaisala 
# # software has applied. In addition, no QC has been applied other than a visual inspection for impossible values or obvious errors.'''.replace('\n','')
# #                          }
            
# #             #### attributes of the coordinates
# #             dsr.cloud_layer.attrs = dict(long_name = "Cloud layer index",
# #                                         units = "1" ,
# #                                         axis = "Z" ,
# #                                         positive = "up",)

# #             dsr.range.attrs = {'long_name': 'distance from ground',
# #                                'units': 'm',
# #                                'standard_name' : 'distance',
# #                                'axis' : 'Z',
# #                                'positive': 'up'}
# #             #### reindex on 1min data and fill any gaps with nan ... every minute of the day has a timestamp!
# #             time = dsr.time.to_pandas()
# #             start = time.iloc[0].date()
# #             end = start + _pd.to_timedelta(1, 'd')
# #             newtime = _pd.date_range(start = start, end = end, freq = _pd.to_timedelta(1, 'minute'), inclusive = 'left')
            
# #             assert(time.iloc[0].date() == time.iloc[-1].date()), f'Fist and last timestamp are not on the same date!! {time.iloc[0].date()} vs {time.iloc[-1].date()}'
# #             dsr['instrumen_reported_time'] = dsr.time.copy() #this ensures we still have the original timestamp for each value
# #             dsr.instrumen_reported_time.attrs = {'long_name': "Instrument reported time",
# #                                                  'comments': 'Timestamp that was reported by the instrument. The "time" coordinate was re-indexed to the full minute (see global comments).'}
# #             dsr = dsr.reindex(time = newtime, method = 'bfill', tolerance = '1min')
            
#             self._product_dataset = ds
                        
            
#         return self._product_dataset

# class Cl51CloudProdProcessor_v1p3(object):
#     def __init__(self, 
#                  p2fl_in = '/nfs/grad/Inst/Ceil/SURFRAD/',
#                  p2fl_out = '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1_{version}',
#                  # version = '0',
#                  hist_file_format = '*_CEILOMETER_1_LEVEL_3*.his',
#                  ignore = [],
#                  verbose = False,
#                  ):
#         """
#         This class aggreates all functions that are needed to create the
#         cl51_cloud_prod_lev1 (formally known as cl51_cloud_prod_lev0) product.
#         I

#         Parameters
#         ----------
#         p2fl_in : str, optional
#             Path to input folder. This folder is expected to have further subfolders, one for each instrument (typically the measurment site, e.g. TBL)
#         p2fl_out : str, optional
#             Output folder. The default is '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1'.
#         hist_file_format : TYPE, optional
#             Format (pattern) of the hist files. As this format sometimes changes it not hardcoded. The default is '*_CEILOMETER_1_LEVEL_3*.his'.
#         ignore : list, optional
#             When there is a folder in the p2fl_in folder that mathes an element of this list it will be ignored.
#         verbose : bool, optional
#             If True there will be more output ... mostly for develpment purposes. The default is False.


#         Returns
#         -------
#         Cl51CloudProdProcessor instance.

#         """
#         self.version = '1.3'
#         self.p2fl_out = _pl.Path(p2fl_out.format(version = self.version))
#         self.ignore = ignore
#         self.p2fl_in = _pl.Path(p2fl_in)
#         self.hist_file_format= hist_file_format
#         self.bl_file_format = 'L2*.nc'
#         self.fn_format_out = '{site}.cl51.cloud_prod.{date}.nc'
#         self.verbose = verbose
#         # self.test = test
        
#         self._workplan = None
    
    
#     @property
#     def workplan(self):
#         if isinstance(self._workplan, type(None)):
#             if self.verbose:
#                 print('creating workplan')
#                 print('=================')
#             def bl2date(row):
#                 if row.path2raw.name.split('_')[-1].split('.')[0].isnumeric():
#                     dt = row.path2raw.name.split('_')[-1].split('.')[0]
#                 else:
#                     dt = row.path2raw.name.split('_')[-2]
#                 return _pd.to_datetime(dt)

#             workplan = _pd.DataFrame()

#             for p2site in self.p2fl_in.glob('*'):
#                 if not p2site.is_dir():
#                     continue
                
#                 if self.verbose:
#                     print(f'\t folder: {p2site}')
#                 if p2site.name in self.ignore:
#                     if self.verbose:
#                         print('\t\t ignore!')
#                     continue
                
#                 ## bl files
#                 p2fl_bl = p2site.joinpath('bl')
#                 df = _pd.DataFrame(p2fl_bl.glob(self.bl_file_format), columns = ['path2raw'])
                
#                 #### make datetime index
#                 df.index = df.apply(bl2date, axis = 1)

#                 ## concat  

#                 df['site'] = p2site.name.lower()
#                 workplan = _pd.concat([workplan, df])

#             workplan.sort_index(inplace=True)
            
#             #### add date collumn
#             workplan['date'] = workplan.apply(lambda row: _pd.to_datetime(row.name.date()), axis = 1)
            
#             ## create path to output file

#             def create_path_out(row):
#                 fname = self.fn_format_out.format(site = row.site, date = row.name.date().__str__().replace('-',''))
#                 p2fn = self.p2fl_out.joinpath(row.site).joinpath(f'{row.name.year}').joinpath(fname)
#                 return p2fn

#             workplan['path2fn_out'] = workplan.apply(create_path_out, axis = 1)
            
#             #### add serial number
#             workplan['serial_no'] = workplan.apply(lambda row: row.path2raw.name.split('_')[-1].replace('.nc',''), axis = 1)
            
#             #### remove row when output file exists

#             workplan = workplan[~ workplan.apply(lambda row: row.path2fn_out.is_file(), axis = 1)]

#             #### remove last day
            
#             workplan = workplan[~(workplan.date == workplan.date.iloc[-1]).values]
            
#             #### add alternative hist files. This is needed if netcdf is corrupt
#             workplan['path2althist_l2'] = workplan.apply(lambda row: row.path2raw.parent.parent.joinpath('hist').joinpath(f'{row.name.year}{row.name.month:02d}_CEILOMETER_1_LEVEL_2_{row.name.day:02d}.his'), axis = 1)
#             workplan['path2althist_l3'] = workplan.apply(lambda row: row.path2raw.parent.parent.joinpath('hist').joinpath(f'{row.name.year}{row.name.month:02d}_CEILOMETER_1_LEVEL_3_DEFAULT_{row.name.day:02d}.his'), axis = 1)
            
#             #### test if bl and hist exist and remove if not
#             #is it enough that bl and hist exists or do we have to check that it covers the entire day? We could execute this before deleting the last day ... would that be saver?
            
#             #workplan = workplan[workplan['path2fn_out'].map(workplan.groupby('path2fn_out').apply(lambda group: group.file_type.eq('bl').any()&group.file_type.eq('hist').any()))]
                        
#             self._workplan = workplan
#             if self.verbose:
#                 print('========= workplan done ==========')
#         return self._workplan
    
#     @workplan.setter
#     def workplan(self, value):
#         self._workplan = value
        
#     def process(self, test=False, 
#                 path2fn_out = None,
#                 generate_missing_folders = False,
#                 error_handling = 'raise',
#                 error_handling_serial='raise',
#                 verbose = False, 
#                 complevel = 4):
#         """
#         Processes the workplan

#         Parameters
#         ----------
#         test : [bool, int], optional
#             Creats several test scenarios:
#                 1: returns firt product dataset, stops processing afterwards
#                 2: as 1 but saves dataset
#         path2fn_out: str, optional, default is None
#             For testing. Excecute only files for this particular output path.
#         error_handling: str, optional, default is 'return'
#             What to do if an error accures during processing of single 
#             retrieval.
#             'raise': will raise the error
#             'return': will not raise an errorj, but will return it. Check 
#                 'errors' key of return dict.
            

#         Returns
#         -------
#         ds : TYPE
#             DESCRIPTION.

#         """
#         if self.verbose:
#             print('processing')
#             print('==========')
#         def make_all_parent_flds(file): 
#             if not file.parent.is_dir():
#                 make_all_parent_flds(file.parent)
#                 file.parent.mkdir(mode = 0o775)
#                 file.parent.chmod(0o775) # for some reason the mode kwarg in the line above did not give the desired result?!?
                
#         no_processed = 0
#         errors = []
#         errror_grp = []
#         out = {}
#         out['start_time'] = _pd.Timestamp(_dt.datetime.now())
#         for p2fnout, poutg in self.workplan.groupby('path2fn_out'): 
#             if verbose:
#                 print(f'\t path2fn_out: {p2fnout} - {poutg}')
#             if not isinstance(path2fn_out, type(None)):
#                 if p2fnout != _pl.Path(path2fn_out):
#                     continue
                              
#             ds = None
#             try:
#                 retriever = Cl51CloudProdRetriever_v1p3(poutg, self.version)
#                 self.tp_retriever = retriever
#                 ####Fixme retriever.check_serial(error_handling=error_handling_serial)
#                 ds =retriever.product_dataset
#                 no_processed += 1
                
#             except Exception as err:
#                 if error_handling == 'raise':
#                     raise
                
#                 #### error handling
#                 elif error_handling == 'return':
#                     errors.append(err)
#                     errror_grp.append(poutg)
#                     # errors.append(sys.exc_info())
#                     # error, error_txt, trace = sys.exc_info()
#                     # tm = ['{}: {}'.format(error.__name__, error_txt.args[0])] + _tb.format_tb(trace)
                
#             if test == 1:
#                 break
            
#             #### Write to disk.        
#             if not isinstance(ds, type(None)):
#                 assert(not p2fnout.is_file()), f'File exists, this should not be possible! path: {p2fnout}'
            
#                 if generate_missing_folders:
#                     make_all_parent_flds(p2fnout)
                
#                 #### encodings
#                 encoding = {}
#                 for var in ds:
#                     encoding[var] = {'zlib': True,     # Enable compression
#                                      'complevel': complevel   # Moderate compression level
#                                     }
#                 encoding['time'] = {'units': "minutes since 1900-01-01 00:00:00 UTC",
#                                    'zlib': True,     # Enable compression
#                                      'complevel': complevel   # Moderate compression level
#                                     }
#                 encoding['instrumen_reported_time'] = {'units': "minutes since 1900-01-01 00:00:00 UTC",
#                                                       'zlib': True,     # Enable compression
#                                                          'complevel': complevel   # Moderate compression level
#                                                         }
                
#                 #### save to netcdf
#                 ds.to_netcdf(path=p2fnout, 
#                              # mode='w', 
#                              encoding = encoding,
#                              format='NETCDF4')
#                 p2fnout.chmod(0o664) # to give the group writing access.
#             else:
#                 if verbose:
#                     print(f'{p2fnout} was not saved as ds is None.')
                
#             if test == 2:
#                 break
        
#         out['errors'] = errors
#         out['errror_grp'] = errror_grp
#         out['last_dataset'] = ds
#         out['no_of_files_processed'] = no_processed
#         out['finish_time'] = _pd.Timestamp(_dt.datetime.now())
#         self._last_processing = out
#         return out
    
#     def get_single_day_from_worplan(self, index = -1, random = False):
#         gl = [g for idc,g in self.workplan.groupby('path2fn_out')]
        
#         if random:
#             index = int(_np.round(_np.random.random(1) * (len(gl)-1)))
    
#         oneday = gl[index]
#         return oneday
    
#     def notify(self, subject = 'cl51cloudprod - status: {status} (clean: {no_of_files_processed}; errors: {no_of_errors})'):
#         config = load_config()
#         assert(config.get('notify', 'email_address') != 'None'), 'No email has been specified, _please do so in ~/.ceilopy/config.ini'
        
#         out = self._last_processing
#         messages = ['run started {}'.format(out['start_time'])]
#         messages.append('run finshed {}'.format(out['finish_time']))
        
#         #### summary
#         no_of_errors = len(out['errors'])
#         no_of_files_processed = out['no_of_files_processed']
            
#         messages.append('\n'.join([f'length of workplan:\t\t{self.workplan.shape[0]}',f'successfully created files:\t{no_of_files_processed}', f'errors encountered:\t\t{no_of_errors}']))
        
#         #### errors
#         if no_of_errors != 0:
#             errs = out['errors']
#             err_types = _np.array([err.__class__.__name__ for err in errs])
            
#             typoferrors = []
#             for toe in _np.unique(err_types):
#                 typoferrors.append({'name': toe, 'times': (err_types == toe).sum(), 'first_err': errs[_np.where(err_types == toe)[0][0]]})
        
#             messages.append('\n'.join(['Errors by type:',] + ['\t{tn}:\t{tt}'.format(tn = toe['name'], tt = toe['times']) for toe in typoferrors]))
#             messages.append('\n=============================================\n'.join(['First traceback for each error type',]+[''.join(_tb.format_tb(toe['first_err'].__traceback__) + [toe['first_err'].__str__(),]) for toe in typoferrors]))
        
#         #### email body
#         message_txt = '\n=========================================================================\n'.join(messages)
#         msg = _MIMEText(message_txt)
        
#         #### subject
#         if no_of_errors ==0:
#             status = 'clean'
#         else:
#             status = 'errors'
#         # subject = f'cl51cloudprod - status: {status} (clean: {no_of_files_processed}; errors: {no_of_errors})'
#         subject = subject.format(status = status, no_of_files_processed = no_of_files_processed, no_of_errors = no_of_errors)
#         address  = config.get('notify', 'email_address')
#         smtp = config.get('notify', 'smtp')
        
#         msg['Subject'] = subject
#         msg['From'] = address
#         msg['To'] = address
        
#         # Send the message via our own SMTP server.
#         s = _smtplib.SMTP(smtp)
#         s.send_message(msg)
#         s.quit()



class Cl51CloudProdRetriever_v1p2():
    def __init__(self, poutg, version,
                 # check_serial = True,
                 ):
        self.version = version
        self.poutg = poutg
#         self.p2fnout = poutg.path2fn_out.unique()[0]
        self._product_dataset = None
        self.get_serial_numbers()
        # if check_serial:
        #     self.check_serial()
        
    def get_serial_numbers(self):
        def get_serial(row):
            # Extract serial numbers from files
            key = row.file_type
            file = row.path2raw
            if key in ['L1', 'L2', 'bl']:
                serial = file.name[-11:-3]  # Extract serial number from L1 filename.
            # elif key == 'L3':
            #     serial = files['L3'][-11:-3]
            elif key in ['H2','H3','hist']:
                h = _pd.read_csv(file, skiprows=1, header=0, sep=',')
                serial = h[' CEILOMETER'][0].strip() # Extract serial number from H2 file.
            else:
                raise KeyError('File type unknown')
            return serial
        
        self.poutg['sn'] = self.poutg.apply(get_serial, axis = 1)
        
    def check_serial(self, error_handling = 'raise'):
        """
        Checks if the serial numbers in all the files are the same. In early 
        measurments the serial number was not stored ... use error_handling to
        deal with occuring errors.

        Parameters
        ----------
        error_handling : str, optional
            How to deal with errors. The default is 'raise'.
            raise: raises occuring errors
            allow_empty: do not raise an error if serial number is not available

        Raises
        ------
        KeyError
            DESCRIPTION.

        Returns
        -------
        serial : TYPE
            DESCRIPTION.

        """
        sn_series = self.poutg['sn'].copy()
        # self.poutg['sn'] = sn_series.copy()
        valid = ['raise', 'allow_empty']
        assert(error_handling in valid), f'error_handling got an unexpected value ({error_handling}. Choose from: {valid})'
        if error_handling == 'allow_empty':
            sn_series = sn_series[sn_series.apply(lambda x: len(x)) != 0]
        if sn_series.unique().shape[0] != 1:
            if len(sn_series[sn_series.apply(lambda x: len(x)) != 0]) != len(sn_series):
                fnj = '\n\t'.join([fn.as_posix() for fn in self.poutg.path2raw])
                raise MissingSerialNumberError(f'At least one of the following files is missing a serial number:\n\t{fnj}')
            raise SerialNumberMissmatchError(f'Serial numbers ({sn_series.unique()}) do not match')
        
    @property
    def product_dataset(self):
        if isinstance(self._product_dataset, type(None)):
            poutg = self.poutg
            ds = read_L1(poutg[poutg.file_type == 'bl'].path2raw, parent = self,
                         ignore1=['name', 'message_type', 'version', 'date_stamp', 'period', 'tilt_angle', 'status_bits', 'profile_scale', 'profile_resolution', 'profile_length'])
            # self.tp_L1 = L1.copy()
            
            #### do we really need to go to pandas and back to xarray again?
            # dfL1 = L1.rcs_910.to_pandas()
            # assert(dfL1.index.duplicated().sum() == 0), "there are duplicates in L1's index, I would think this should be happening. if it does un-comment the following line"
            # dfL1 = dfL1[~dfL1.index.duplicated(keep='first')]
            # assert(dfL1.index.duplicated().sum() == 0), "haaaa, duplicates should all be gone"

            # his3 = read_level3_hist(poutg[poutg.file_type == 'hist'].path2raw, parent = self)

            ##### Clean and resample to 36s ###########################################    
            # resample to 36s even though L1 files are already at 36 sec because the 
            # time intervals on the L1 files from BL-View are sometimes off by a second. 
            # Keeping the limit at 1 step prevents the resample from repeating a nearby
            # data point to fill in large gaps in the data.
            # dfL1 = dfL1.resample('36S').nearest(limit=1)

            # The .his files are originally at 16 sec.
            # self.tp_his3 = his3.copy()
            # his3 = his3[his3.index.notna()] # sometimes the last row (or others?) has a NaT index, this removes it
            # his3 = his3.resample('36S').nearest(limit=1)

            # Do this to fill in an gaps in the data with nans to build com_plete days.
            # Create a date range of a com_plete day with 36s intervals with no gaps.
            # do something similar to put it onto 1 minut
            # day = _pd.date_range(dfL1.index[0].floor('D'), dfL1.index[-1].ceil('D'),freq='36S')
            # df = _pd.DataFrame(index=day[:-1])
            # df.index.rename('time', inplace=True)

            # Merge the date range from above with the dataframes to fill in any gaps 
            # left to com_plete a whole day of 36s intervals.
            # dfx = _pd.merge_ordered(df, dfL1, on='time') # For the L1 file
            # dfx.set_index('time', drop = True, inplace = True)

            # dfhis = _pd.merge_ordered(df, his3, on='time')
            # dfhis.set_index('time', drop = True, inplace = True)

            ##### Build the Variables and attributes ##################################
            
            
            # var = {} # Create empty variable dictionary.

            # L1 file
            # var['backscatter_profile']=(['time','range'], _np.float32(dfx.values),
            #                         {'long_name':'2-D ceilometer signal backscatter profile.',
            #                          'units':'10e-9 m^-1 sr^-1',
            #                          'comments':'Range-corrected-scattering'})

            # # Level3 .his files.
            # var['cloud_status']=(['time'], _np.float32(dfhis['CLOUD_STATUS']),
            #                      {'long_name':'Cloud detection status.',
            #                       'units':'1',
            #                       'flag_values':_np.float32([0,1,2,3,4]),
            #                       'flag_0':'No significant backscatter.',
            #                       'flag_1':'One cloud layer detected.',
            #                       'flag_2':'Two cloud layers detected.',
            #                       'flag_3':'Three cloud layers detected.',
            #                       'flag_4':'''Full obscuration/vertical visibility mode. First cloud_base will report vertical visibility and second cloud_base will report highest signal''',
            #                       'comments':'''When cloud_status=4 there is an optically thick cloud that obscures the signal. Therefore it is not possible to discern additional cloud layers above it so the vertical visibility and highest signal are reported instead.'''})

            # var['cloud_base']=(['time','cloud_layer'],
            #                    _np.float32(dfhis[['CLOUD_1','CLOUD_2','CLOUD_3']].values.tolist()),
            #                    {'long_name':'Cloud base heights.',
            #                     'units':'m',
            #                     'comments':'''A 2D array containing all three cloud bases at each timestep. -999 if no significant signal.''',
            #                     'cloud_base_1':'''First cloud base height or vertical visibility if cloud_status=4''',
            #                     'cloud_base_2':'''Second cloud base height or highest received signal if cloud_status=4''',
            #                     'cloud_base_3':'Third cloud base height'})   


            # ##### Create dataset. #####################################################
            # ds = _xr.Dataset(attrs = {'title':'Ceilometer cloud product',
            #                          'version':self.version,
            #                          'institution':'NOAA/GML/GRAD',
            #                          'author':'hagen.telg@noaa.gov',
            #                          'source':'Vaisala CL51 ceilometer',
            #                          'serial_number': poutg.sn.unique()[0],
            #                          'input_files': ', '.join([fn.name for fn in poutg.path2raw]),
            #                          'Conventions':'CF-1.8',
            #                          'comments':'''The data has not undergone any processing other than what the Vaisala software has applied. In addition, no QC has been applied other than a visual inspection for impossible values or obvious errors.'''
            #                          },
            #                 coords = {'time':('time', df.index.values,
            #                                   {'long_name':'Time in UTC',
            #                                   'comments':'36 second message interval'}),
            #                           'range':('range', L1['range'].values,
            #                                   {'long_name':'Vertical height bins',
            #                                   'units':'meters'}),                  
            #                           'cloud_layer':('cloud_layer', _np.array([1,2,3], dtype = _np.int8),
            #                                          {'long_name':'cloud layer',
            #                                           'units':'1'}),
            #                           },
            #                 data_vars = var
            #                 )
            
            # self._product_dataset = ds
            
            ####FIXME : I don't see why this should not be in here and not in the L1 reader?!?
            
            dsr = _xr.Dataset()
            dsr['backscatter_profile']=ds.rcs_910.astype(_np.float32)
            
            dsr.backscatter_profile.attrs = {'long_name':'2-D ceilometer signal backscatter profile.',
                                             'units':'10e-9 m^-1 sr^-1',
                                             'comments':'Range-corrected-scattering'}
            
            dsr['cloud_status'] = ds.cloud_status.isel(cloud_statusDim = 0).astype(_np.float32)
            
            dsr.cloud_status.attrs = {'long_name':'Cloud detection status.',
                                              'units':'1',
                                              'flag_values':'0,1,2,3,4',
                                              'flag_0':'No significant backscatter.',
                                              'flag_1':'One cloud layer detected.',
                                              'flag_2':'Two cloud layers detected.',
                                              'flag_3':'Three cloud layers detected.',
                                              'flag_4':'''Full obscuration/vertical visibility mode. First cloud_base will report vertical visibility and second cloud_base will report highest signal''',
                                              'comments':'''When cloud_status=4 there is an optically thick cloud that obscures the signal. Therefore it is not possible to discern additional cloud layers above it so the vertical visibility and highest signal are reported instead.'''}
            
            dst = ds.cloud_data.rename({'cloud_dataDim': 'cloud_layer'})
            dst = dst.assign_coords(cloud_layer = _np.int8(dst.cloud_layer + 1))
            dst = dst.where(dst > 0) #for The L1 it looks like 0 is invalid
            dsr['cloud_data'] = dst.astype(_np.float32)
            
            dsr.cloud_data.attrs = {'long_name':'Cloud base heights.',
                                            'units':'m',
                                            'comments':'''A 2D array containing all three cloud bases at each timestep. -999 if no significant signal.''',
                                            'cloud_base_1':'''First cloud base height or vertical visibility if cloud_status=4''',
                                            'cloud_base_2':'''Second cloud base height or highest received signal if cloud_status=4''',
                                            'cloud_base_3':'Third cloud base height'}
            
            dsr.attrs = {'title':'Ceilometer cloud product',
                         'version':self.version,
                         'institution':'NOAA/GML/GRAD',
                         'author':'hagen.telg@noaa.gov',
                         'source':'Vaisala CL51 ceilometer',
                         'serial_number': poutg.sn.unique()[0],
                         'input_files': ', '.join([fn.name for fn in poutg.path2raw]),
                         'Conventions':'CF-1.8',
                         'comments':''' The "time" coordinate was re-indexed to the full minute by back-filling the nearest 
valid data value within the following minute.The data values have not undergone any processing other than what the Vaisala 
software has applied. In addition, no QC has been applied other than a visual inspection for impossible values or obvious errors.'''.replace('\n','')
                         }
            
            #### attributes of the coordinates
            dsr.cloud_layer.attrs = dict(long_name = "Cloud layer index",
                                        units = "1" ,
                                        axis = "Z" ,
                                        positive = "up",)

            dsr.range.attrs = {'long_name': 'Distance from ground',
                               'units': 'm',
                               'standard_name' : 'distance',
                               'axis' : 'Z',
                               'positive': 'up'}
            #### reindex on 1min data and fill any gaps with nan ... every minute of the day has a timestamp!
            time = dsr.time.to_pandas()
            start = time.iloc[0].date()
            end = start + _pd.to_timedelta(1, 'd')
            newtime = _pd.date_range(start = start, end = end, freq = _pd.to_timedelta(1, 'minute'), inclusive = 'left')
            
            assert(time.iloc[0].date() == time.iloc[-1].date()), f'Fist and last timestamp are not on the same date!! {time.iloc[0].date()} vs {time.iloc[-1].date()}'
            dsr['instrumen_reported_time'] = dsr.time.copy() #this ensures we still have the original timestamp for each value
            dsr.instrumen_reported_time.attrs = {'long_name': "Instrument reported time",
                                                 'comments': 'Timestamp that was reported by the instrument. The "time" coordinate was re-indexed to the full minute (see global comments).'}
            dsr = dsr.reindex(time = newtime, method = 'bfill', tolerance = '1min')
            
            self._product_dataset = dsr
                        
            
        return self._product_dataset

class Cl51CloudProdProcessor_v1p2(object):
    def __init__(self, 
                 p2fl_in = '/nfs/grad/Inst/Ceil/SURFRAD/',
                 p2fl_out = '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1',
                 version = '0',
                 hist_file_format = '*_CEILOMETER_1_LEVEL_3*.his',
                 ignore = [],
                 verbose = False,
                 ):
        """
        This class aggreates all functions that are needed to create the cl51_cloud_prod_lev1 (formally known as cl51_cloud_prod_lev0) product.
        

        Parameters
        ----------
        p2fl_in : str, optional
            Path to input folder. This folder is expected to have further subfolders, one for each instrument (typically the measurment site, e.g. TBL)
        p2fl_out : str, optional
            Output folder. The default is '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1'.
        hist_file_format : TYPE, optional
            Format (pattern) of the hist files. As this format sometimes changes it not hardcoded. The default is '*_CEILOMETER_1_LEVEL_3*.his'.
        ignore : list, optional
            When there is a folder in the p2fl_in folder that mathes an element of this list it will be ignored.
        verbose : bool, optional
            If True there will be more output ... mostly for develpment purposes. The default is False.


        Returns
        -------
        Cl51CloudProdProcessor instance.

        """
        self.version = version
        self.ignore = ignore
        self.p2fl_in = _pl.Path(p2fl_in)
        self.hist_file_format= hist_file_format
        self.bl_file_format = 'L1*.nc'
        self.p2fl_out = _pl.Path(p2fl_out)
        self.fn_format_out = '{site}.cl51.cloud_prod.{date}.nc'
        self.verbose = verbose
        # self.test = test
        
        self._workplan = None
    
    
    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
            if self.verbose:
                print('creating workplan')
                print('=================')
            def bl2date(row):
                if row.path2raw.name.split('_')[-1].split('.')[0].isnumeric():
                    dt = row.path2raw.name.split('_')[-1].split('.')[0]
                else:
                    dt = row.path2raw.name.split('_')[-2]
                return _pd.to_datetime(dt)

            workplan = _pd.DataFrame()

            for p2site in self.p2fl_in.glob('*'):
                if not p2site.is_dir():
                    continue
                
                if self.verbose:
                    print(f'\t folder: {p2site}')
                if p2site.name in self.ignore:
                    if self.verbose:
                        print('\t\t ignore!')
                    continue

                # get the different file types
                ## hist files

                # p2fl_hist = p2site.joinpath('hist')
                # df = _pd.DataFrame(p2fl_hist.glob(self.hist_file_format), columns = ['path2raw'])
                # df['file_type'] = 'hist'
                # df.index = df.apply(lambda row: _pd.to_datetime('{}{}'.format(row.path2raw.name.split('_')[0], row.path2raw.name.split('_')[-1].split('.')[0])), axis = 1)
                # df_hist = df
                
                ## bl files
                p2fl_bl = p2site.joinpath('bl')
                df = _pd.DataFrame(p2fl_bl.glob(self.bl_file_format), columns = ['path2raw'])
                df['file_type'] = 'bl'
                df.index = df.apply(bl2date, axis = 1)
                # df_bl = df

                ## concat  

                # df = _pd.concat([df_hist, df_bl], sort = True)

                # some more collumns
                ## site
                df['site'] = p2site.name.lower()
                # add to workplan
                workplan = _pd.concat([workplan, df])
                

            workplan.sort_index(inplace=True)
            
            ## add datetime collumn
            workplan['date'] = workplan.apply(lambda row: _pd.to_datetime(row.name.date()), axis = 1)
            
            ## create path to output file

            def create_path_out(row):
                fname = self.fn_format_out.format(site = row.site, date = row.name.date().__str__().replace('-',''))
                p2fn = self.p2fl_out.joinpath(row.site).joinpath(f'{row.name.year}').joinpath(fname)
                return p2fn

            workplan['path2fn_out'] = workplan.apply(create_path_out, axis = 1)
            
            #### remove row when output file exists

            workplan = workplan[~ workplan.apply(lambda row: row.path2fn_out.is_file(), axis = 1)]

            #### remove last day
            
            workplan = workplan[~(workplan.date == workplan.date.iloc[-1]).values]
            
            #### test if bl and hist exist and remove if not
            #is it enough that bl and hist exists or do we have to check that it covers the entire day? We could execute this before deleting the last day ... would that be saver?
            
            #workplan = workplan[workplan['path2fn_out'].map(workplan.groupby('path2fn_out').apply(lambda group: group.file_type.eq('bl').any()&group.file_type.eq('hist').any()))]
                        
            self._workplan = workplan
            if self.verbose:
                print('========= workplan done ==========')
        return self._workplan
    
    @workplan.setter
    def workplan(self, value):
        self._workplan = value
        
    def process(self, test=False, 
                path2fn_out = None,
                generate_missing_folders = False,
                error_handling = 'raise',
                error_handling_serial='raise',
                verbose = False, 
                complevel = 4):
        """
        Processes the workplan

        Parameters
        ----------
        test : [bool, int], optional
            Creats several test scenarios:
                1: returns firt product dataset, stops processing afterwards
                2: as 1 but saves dataset
        path2fn_out: str, optional, default is None
            For testing. Excecute only files for this particular output path.
        error_handling: str, optional, default is 'return'
            What to do if an error accures during processing of single 
            retrieval.
            'raise': will raise the error
            'return': will not raise an errorj, but will return it. Check 
                'errors' key of return dict.
            

        Returns
        -------
        ds : TYPE
            DESCRIPTION.

        """
        if self.verbose:
            print('processing')
            print('==========')
        def make_all_parent_flds(file): 
            if not file.parent.is_dir():
                make_all_parent_flds(file.parent)
                file.parent.mkdir(mode = 0o775)
                file.parent.chmod(0o775) # for some reason the mode kwarg in the line above did not give the desired result?!?
                
        no_processed = 0
        errors = []
        errror_grp = []
        out = {}
        out['start_time'] = _pd.Timestamp(_dt.datetime.now())
        for p2fnout, poutg in self.workplan.groupby('path2fn_out'): 
            if verbose:
                print(f'\t path2fn_out: {p2fnout} - {poutg}')
            if not isinstance(path2fn_out, type(None)):
                if p2fnout != _pl.Path(path2fn_out):
                    continue
                              
            ds = None
            try:
                retriever = Cl51CloudProdRetriever_v1p2(poutg, self.version)
                self.tp_retriever = retriever
                retriever.check_serial(error_handling=error_handling_serial)
                ds =retriever.product_dataset
                no_processed += 1
                
            except Exception as err:
                if error_handling == 'raise':
                    raise
                
                #### error handling
                elif error_handling == 'return':
                    errors.append(err)
                    errror_grp.append(poutg)
                    # errors.append(sys.exc_info())
                    # error, error_txt, trace = sys.exc_info()
                    # tm = ['{}: {}'.format(error.__name__, error_txt.args[0])] + _tb.format_tb(trace)
                
            if test == 1:
                break
            
            #### Write to disk.        
            if not isinstance(ds, type(None)):
                assert(not p2fnout.is_file()), f'File exists, this should not be possible! path: {p2fnout}'
            
                if generate_missing_folders:
                    make_all_parent_flds(p2fnout)
                
                #### encodings
                encoding = {}
                for var in ds:
                    encoding[var] = {'zlib': True,     # Enable compression
                                     'complevel': complevel   # Moderate compression level
                                    }
                encoding['time'] = {'units': "minutes since 1900-01-01 00:00:00 UTC",
                                   'zlib': True,     # Enable compression
                                     'complevel': complevel   # Moderate compression level
                                    }
                encoding['instrumen_reported_time'] = {'units': "minutes since 1900-01-01 00:00:00 UTC",
                                                      'zlib': True,     # Enable compression
                                                         'complevel': complevel   # Moderate compression level
                                                        }
                
                #### save to netcdf
                ds.to_netcdf(path=p2fnout, 
                             # mode='w', 
                             encoding = encoding,
                             format='NETCDF4')
                p2fnout.chmod(0o664) # to give the group writing access.
            else:
                if verbose:
                    print(f'{p2fnout} was not saved as ds is None.')
                
            if test == 2:
                break
        
        out['errors'] = errors
        out['errror_grp'] = errror_grp
        out['last_dataset'] = ds
        out['no_of_files_processed'] = no_processed
        out['finish_time'] = _pd.Timestamp(_dt.datetime.now())
        self._last_processing = out
        return out
    
    def get_single_day_from_worplan(self, index = -1, random = False):
        gl = [g for idc,g in self.workplan.groupby('path2fn_out')]
        
        if random:
            index = int(_np.round(_np.random.random(1) * (len(gl)-1)))
    
        oneday = gl[index]
        return oneday
    
    def notify(self, subject = 'cl51cloudprod - status: {status} (clean: {no_of_files_processed}; errors: {no_of_errors})'):
        config = load_config()
        assert(config.get('notify', 'email_address') != 'None'), 'No email has been specified, _please do so in ~/.ceilopy/config.ini'
        
        out = self._last_processing
        messages = ['run started {}'.format(out['start_time'])]
        messages.append('run finshed {}'.format(out['finish_time']))
        
        #### summary
        no_of_errors = len(out['errors'])
        no_of_files_processed = out['no_of_files_processed']
            
        messages.append('\n'.join([f'length of workplan:\t\t{self.workplan.shape[0]}',f'successfully created files:\t{no_of_files_processed}', f'errors encountered:\t\t{no_of_errors}']))
        
        #### errors
        if no_of_errors != 0:
            errs = out['errors']
            err_types = _np.array([err.__class__.__name__ for err in errs])
            
            typoferrors = []
            for toe in _np.unique(err_types):
                typoferrors.append({'name': toe, 'times': (err_types == toe).sum(), 'first_err': errs[_np.where(err_types == toe)[0][0]]})
        
            messages.append('\n'.join(['Errors by type:',] + ['\t{tn}:\t{tt}'.format(tn = toe['name'], tt = toe['times']) for toe in typoferrors]))
            messages.append('\n=============================================\n'.join(['First traceback for each error type',]+[''.join(_tb.format_tb(toe['first_err'].__traceback__) + [toe['first_err'].__str__(),]) for toe in typoferrors]))
        
        #### email body
        message_txt = '\n=========================================================================\n'.join(messages)
        msg = _MIMEText(message_txt)
        
        #### subject
        if no_of_errors ==0:
            status = 'clean'
        else:
            status = 'errors'
        # subject = f'cl51cloudprod - status: {status} (clean: {no_of_files_processed}; errors: {no_of_errors})'
        subject = subject.format(status = status, no_of_files_processed = no_of_files_processed, no_of_errors = no_of_errors)
        address  = config.get('notify', 'email_address')
        smtp = config.get('notify', 'smtp')
        
        msg['Subject'] = subject
        msg['From'] = address
        msg['To'] = address
        
        # Send the message via our own SMTP server.
        s = _smtplib.SMTP(smtp)
        s.send_message(msg)
        s.quit()
        
#######
#### This is Level 1
###FIXME dont want to change the classname because it is still used by the operational SURFRAD retrieval  

class Cl51CloudProdRetriever():
    def __init__(self, poutg, version,
                 # check_serial = True,
                 ):
        self.version = version
        self.poutg = poutg
#         self.p2fnout = poutg.path2fn_out.unique()[0]
        self._product_dataset = None
        self.get_serial_numbers()
        # if check_serial:
        #     self.check_serial()
        
    def get_serial_numbers(self):
        def get_serial(row):
            # Extract serial numbers from files
            key = row.file_type
            file = row.path2raw
            if key in ['L1', 'L2', 'bl']:
                serial = file.name[-11:-3]  # Extract serial number from L1 filename.
            # elif key == 'L3':
            #     serial = files['L3'][-11:-3]
            elif key in ['H2','H3','hist']:
                h = _pd.read_csv(file, skiprows=1, header=0, sep=',')
                serial = h[' CEILOMETER'][0].strip() # Extract serial number from H2 file.
            else:
                raise KeyError('File type unknown')
            return serial
        
        self.poutg['sn'] = self.poutg.apply(get_serial, axis = 1)
        
    def check_serial(self, error_handling = 'raise'):
        """
        Checks if the serial numbers in all the files are the same. In early 
        measurments the serial number was not stored ... use error_handling to
        deal with occuring errors.

        Parameters
        ----------
        error_handling : str, optional
            How to deal with errors. The default is 'raise'.
            raise: raises occuring errors
            allow_empty: do not raise an error if serial number is not available

        Raises
        ------
        KeyError
            DESCRIPTION.

        Returns
        -------
        serial : TYPE
            DESCRIPTION.

        """
        sn_series = self.poutg['sn'].copy()
        # self.poutg['sn'] = sn_series.copy()
        valid = ['raise', 'allow_empty']
        assert(error_handling in valid), f'error_handling got an unexpected value ({error_handling}. Choose from: {valid})'
        if error_handling == 'allow_empty':
            sn_series = sn_series[sn_series.apply(lambda x: len(x)) != 0]
        if sn_series.unique().shape[0] != 1:
            if len(sn_series[sn_series.apply(lambda x: len(x)) != 0]) != len(sn_series):
                fnj = '\n\t'.join([fn.as_posix() for fn in self.poutg.path2raw])
                raise MissingSerialNumberError(f'At least one of the following files is missing a serial number:\n\t{fnj}')
            raise SerialNumberMissmatchError(f'Serial numbers ({sn_series.unique()}) do not match')
        
    @property
    def product_dataset(self):
        if isinstance(self._product_dataset, type(None)):
            poutg = self.poutg
            L1 = read_L1(poutg[poutg.file_type == 'bl'].path2raw, parent = self)
            self.tp_L1 = L1.copy()
            dfL1 = L1.rcs_910.to_pandas()
            # assert(dfL1.index.duplicated().sum() == 0), "there are duplicates in L1's index, I would think this should be happening. if it does un-comment the following line"
            dfL1 = dfL1[~dfL1.index.duplicated(keep='first')]
            assert(dfL1.index.duplicated().sum() == 0), "haaaa, duplicates should all be gone"

            his3 = read_level3_hist(poutg[poutg.file_type == 'hist'].path2raw, parent = self)

            ##### Clean and resample to 36s ###########################################    
            # resample to 36s even though L1 files are already at 36 sec because the 
            # time intervals on the L1 files from BL-View are sometimes off by a second. 
            # Keeping the limit at 1 step prevents the resample from repeating a nearby
            # data point to fill in large gaps in the data.
            dfL1 = dfL1.resample('36S').nearest(limit=1)

            # The .his files are originally at 16 sec.
            self.tp_his3 = his3.copy()
            his3 = his3[his3.index.notna()] # sometimes the last row (or others?) has a NaT index, this removes it
            his3 = his3.resample('36S').nearest(limit=1)

            # Do this to fill in an gaps in the data with nans to build com_plete days.
            # Create a date range of a com_plete day with 36s intervals with no gaps.
            day = _pd.date_range(dfL1.index[0].floor('D'), dfL1.index[-1].ceil('D'),freq='36S')
            df = _pd.DataFrame(index=day[:-1])
            df.index.rename('time', inplace=True)

            # Merge the date range from above with the dataframes to fill in any gaps 
            # left to com_plete a whole day of 36s intervals.
            dfx = _pd.merge_ordered(df, dfL1, on='time') # For the L1 file
            dfx.set_index('time', drop = True, inplace = True)

            dfhis = _pd.merge_ordered(df, his3, on='time')
            dfhis.set_index('time', drop = True, inplace = True)

            ##### Build the Variables and attributes ##################################
            var = {} # Create empty variable dictionary.

            # L1 file
            var['backscatter_profile']=(['time','range'], _np.float32(dfx.values),
                                    {'long_name':'2-D ceilometer signal backscatter profile.',
                                     'units':'10e-9 m^-1 sr^-1',
                                     'comments':'Range-corrected-scattering'})

            # Level3 .his files.
            var['cloud_status']=(['time'], _np.float32(dfhis['CLOUD_STATUS']),
                                 {'long_name':'Cloud detection status.',
                                  'units':'1',
                                  'flag_values':_np.float32([0,1,2,3,4]),
                                  'flag_0':'No significant backscatter.',
                                  'flag_1':'One cloud layer detected.',
                                  'flag_2':'Two cloud layers detected.',
                                  'flag_3':'Three cloud layers detected.',
                                  'flag_4':'''Full obscuration/vertical visibility mode. First cloud_base will report vertical visibility and second cloud_base will report highest signal''',
                                  'comments':'''When cloud_status=4 there is an optically thick cloud that obscures the signal. Therefore it is not possible to discern additional cloud layers above it so the vertical visibility and highest signal are reported instead.'''})

            var['cloud_base']=(['time','cloud_layer'],
                               _np.float32(dfhis[['CLOUD_1','CLOUD_2','CLOUD_3']].values.tolist()),
                               {'long_name':'Cloud base heights.',
                                'units':'m',
                                'comments':'''A 2D array containing all three cloud bases at each timestep. -999 if no significant signal.''',
                                'cloud_base_1':'''First cloud base height or vertical visibility if cloud_status=4''',
                                'cloud_base_2':'''Second cloud base height or highest received signal if cloud_status=4''',
                                'cloud_base_3':'Third cloud base height'})   


            ##### Create dataset. #####################################################
            ds = _xr.Dataset(attrs = {'title':'Ceilometer cloud product',
                                     'version':self.version,
                                     'institution':'NOAA/GML/GRAD',
                                     'author':'hagen.telg@noaa.gov',
                                     'source':'Vaisala CL51 ceilometer',
                                     'serial_number': poutg.sn.unique()[0],
                                     'input_files': ', '.join([fn.name for fn in poutg.path2raw]),
                                     'Conventions':'CF-1.8',
                                     'comments':'''The data has not undergone any processing other than what the Vaisala software has applied. In addition, no QC has been applied other than a visual inspection for impossible values or obvious errors.'''
                                     },
                            coords = {'time':('time', df.index.values,
                                              {'long_name':'Time in UTC',
                                              'comments':'36 second message interval'}),
                                      'range':('range', L1['range'].values,
                                              {'long_name':'Vertical height bins',
                                              'units':'meters'}),                  
                                      'cloud_layer':('cloud_layer', _np.array([1,2,3], dtype = _np.int8),
                                                     {'long_name':'cloud layer',
                                                      'units':'1'}),
                                      },
                            data_vars = var
                            )
            
            self._product_dataset = ds
            
        return self._product_dataset

class Cl51CloudProdProcessor(object):
    def __init__(self, 
                 p2fl_in = '/nfs/grad/Inst/Ceil/SURFRAD/',
                 p2fl_out = '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1',
                 version = '0',
                 hist_file_format = '*_CEILOMETER_1_LEVEL_3*.his',
                 ignore = [],
                 verbose = False,
                 ):
        """
        This class aggreates all functions that are needed to create the cl51_cloud_prod_lev1 (formally known as cl51_cloud_prod_lev0) product.
        

        Parameters
        ----------
        p2fl_in : str, optional
            Path to input folder. This folder is expected to have further subfolders, one for each instrument (typically the measurment site, e.g. TBL)
        p2fl_out : str, optional
            Output folder. The default is '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1'.
        hist_file_format : TYPE, optional
            Format (pattern) of the hist files. As this format sometimes changes it not hardcoded. The default is '*_CEILOMETER_1_LEVEL_3*.his'.
        ignore : list, optional
            When there is a folder in the p2fl_in folder that mathes an element of this list it will be ignored.
        verbose : bool, optional
            If True there will be more output ... mostly for develpment purposes. The default is False.


        Returns
        -------
        Cl51CloudProdProcessor instance.

        """
        self.version = version
        self.ignore = ignore
        self.p2fl_in = _pl.Path(p2fl_in)
        self.hist_file_format= hist_file_format
        self.bl_file_format = 'L1*.nc'
        self.p2fl_out = _pl.Path(p2fl_out)
        self.fn_format_out = '{site}.cl51.cloud_prod.{date}.nc'
        self.verbose = verbose
        # self.test = test
        
        self._workplan = None
        
    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
            if self.verbose:
                print('creating workplan')
                print('=================')
            def bl2date(row):
                if row.path2raw.name.split('_')[-1].split('.')[0].isnumeric():
                    dt = row.path2raw.name.split('_')[-1].split('.')[0]
                else:
                    dt = row.path2raw.name.split('_')[-2]
                return _pd.to_datetime(dt)

            workplan = _pd.DataFrame()

            for p2site in self.p2fl_in.glob('*'):
                if not p2site.is_dir():
                    continue
                
                if self.verbose:
                    print(f'\t folder: {p2site}')
                if p2site.name in self.ignore:
                    if self.verbose:
                        print('\t\t ignore!')
                    continue

                # get the different file types
                ## hist files

                p2fl_hist = p2site.joinpath('hist')
                df = _pd.DataFrame(p2fl_hist.glob(self.hist_file_format), columns = ['path2raw'])
                df['file_type'] = 'hist'
                df.index = df.apply(lambda row: _pd.to_datetime('{}{}'.format(row.path2raw.name.split('_')[0], row.path2raw.name.split('_')[-1].split('.')[0])), axis = 1)
                df_hist = df

                ## bl files
                p2fl_bl = p2site.joinpath('bl')
                df = _pd.DataFrame(p2fl_bl.glob(self.bl_file_format), columns = ['path2raw'])
                df['file_type'] = 'bl'
                df.index = df.apply(bl2date, axis = 1)
                df_bl = df

                ## concat  

                df = _pd.concat([df_hist, df_bl], sort = True)

                # some more collumns
                ## site
                df['site'] = p2site.name.lower()

                # add to workplan
                workplan = _pd.concat([workplan, df])


            workplan.sort_index(inplace=True)
            
            ## add datetime collumn
            workplan['date'] = workplan.apply(lambda row: _pd.to_datetime(row.name.date()), axis = 1)
            
            ## create path to output file

            def create_path_out(row):
                fname = self.fn_format_out.format(site = row.site, date = row.name.date().__str__().replace('-',''))
                p2fn = self.p2fl_out.joinpath(row.site).joinpath(f'{row.name.year}').joinpath(fname)
                return p2fn

            workplan['path2fn_out'] = workplan.apply(create_path_out, axis = 1)
            
            #### remove row when output file exists

            workplan = workplan[~ workplan.apply(lambda row: row.path2fn_out.is_file(), axis = 1)]

            #### remove last day
            
            workplan = workplan[~(workplan.date == workplan.date.iloc[-1]).values]
            
            #### test if bl and hist exist and remove if not
            #is it enough that bl and hist exists or do we have to check that it covers the entire day? We could execute this before deleting the last day ... would that be saver?
            
            workplan = workplan[workplan['path2fn_out'].map(workplan.groupby('path2fn_out').apply(lambda group: group.file_type.eq('bl').any()&group.file_type.eq('hist').any()))]
                        
            self._workplan = workplan
            if self.verbose:
                print('========= workplan done ==========')
        return self._workplan
    
    @workplan.setter
    def workplan(self, value):
        self._workplan = value
        
    def process(self, test=False, 
                path2fn_out = None,
                generate_missing_folders = False,
                error_handling = 'raise',
                error_handling_serial='raise',
                verbose = False):
        """
        Processes the workplan

        Parameters
        ----------
        test : [bool, int], optional
            Creats several test scenarios:
                1: returns firt product dataset, stops processing afterwards
                2: as 1 but saves dataset
        path2fn_out: str, optional, default is None
            For testing. Excecute only files for this particular output path.
        error_handling: str, optional, default is 'return'
            What to do if an error accures during processing of single 
            retrieval.
            'raise': will raise the error
            'return': will not raise an errorj, but will return it. Check 
                'errors' key of return dict.
            

        Returns
        -------
        ds : TYPE
            DESCRIPTION.

        """
        if self.verbose:
            print('processing')
            print('==========')
        def make_all_parent_flds(file): 
            if not file.parent.is_dir():
                make_all_parent_flds(file.parent)
                file.parent.mkdir(mode = 0o775)
                file.parent.chmod(0o775) # for some reason the mode kwarg in the line above did not give the desired result?!?
                
        no_processed = 0
        errors = []
        errror_grp = []
        out = {}
        out['start_time'] = _pd.Timestamp(_dt.datetime.now())
        for p2fnout, poutg in self.workplan.groupby('path2fn_out'): 
            if verbose:
                print(f'\t path2fn_out: {p2fnout} - {poutg}')
            if not isinstance(path2fn_out, type(None)):
                if p2fnout != _pl.Path(path2fn_out):
                    continue
                              
            ds = None
            try:
                retriever = Cl51CloudProdRetriever(poutg, self.version)
                self.tp_retriever = retriever
                retriever.check_serial(error_handling=error_handling_serial)
                ds =retriever.product_dataset
                no_processed += 1
                
            except Exception as err:
                if error_handling == 'raise':
                    raise
                
                #### error handling
                elif error_handling == 'return':
                    errors.append(err)
                    errror_grp.append(poutg)
                    # errors.append(sys.exc_info())
                    # error, error_txt, trace = sys.exc_info()
                    # tm = ['{}: {}'.format(error.__name__, error_txt.args[0])] + _tb.format_tb(trace)
                
            if test == 1:
                break
            
            #### Write to disk.        
            if not isinstance(ds, type(None)):
                assert(not p2fnout.is_file()), f'File exists, this should not be possible! path: {p2fnout}'
            
                if generate_missing_folders:
                    make_all_parent_flds(p2fnout)
                    
                ds.to_netcdf(path=p2fnout, mode='w', format='NETCDF4')
                p2fnout.chmod(0o664) # to give the group writing access.
            else:
                if verbose:
                    print(f'{p2fnout} was not saved as ds is None.')
                
            if test == 2:
                break
        
        out['errors'] = errors
        out['errror_grp'] = errror_grp
        out['last_dataset'] = ds
        out['no_of_files_processed'] = no_processed
        out['finish_time'] = _pd.Timestamp(_dt.datetime.now())
        self._last_processing = out
        return out
    
    def get_single_day_from_worplan(self, index = -1, random = False):
        gl = [g for idc,g in self.workplan.groupby('path2fn_out')]
        
        if random:
            index = int(_np.round(_np.random.random(1) * (len(gl)-1)))
    
        oneday = gl[index]
        return oneday
    
    def notify(self, subject = 'cl51cloudprod - status: {status} (clean: {no_of_files_processed}; errors: {no_of_errors})'):
        config = load_config()
        assert(config.get('notify', 'email_address') != 'None'), 'No email has been specified, _please do so in ~/.ceilopy/config.ini'
        
        out = self._last_processing
        messages = ['run started {}'.format(out['start_time'])]
        messages.append('run finshed {}'.format(out['finish_time']))
        
        #### summary
        no_of_errors = len(out['errors'])
        no_of_files_processed = out['no_of_files_processed']
            
        messages.append('\n'.join([f'length of workplan:\t\t{self.workplan.shape[0]}',f'successfully created files:\t{no_of_files_processed}', f'errors encountered:\t\t{no_of_errors}']))
        
        #### errors
        if no_of_errors != 0:
            errs = out['errors']
            err_types = _np.array([err.__class__.__name__ for err in errs])
            
            typoferrors = []
            for toe in _np.unique(err_types):
                typoferrors.append({'name': toe, 'times': (err_types == toe).sum(), 'first_err': errs[_np.where(err_types == toe)[0][0]]})
        
            messages.append('\n'.join(['Errors by type:',] + ['\t{tn}:\t{tt}'.format(tn = toe['name'], tt = toe['times']) for toe in typoferrors]))
            messages.append('\n=============================================\n'.join(['First traceback for each error type',]+[''.join(_tb.format_tb(toe['first_err'].__traceback__) + [toe['first_err'].__str__(),]) for toe in typoferrors]))
        
        #### email body
        message_txt = '\n=========================================================================\n'.join(messages)
        msg = _MIMEText(message_txt)
        
        #### subject
        if no_of_errors ==0:
            status = 'clean'
        else:
            status = 'errors'
        # subject = f'cl51cloudprod - status: {status} (clean: {no_of_files_processed}; errors: {no_of_errors})'
        subject = subject.format(status = status, no_of_files_processed = no_of_files_processed, no_of_errors = no_of_errors)
        address  = config.get('notify', 'email_address')
        smtp = config.get('notify', 'smtp')
        
        msg['Subject'] = subject
        msg['From'] = address
        msg['To'] = address
        
        # Send the message via our own SMTP server.
        s = _smtplib.SMTP(smtp)
        s.send_message(msg)
        s.quit()
        
