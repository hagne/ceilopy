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
import magic as _magic



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

def read_L1(file, parent = None):
    if isinstance(file, (str, _pl.Path)):
        file = [file]
        
    assert(isinstance(file, (_pd.Series,list, _np.array))), f'File type not recognized: {type(file)}'
    
    ignore1 = ['name','message_type','version','date_stamp',
               'period','tilt_angle',
              'cloud_status','cloud_data','status_bits','profile_scale',
              'profile_resolution','profile_length']
    
    if not _np.all([_magic.from_file(fn.as_posix()) == 'Hierarchical Data Format (version 5) data' for fn in file]):
        fnc = '\n\t'.join([fn.as_posix() for fn in file])
        raise CorruptFileError(f'At least one of the following can not be identified as a netcdf file: \n\t {fnc}')
        
    L1 = _xr.open_mfdataset(file, concat_dim = 'timeDim', drop_variables=ignore1)
    L1 = L1.assign_coords(time = _pd.to_datetime(L1.time.values, unit = 's'))
    for var in L1.variables:
        if 'timeDim' in L1[var].dims:
            L1[var] = L1[var].swap_dims({'timeDim':'time'})
    return L1

# read hist file
##### Read Level3 hist files. #############################################
def read_level3_hist(file, parent = None):
    def read_file(fn):
        cols = ['CREATEDATE',' CEILOMETER',' CLOUD_STATUS',' CLOUD_1',' CLOUD_2',
                ' CLOUD_3'] # What columns to keep.
        his3 = _pd.read_csv(fn, skiprows=1, header=0, sep=',',
                           na_values='-9999', index_col=0, parse_dates=True,
                           infer_datetime_format=True, usecols=cols)
        his3.index.rename('time', inplace=True)  
        his3.columns = [col.strip() for col in his3.columns]
        return his3
    
    if isinstance(file, (str, _pl.Path)):
        file = [file]
        
    assert(isinstance(file, (_pd.Series,list, _np.array))), f'File type not recognized: {type(file)}'
    df = _pd.concat([read_file(fn) for fn in file], sort = True)
    #### testpoint
    parent.tp_dfcc = df.copy()
    # assert(df.index.duplicated().sum() == 0), 'There are duplicates in the hist file ... I would think this should not be happening. if it does un-comment the following line'
    df = df[~df.index.duplicated(keep='first')] # Remove duplicates
    return df
        
class Cl51CloudProdRetriever():
    def __init__(self, poutg, 
                 # check_serial = True,
                 ):
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
            
            L1 = read_L1(poutg[poutg.file_type == 'bl'].path2raw)

            dfL1 = L1.rcs_910.to_pandas()
            assert(dfL1.index.duplicated().sum() == 0), "there are duplicates in L1's index, I would think this should be happening. if it does un-comment the following line"
            # dfL1 = dfL1[~dfL1.index.duplicated(keep='first')]

            his3 = read_level3_hist(poutg[poutg.file_type == 'hist'].path2raw, parent = self)

            ##### Clean and resample to 36s ###########################################    
            # resample to 36s even though L1 files are already at 36 sec because the 
            # time intervals on the L1 files from BL-View are sometimes off by a second. 
            # Keeping the limit at 1 step prevents the resample from repeating a nearby
            # data point to fill in large gaps in the data.
            dfL1 = dfL1.resample('36S').nearest(limit=1)

            # The .his files are originally at 16 sec.
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
                                  'flag_4':'''Full obscuration/vertical visibility mode.
                                            First cloud_base will report vertical visibility
                                            and second cloud_base will report highest signal''',
                                  'comments':'''When cloud_status=4 there is an optically
                                              thick cloud that obscures the signal. Therefore
                                              it is not possible to discern additional cloud
                                              layers above it so the vertical visibility and
                                              highest signal are reported instead.'''})

            var['cloud_base']=(['time','cloud_layer'],
                               _np.float32(dfhis[['CLOUD_1','CLOUD_2','CLOUD_3']].values.tolist()),
                               {'long_name':'Cloud base heights.',
                                'units':'m',
                                'comments':'''A 2D array containing all three cloud bases
                                            at each timestep. -999 if no significant signal.''',
                                'cloud_base_1':'''First cloud base height or vertical visibility
                                                if cloud_status=4''',
                                'cloud_base_2':'''Second cloud base height or highest received
                                                signal if cloud_status=4''',
                                'cloud_base_3':'Third cloud base height'})   


            ##### Create dataset. #####################################################
            ds = _xr.Dataset(attrs = {'title':'Ceilometer cloud product',
                                     'version':'1.0',
                                     'institution':'NOAA/GML/GRAD',
                                     'author':'christian.herrera@noaa.gov',
                                     'source':'Vaisala CL51 ceilometer',
                                     'serial_number': poutg.sn.unique()[0],
                                     'i_nput_files': [fn.name for fn in poutg.path2raw],
                                     'Conventions':'CF-1.8',
                                     'comments':'''The data has not undergone any processing
                                         other than what the Vaisala software has ap_plied. 
                                         In addition, no QC has been ap_plied other than 
                                         a visual inspection for impossible values or 
                                         obvious errors.'''
                                     },
                            coords = {'time':('time', df.index.values,
                                              {'long_name':'Time in UTC',
                                              'comments':'36 second message interval'}),
                                      'range':('range', L1['range'].values,
                                              {'long_name':'Vertical height bins',
                                              'units':'meters'}),                  
                                      'cloud_layer':('cloud_layer', _np.array([1,2,3]),
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
                 p2fl_out = '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev0', 
                 ):
        
        self.p2fl_in = _pl.Path(p2fl_in)
        self.hist_file_format= '*LEVEL_3*.his'
        self.bl_file_format = 'L1*.nc'
        self.p2fl_out = _pl.Path(p2fl_out)
        self.fn_format_out = '{site}.cl51.cloud_prod.{date}.nc'
        # self.test = test
        
        self._workplan = None
        
    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
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
            #     if p2site.name != 'TBL':
            #         continue
            #     print(p2site)

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
        return self._workplan
    
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
        def make_all_parent_flds(file): 
            if not file.parent.is_dir():
                make_all_parent_flds(file.parent)
                file.parent.mkdir(mode = 0o775)
                file.parent.chmod(0o775) # for some reason the mode kwarg in the line above did not give the desired result?!?
                
        no_processed = 0
        errors = []
        out = {}
        out['start_time'] = _pd.Timestamp(_dt.datetime.now())
        for p2fnout, poutg in self.workplan.groupby('path2fn_out'): 
            if not isinstance(path2fn_out, type(None)):
                if p2fnout != _pl.Path(path2fn_out):
                    continue
                              
            ds = None
            try:
                retriever = Cl51CloudProdRetriever(poutg)
                retriever.check_serial(error_handling=error_handling_serial)
                ds =retriever.product_dataset
                no_processed += 1
                
            except Exception as err:
                if error_handling == 'raise':
                    raise
                
                #### error handling
                elif error_handling == 'return':
                    errors.append(err)
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
    
    def notify(self):
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
        subject = f'cl51cloudprod - status: {status} (clean: {no_of_files_processed}; errors: {no_of_errors})'
        address  = config.get('notify', 'email_address')
        smtp = config.get('notify', 'smtp')
        
        msg['Subject'] = subject
        msg['From'] = address
        msg['To'] = address
        
        # Send the message via our own SMTP server.
        s = _smtplib.SMTP(smtp)
        s.send_message(msg)
        s.quit()
        
