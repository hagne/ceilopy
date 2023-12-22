#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 12:09:24 2023

@author: hagen
"""
import xarray as _xr
import pandas as _pd
import numpy as _np
import pathlib as _pl
import traceback as _tb
import datetime as _dt
from email.mime.text import MIMEText as _MIMEText
import smtplib as _smtplib
# import pathlib as __pl
# import configparser as _cp
# import magic as _magic
# import warnings as _warnings
import ceilopy.file_io as fio
import ceilopy.ceilolab as ceilolab
import matplotlib.pyplot as _plt
import gc as _gc

class Cl51CloudProdRetriever_v1p3():
    def __init__(self, poutg, version,reporter,
                 # check_serial = True,
                 ):
        """
        There was a time when the product was a aggregate of different functions, hist and netcdf. 
        In todays version this might not be needed anymore... Consider consoidating this into the processor class unless?!?

        Parameters
        ----------
        poutg : TYPE
            DESCRIPTION.
        version : TYPE
            DESCRIPTION.
        # check_serial : TYPE, optional
            DESCRIPTION. The default is True.
         : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.version = version
        self.poutg = poutg
#         self.p2fnout = poutg.path2fn_out.unique()[0]
        self._product_dataset = None
        self.reporter = reporter
        # self.get_serial_numbers()
        # if check_serial:
        #     self.check_serial()
        
    # def get_serial_numbers(self):
    #     def get_serial(row):
    #         # Extract serial numbers from files
    #         key = row.file_type
    #         file = row.path2raw
    #         if key in ['L1', 'L2', 'bl']:
    #             serial = file.name[-11:-3]  # Extract serial number from L1 filename.
    #         # elif key == 'L3':
    #         #     serial = files['L3'][-11:-3]
    #         elif key in ['H2','H3','hist']:
    #             h = _pd.read_csv(file, skiprows=1, header=0, sep=',')
    #             serial = h[' CEILOMETER'][0].strip() # Extract serial number from H2 file.
    #         else:
    #             raise KeyError('File type unknown')
    #         return serial
        
    #     self.poutg['sn'] = self.poutg.apply(get_serial, axis = 1)
        
    # def check_serial(self, error_handling = 'raise'):
    #     """
    #     Currently not used anyChecks if the serial numbers in all the files are the same. In early 
    #     measurments the serial number was not stored ... use error_handling to
    #     deal with occuring errors.

    #     Parameters
    #     ----------
    #     error_handling : str, optional
    #         How to deal with errors. The default is 'raise'.
    #         raise: raises occuring errors
    #         allow_empty: do not raise an error if serial number is not available

    #     Raises
    #     ------
    #     KeyError
    #         DESCRIPTION.

    #     Returns
    #     -------
    #     serial : TYPE
    #         DESCRIPTION.

    #     """
    #     sn_series = self.poutg['sn'].copy()
    #     # self.poutg['sn'] = sn_series.copy()
    #     valid = ['raise', 'allow_empty']
    #     assert(error_handling in valid), f'error_handling got an unexpected value ({error_handling}. Choose from: {valid})'
    #     if error_handling == 'allow_empty':
    #         sn_series = sn_series[sn_series.apply(lambda x: len(x)) != 0]
    #     if sn_series.unique().shape[0] != 1:
    #         if len(sn_series[sn_series.apply(lambda x: len(x)) != 0]) != len(sn_series):
    #             fnj = '\n\t'.join([fn.as_posix() for fn in self.poutg.path2raw])
    #             raise MissingSerialNumberError(f'At least one of the following files is missing a serial number:\n\t{fnj}')
    #         raise SerialNumberMissmatchError(f'Serial numbers ({sn_series.unique()}) do not match')
        
    @property
    def product_dataset(self):
        if isinstance(self._product_dataset, type(None)):
            poutg = self.poutg
            
            try:
                dsl = []
                input_files = []
                for p in poutg.path2raw:
                    dst = fio.read_netcdf_level2(p)
                    dsl.append(dst.dataset)
                    input_files.append(p.name)
                ds = ceilolab.CeilometerData(_xr.concat(dsl, 'time'))
            except Exception as e:                
                if isinstance(e, OSError): #netcdf is not readable
                    if not (e.args[0] == -101) and (e.args[1] == 'NetCDF: HDF error'):
                        raise
                    # print('test not just for the error but also the argument of the error ... let it raise one time to get the message')
                    # raise
                    # pass
                
                # older L2 nc files have a strange file format, for now, just use hist files instead
                # valid for next exceptions
                elif isinstance(e, RuntimeError):
                    if e.args[0] != "Failed to decode variable 'name': NetCDF: HDF error":
                        raise
                elif isinstance(e, ValueError):
                    if e.args[0] == "Dimensions {'cloud_statusDim'} do not exist. Expected one or more of ('time', 'cloud_status')":
                        pass
                    elif "did not find a match in any of xarray's currently installed IO backends" in e.args[0]:
                        pass
                    else:
                        raise
                else:
                    raise
                
                if not isinstance(self.reporter, type(None)): self.reporter.warnings_increment()
                
                dsl = []
                input_files = []
                
                #### load hist lev 2 file
                files = poutg.path2althist_l2.unique()
                assert(files.shape[0] ==1), 'I thought there is always just 1 hist file a day, even if there is a reboot ... unlike the netcdffiles'
                dsh2 = fio.read_hist_level2(files[0])
                input_files.append(files[0].name)
                
                #### load hist lev 3 file
                files = poutg.path2althist_l3.unique()
                assert(files.shape[0] ==1), 'I thought there is always just 1 hist file a day, even if there is a reboot ... unlike the netcdffiles'
                dsh3 = fio.read_hist_level3(files[0])
                input_files.append(files[0].name)    
                
                ds = _xr.merge([dsh2.dataset, dsh3.dataset])
                ds = ceilolab.CeilometerData(ds)
                # raise

            ds = ds.decorate_dataset(version=self.version, site_name = poutg.site.iloc[0], serial_no=poutg.serial_no.iloc[0], input_files=input_files)
            ds = ds.remove_timeinconsitancies()
            
            if len(ds.log) > 0:
                print('--------')
                print(poutg.iloc[0].name)
                print(ds.log)
                
            ds = ds.reindex()
            
            
            
            
            
#             ds = read_L1(poutg[poutg.file_type == 'bl'].path2raw, parent = self,
#                          ignore1=['name', 'message_type', 'version', 'date_stamp', 'period', 'tilt_angle', 'status_bits', 'profile_scale', 'profile_resolution', 'profile_length'])

            
#             dsr = _xr.Dataset()
#             dsr['backscatter_profile']=ds.rcs_910.astype(_np.float32)
            
#             dsr.backscatter_profile.attrs = {'long_name':'2-D ceilometer signal backscatter profile.',
#                                              'units':'10e-9 m^-1 sr^-1',
#                                              'comments':'Range-corrected-scattering'}
            
#             dsr['cloud_status'] = ds.cloud_status.isel(cloud_statusDim = 0).astype(_np.float32)
            
#             dsr.cloud_status.attrs = {'long_name':'Cloud detection status.',
#                                               'units':'1',
#                                               'flag_values':'0,1,2,3,4',
#                                               'flag_0':'No significant backscatter.',
#                                               'flag_1':'One cloud layer detected.',
#                                               'flag_2':'Two cloud layers detected.',
#                                               'flag_3':'Three cloud layers detected.',
#                                               'flag_4':'''Full obscuration/vertical visibility mode. First cloud_base will report vertical visibility and second cloud_base will report highest signal''',
#                                               'comments':'''When cloud_status=4 there is an optically thick cloud that obscures the signal. Therefore it is not possible to discern additional cloud layers above it so the vertical visibility and highest signal are reported instead.'''}
            
#             dst = ds.cloud_data.rename({'cloud_dataDim': 'cloud_layer'})
#             dst = dst.assign_coords(cloud_layer = _np.int8(dst.cloud_layer + 1))
#             dst = dst.where(dst > 0) #for The L1 it looks like 0 is invalid
#             dsr['cloud_data'] = dst.astype(_np.float32)
            
#             dsr.cloud_data.attrs = {'long_name':'Cloud base heights.',
#                                             'units':'m',
#                                             'comments':'''A 2D array containing all three cloud bases at each timestep. -999 if no significant signal.''',
#                                             'cloud_base_1':'''First cloud base height or vertical visibility if cloud_status=4''',
#                                             'cloud_base_2':'''Second cloud base height or highest received signal if cloud_status=4''',
#                                             'cloud_base_3':'Third cloud base height'}
            
#             dsr.attrs = {'title':'Ceilometer cloud product',
#                          'version':self.version,
#                          'institution':'NOAA/GML/GRAD',
#                          'author':'hagen.telg@noaa.gov',
#                          'source':'Vaisala CL51 ceilometer',
#                          'serial_number': poutg.sn.unique()[0],
#                          'input_files': ', '.join([fn.name for fn in poutg.path2raw]),
#                          'Conventions':'CF-1.8',
#                          'comments':''' The "time" coordinate was re-indexed to the full minute by back-filling the nearest 
# valid data value within the following minute.The data values have not undergone any processing other than what the Vaisala 
# software has applied. In addition, no QC has been applied other than a visual inspection for impossible values or obvious errors.'''.replace('\n','')
#                          }
            
#             #### attributes of the coordinates
#             dsr.cloud_layer.attrs = dict(long_name = "Cloud layer index",
#                                         units = "1" ,
#                                         axis = "Z" ,
#                                         positive = "up",)

#             dsr.range.attrs = {'long_name': 'distance from ground',
#                                'units': 'm',
#                                'standard_name' : 'distance',
#                                'axis' : 'Z',
#                                'positive': 'up'}
#             #### reindex on 1min data and fill any gaps with nan ... every minute of the day has a timestamp!
#             time = dsr.time.to_pandas()
#             start = time.iloc[0].date()
#             end = start + _pd.to_timedelta(1, 'd')
#             newtime = _pd.date_range(start = start, end = end, freq = _pd.to_timedelta(1, 'minute'), inclusive = 'left')
            
#             assert(time.iloc[0].date() == time.iloc[-1].date()), f'Fist and last timestamp are not on the same date!! {time.iloc[0].date()} vs {time.iloc[-1].date()}'
#             dsr['instrumen_reported_time'] = dsr.time.copy() #this ensures we still have the original timestamp for each value
#             dsr.instrumen_reported_time.attrs = {'long_name': "Instrument reported time",
#                                                  'comments': 'Timestamp that was reported by the instrument. The "time" coordinate was re-indexed to the full minute (see global comments).'}
#             dsr = dsr.reindex(time = newtime, method = 'bfill', tolerance = '1min')
            
            self._product_dataset = ds
                        
            
        return self._product_dataset

class Cl51CloudProdProcessor_v1p3(object):
    def __init__(self, 
                 p2fl_in = '/nfs/grad/Inst/Ceil/SURFRAD/',
                 p2fl_out = '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1_{version}',
                 p2fl_quicklooks = None,
                # file_type = 'bl',
                 # version = '0',
                 hist_file_format = '*_CEILOMETER_1_LEVEL_3*.his',
                  ignore = [],
                 # create_quicklooks = False,
                 reporter = None,
                 verbose = False,
                 ):
        """
        This class aggreates all functions that are needed to create the
        cl51_cloud_prod_lev1 (formally known as cl51_cloud_prod_lev0) product.
        I

        Parameters
        ----------
        p2fl_in : str
            Path to input folder. This folder is expected to have further subfolders, one for each instrument (typically the measurment site, e.g. TBL)
        p2fl_out : str
            Output folder. The default is '/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1'.
        p2fl_quicklooks: str, optional
            Path to quicklooks. If None no quicklooks will be generated. Default is None.
        hist_file_format : TYPE, optional
            Format (pattern) of the hist files. As this format sometimes changes it not hardcoded. The default is '*_CEILOMETER_1_LEVEL_3*.his'.
        file_type: str, optional ['bl', 'hist']
            What filetype to try first. Currently when 'hist' is selected 'bl' is not considered at all (it is unusual that the bl file exists but not the hist file.)
        ignore : list, optional
            When there is a folder in the p2fl_in folder that mathes an element of this list it will be ignored.
        verbose : bool, optional
            If True there will be more output ... mostly for develpment purposes. The default is False.


        Returns
        -------
        Cl51CloudProdProcessor instance.

        """
        self.version = '1.3.1'
        self.reporter = reporter
        self.p2fl_out = _pl.Path(p2fl_out.format(version = self.version))
        self.p2fl_quicklooks = _pl.Path(p2fl_quicklooks.format(version = self.version))
        
        # assert(file_type in ['bl', 'hist']), f'file_type neither hist, nor bl. Is {file_type}'
        # self.file_type = file_type
        self.ignore = ignore
        self.p2fl_in = _pl.Path(p2fl_in)
        self.hist_file_format= hist_file_format
        self.bl_file_format = 'L2*.nc'
        self.fn_format_out = '{site}.cl51.cloud_prod.{date}.nc'
        
        if isinstance(p2fl_quicklooks, type(None)):
            create_quicklooks = False
        else:
            create_quicklooks = True
        self.create_quicklooks = create_quicklooks
        
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
                
                ## bl files
                p2fl_bl = p2site.joinpath('bl')
                df = _pd.DataFrame(p2fl_bl.glob(self.bl_file_format), columns = ['path2raw'])
                
                #### make datetime index
                df.index = df.apply(bl2date, axis = 1)

                ## concat  

                df['site'] = p2site.name.lower()
                workplan = _pd.concat([workplan, df])

            workplan.sort_index(inplace=True)
            
            #### add date collumn
            workplan['date'] = workplan.apply(lambda row: _pd.to_datetime(row.name.date()), axis = 1)
            
            ## create path to output file

            def create_path_out(row):
                fname = self.fn_format_out.format(site = row.site, date = row.name.date().__str__().replace('-',''))
                p2fn = self.p2fl_out.joinpath(row.site).joinpath(f'{row.name.year}').joinpath(fname)
                return p2fn

            workplan['path2fn_out'] = workplan.apply(create_path_out, axis = 1)
            
            if self.create_quicklooks:
                def create_path2plots(row):
                    fname = self.fn_format_out.replace('.nc', '.png').format(site = row.site, date = row.name.date().__str__().replace('-',''))
                    p2fn = self.p2fl_quicklooks.joinpath(row.site).joinpath(f'{row.name.year}').joinpath(fname)
                    return p2fn

                workplan['path2quicklooks'] = workplan.apply(create_path2plots, axis = 1)
                
            
            #### add serial number
            workplan['serial_no'] = workplan.apply(lambda row: row.path2raw.name.split('_')[-1].replace('.nc',''), axis = 1)
            
            #### remove row when output file exists

            workplan = workplan[~ workplan.apply(lambda row: row.path2fn_out.is_file(), axis = 1)]

            #### remove last day
            
            workplan = workplan[~(workplan.date == workplan.date.iloc[-1]).values]
            
            #### add alternative hist files. This is needed if netcdf is corrupt
            workplan['path2althist_l2'] = workplan.apply(lambda row: row.path2raw.parent.parent.joinpath('hist').joinpath(f'{row.name.year}{row.name.month:02d}_CEILOMETER_1_LEVEL_2_{row.name.day:02d}.his'), axis = 1)
            workplan['path2althist_l3'] = workplan.apply(lambda row: row.path2raw.parent.parent.joinpath('hist').joinpath(f'{row.name.year}{row.name.month:02d}_CEILOMETER_1_LEVEL_3_DEFAULT_{row.name.day:02d}.his'), axis = 1)
            
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
                error_handling_missing_level3 = 'raise',
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
        if self.workplan.shape[0]  == 0:
            print('workplan is empty!! DONE!')
            return
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
                retriever = Cl51CloudProdRetriever_v1p3(poutg, self.version, self.reporter)
                self.tp_retriever = retriever
                ####Fixme retriever.check_serial(error_handling=error_handling_serial)
                ds =retriever.product_dataset.dataset
                no_processed += 1
                if not isinstance(self.reporter, type(None)): self.reporter.clean_increment()
                    
                
            except Exception as err:
                #### TODO: there are cases where the level and only the level 3 files are missing, what happens If I still generate a file just with all those values empty
                if isinstance(err, FileNotFoundError):
                    if error_handling_missing_level3 =='return':
                        errors.append(err)
                        errror_grp.append(poutg)
                    else:
                        raise
                    
                    
                else:
                    if error_handling == 'raise':
                        raise
                    
                    #### error handling
                    elif error_handling == 'return':
                        errors.append(err)
                        errror_grp.append(poutg)
                        # errors.append(sys.exc_info())
                        # error, error_txt, trace = sys.exc_info()
                        # tm = ['{}: {}'.format(error.__name__, error_txt.args[0])] + _tb.format_tb(trace)
                if not isinstance(self.reporter, type(None)): self.reporter.errors_increment()
                
            if test == 1:
                break
            
            #### Write to disk.        
            if not isinstance(ds, type(None)):
                assert(not p2fnout.is_file()), f'File exists, this should not be possible! path: {p2fnout}'
            
                if generate_missing_folders:
                    make_all_parent_flds(p2fnout)
                    if self.create_quicklooks:
                        make_all_parent_flds(poutg.path2quicklooks.iloc[0])
                
                #### encodings
                encoding = {}
                for var in ds:
                    encoding[var] = {'zlib': True,     # Enable compression
                                     'complevel': complevel   # Moderate compression level
                                    }
                encoding['time'] = {'units': "seconds since 1900-01-01 00:00:00 UTC",
                                   'zlib': True,     # Enable compression
                                     'complevel': complevel   # Moderate compression level
                                    }
                encoding['instrumen_reported_time'] = {'units': "seconds since 1900-01-01 00:00:00 UTC",
                                                      'zlib': True,     # Enable compression
                                                         'complevel': complevel   # Moderate compression level
                                                        }
                
                #### save to netcdf
                ds.to_netcdf(path=p2fnout, 
                             # mode='w', 
                             encoding = encoding,
                             format='NETCDF4')
                p2fnout.chmod(0o664) # to give the group writing access.
                
                #### plot quicklook
                if self.create_quicklooks:
                    outql = retriever.product_dataset.plot_quicklooks()
                    f = outql[0]
                    f.savefig(poutg.path2quicklooks.iloc[0],bbox_inches = 'tight', dpi = 150)
                    f.clear()
                    _plt.close('all') # this will close all figures, which will otherwise burst the memory
                    _gc.collect()
                print('.', end = '', flush=True)
            else:
                if verbose:
                    print(f'{p2fnout} was not saved as ds is None.')
                print('|', end = '', flush=True)
                
            if test == 2:
                break
            if not isinstance(self.reporter, type(None)): self.reporter.log()
            
            
        
        out['errors'] = errors
        out['errror_grp'] = errror_grp
        out['last_dataset'] = retriever.product_dataset
        out['no_of_files_processed'] = no_processed
        out['finish_time'] = _pd.Timestamp(_dt.datetime.now())
        self._last_processing = out
        if not isinstance(self.reporter, type(None)): self.reporter.log(overwrite_reporting_frequency=True)
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


