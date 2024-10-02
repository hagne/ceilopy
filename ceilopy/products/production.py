#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:09:05 2023

@author: htelg
"""


import ceilopy.products.cl51_cloud_prod_lev1_v1p3 as cl51l1v1p3
import subprocess
import productomator.lab as prolab
import pandas as pd


def surfrad(start = None, end = None, lastdays = 14, reporter_name = 'Cl51CloudProd_surfrad'):
    """
    Takes care of the cl52l1 production for surfrad. For manual execution check 
    products/grad/surfrad/ceilometer/cl51_cloud_prod_lev1_v1p3.ipynb

    Returns
    -------
    None.

    """
    print('catchup')
    # print(pd.Timestamp.now())
    # print(f'Starting {reporter_name} production.')
    reporter = prolab.Reporter(reporter_name, 
                               # log_folder='/export/htelg/tmp/', 
                               reporting_frequency=(3,'h'))

    p2fl_in = '/nfs/grad/Inst/Ceil/SURFRAD/' #'/nfs/grad/gradobs/raw/short_term/sailsplash/Ceil/'
    p2fl_out = '/nfs/grad/surfrad/ceilometer/cl51_cloud_prod_lev1/v{version}'# '/nfs/grad/gradobs/scaled/short_term/sailsplash/ceilometer/cl51_cloud_prod_lev1_v{version}/'
    p2fl_quicklooks = '/nfs/grad/surfrad/quicklooks/ceilometer/cl51_cloud_prod_lev1/v{version}'

    cpp = cl51l1v1p3.Cl51CloudProdProcessor_v1p3(ignore=['plots'], 
                                    # version = version,
                                    # verbose=True,
                                    p2fl_in=p2fl_in,
                                    p2fl_out=p2fl_out,
                                    p2fl_quicklooks=p2fl_quicklooks,
                                    reporter = reporter,
                                    )
    #trunc = '2020-01-01'
    if not isinstance(lastdays, type(None)):
        start = pd.Timestamp.now() - pd.to_timedelta(lastdays, 'd')
    
    trunc = (start, end)
    print(f"processing between from {start} to {end}")
    cpp.workplan = cpp.workplan.truncate(*trunc)

    print(f'Number of files to be processed: {cpp.workplan.shape[0]}')
    cpp.workplan = cpp.workplan.sample(frac=1)
    cpp.process(generate_missing_folders=True,
                error_handling='return',
                error_handling_missing_level3='return',
                # test=1
                )
    print('finished processing')
    #### do the rsync
    # Define the source and destination for rsync
    source = cpp.p2fl_out
    destination = f'/nfs/iftp/aftp/g-rad/surfrad/ceilometer/cl51_cloud_prod_lev1/'

    # Construct the rsync command
    rsync_command = ["rsync", "-av", source, destination]
    # Execute the rsync command
    print('starting rsync to ftp', end=' ... ', flush=True)
    try:
        # result = 
        subprocess.run(rsync_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('done')
        # Output the result
        # print("STDOUT:", result.stdout.decode())
        # print("STDERR:", result.stderr.decode())
    except subprocess.CalledProcessError as e:
        # print(f"Rsync command failed with exit code {e.returncode}")
        # print(e.stderr.decode())
        reporter.errors_increment(20)
        # reporter.wrapup()
        raise

    # print(f'clean: {reporter.clean}')
    # print(f'warnings: {reporter.warnings}')
    # print(f'errors: {reporter.errors}')
    # print(pd.Timestamp.now())
    reporter.wrapup()
    print('====================================================')
    return

def wfip3(start = None, end = None, lastdays = 14, reporter_name = 'Cl51CloudProd_wfip3'):
    """
    Takes care of the cl52l1 production for surfrad. For manual execution check 
    products/grad/surfrad/ceilometer/cl51_cloud_prod_lev1_v1p3.ipynb

    Returns
    -------
    None.

    """
    print(f'Starting {reporter_name} production.')
    # print(pd.Timestamp.now())
    # print(f'Starting {reporter_name} production.')
    reporter = prolab.Reporter(reporter_name, 
                               # log_folder='/export/htelg/tmp/', 
                               reporting_frequency=(3,'h'))

    p2fl_in = '/nfs/grad/gradobs/raw/short_term/wfip3/Ceil/' #'/nfs/grad/gradobs/raw/short_term/sailsplash/Ceil/'
    p2fl_out = '/nfs/grad/gradobs/scaled/short_term/wfip3/Ceil/cl51_cloud_prod_lev1/v{version}'# '/nfs/grad/gradobs/scaled/short_term/sailsplash/ceilometer/cl51_cloud_prod_lev1_v{version}/'
    p2fl_quicklooks = '/nfs/grad/gradobs/quicklooks/short_term/wfip3/Ceil/cl51_cloud_prod_lev1/v{version}' #

    cpp = cl51l1v1p3.Cl51CloudProdProcessor_v1p3(ignore=['plots'], 
                                    # version = version,
                                    # verbose=True,
                                    p2fl_in=p2fl_in,
                                    p2fl_out=p2fl_out,
                                    p2fl_quicklooks=p2fl_quicklooks,
                                    reporter = reporter,
                                    )
    #trunc = '2020-01-01'
    if not isinstance(lastdays, type(None)):
        start = pd.Timestamp.now() - pd.to_timedelta(lastdays, 'd')
    
    trunc = (start, end)
    print(f"processing between from {start} to {end}")
    cpp.workplan = cpp.workplan.truncate(*trunc)

    print(f'Number of files to be processed: {cpp.workplan.shape[0]}')
    cpp.workplan = cpp.workplan.sample(frac=1)
    cpp.process(generate_missing_folders=True,
                error_handling='return',
                error_handling_missing_level3='return',
                # test=1
                )
    print('finished processing')
    reporter.wrapup()
    print('====================================================')
    return



def run():
    surfrad()
    wfip3()


if __name__ == '__main__':
    surfrad(start=None, end='2023-12-08', lastdays=None, reporter_name = 'Cl51CloudProd_surfrad_catchup')
    #wfip3(lastdays = None)