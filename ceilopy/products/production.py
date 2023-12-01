#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:09:05 2023

@author: htelg
"""


import ceilopy.products.cl52_cloud_prod_lev1_v1p3 as cl52l1v1p3

def surfrad():
    p2fl_in = '/nfs/grad/Inst/Ceil/SURFRAD/' #'/nfs/grad/gradobs/raw/short_term/sailsplash/Ceil/'
    p2fl_out = '/nfs/grad/surfrad/ceilometer/cl51_cloud_prod_lev1/v{version}'# '/nfs/grad/gradobs/scaled/short_term/sailsplash/ceilometer/cl51_cloud_prod_lev1_v{version}/'
    p2fl_quicklooks = '/nfs/grad/surfrad/quicklooks/ceilometer/cl51_cloud_prod_lev1/v{version}'
    cpp = cl52l1v1p3.Cl51CloudProdProcessor_v1p3(ignore=['plots'], 
                                    # version = version,
                                    # verbose=True,
                                    p2fl_in=p2fl_in,
                                    p2fl_out=p2fl_out,
                                    p2fl_quicklooks=p2fl_quicklooks
                                    )
    
    cpp.workplan = cpp.workplan.truncate('2023-11-01')
    # out = 
    cpp.process(generate_missing_folders=True, 
                      error_handling='return',
                      # test=1
                     )

def run():
    surfrad()
