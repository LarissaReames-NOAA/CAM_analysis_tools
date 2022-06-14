import numpy as np
import netCDF4 as ncdf
import matplotlib as mlab
import matplotlib.pyplot as plt
import xarray as xr
import glob as glob
import os as os
import sys as sys
from file_tools import *
from functools import partial
import time
import multiprocessing
from contextlib import contextmanager
from zipfile import *
import io

#ECONUS Domain
sw_corner = (26.0, -105.0)
ne_corner = (49.0, -70.0)

dates = ("202205120000","202205190000","202205270000","202205300000")

hours = np.array([ [ '%02d00'%x for x in np.append(np.arange(21,24),np.arange(0,4))],
                   [ '%02d00'%x for x in np.append(np.arange(18,24),np.arange(0,4))],
                   [ '%02d00'%x for x in np.append(np.arange(13,24),np.arange(0,2))],
                   [ '%02d00'%x for x in np.append(np.arange(18,24),np.arange(0,4))] ])

mrms_refl_vars = ["EchoTop_30","EchoTop_50","EchoTop_60",'VIL',"MergedReflectivityQC_01.00",
             "MergedReflectivityQComposite","HeightCompositeReflectivity","Reflectivity_-10C"]

mrms_rot_vars = ["RotationTrackML60min","RotationTrack60min"]
mrms_vars = mrms_rot_vars

for d,date in enumerate(dates):
    print(date)

    # RRFS2 dir

    #rrfs_dir   = '/Users/Louis.Wicker/CAM_Case_Studies/20210526/RRFS/ctrl'
    #rrfs_files = sorted(glob.glob(rrfs_dir+'/PRSLEV.Grb*'))
    rrfs_dir2 = '/scratch1/BMC/gsd-fv3/Larissa.Reames/%s00'%date
    rrfs_files2 = sorted(glob.glob(rrfs_dir2+'/*.zip'))

    #MRMS_dir
    mrms_dir = '/scratch1/BMC/gsd-fv3/Larissa.Reames/%s'%date
    mrms_files = sorted(glob.glob(mrms_dir+'/*.zip'))

    with poolcontext(processes=6) as pool:
        t0 = time.time()
        results = pool.map(partial(regrid_mrms, date=date, hours=hours[d], fins = mrms_files,
                                   ftarget='/scratch1/BMC/gsd-fv3/Larissa.Reames/econus_gridspec.nc',
                                   weights_file='/scratch1/BMC/gsd-fv3/Larissa.Reames/rotmrms_to_hrrr_econus_weights.nc',
                                  outdir=mrms_dir,mrms_vars=mrms_vars), np.arange(len(mrms_files)))
        t1 = time.time()
