import numpy as np
import netCDF4 as nc
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import warnings
import multiprocessing
import ctypes as c
from datetime import datetime, time, timedelta
from scipy.interpolate import interp1d
from optparse import OptionParser
import xarray as xr
from contextlib import contextmanager
from scipy.spatial import cKDTree
from functools import partial
import time
import os

parser = OptionParser()
parser.add_option("-d", "--date", dest="date", type="string", default= None, \
                                  help="Name of netCDF file from 2d run")
parser.add_option("-s", "--start", dest="fstart", type="int", default= None, \
                                  help="Forecast time start")
parser.add_option("-e", "--end", dest="fend", type="int", default=None, \
                                  help="Forecast time end")

(options, args) = parser.parse_args()

if options.date == None or options.fstart == None or options.fend == None:
    print()
    parser.print_help()
    print()
    sys.exit(1)

###### SET FILE NAMES AND CONSTANTS ######
date = options.date #"2020050300"
utc = 18  # Initailization Time
ng = 4  # Number of groups to compare (e.g., WRF & FV3)
ne = 1  # Number of ensemble members in each group of files
fstart = options.fstart #16  # Forecast time start
fend = options.fend #21    # Forecast time end
nf = fend - fstart + 1
nth = 3

group_names = ('HRRR','NAM', 'CTRL', 'NORD3')  # Name of the groups to compare
top_dir = "/work/larissa.reames"
case_dir = "%s/%s"%(top_dir,date)
filter_top_dir = "/work/wicker/CAM_case_studies"
filter_case_dir = "%s/%s"%(filter_top_dir,date)

#filter_file_names = ''
#updraft_file_name = "%s/%s_track_data.nc"%(case_dir,date)
#fcst_type = 0
ens = 0
filter_lvls = [0,12,8]

#/work/wicker/CAM_case_studies/2019071918/ctrl_W_08
#print(updraft_file_name)
z_interp = np.arange(0.0,15000.0,250.0)
nz = z_interp.shape[0]
for g in np.arange(ng):
    updraft_file_name = "%s/%s_%s_track_data.nc"%(case_dir,date,group_names[g])
    fobj = xr.open_dataset(updraft_file_name)
    obj_x = fobj.OBJ_LOCS_X[ens,:,:,:,:] #time,thresh,obj,obj-point
    obj_y = fobj.OBJ_LOCS_Y[ens,:,:,:,:]
    for j,fl in enumerate(filter_lvls):
        print("Begin reading in data for group %s" % group_names[g])
        if (j > 0):
            stream = os.popen("seq -f "+case_dir+"/"+group_names[g].lower()+"_W_%02d"%(fl)+"/filtered_"+date+"00_F%02g.nc"+" %02d %02d"%(fstart,fend))
        else:
            stream = os.popen("seq -f "+case_dir+"/"+group_names[g].lower()+"/region_"+date+"00_F%02g.nc"+" %02d %02d"%(fstart,fend))
        fs = []
        for line in stream.readlines() :
            fs.append(line.strip())
        fs.sort()
        print(fs)
        for n,fname in enumerate(fs):
            print(fname)
            f = xr.open_dataset(fname)
            w = f.W
            z = f.GPH-f.HGT
            nz_cur,_,_ = np.shape(w)
            if (n == 0 and j==0 and g==0): w_objs = np.full((ng,3,nth,nf,100,nz),np.nan)

            for t in np.arange(3):
                obj_x_cur = np.where(obj_x[n,t,...]<0,np.nan, obj_x[n,t,...])
                obj_y_cur = np.where(obj_y[n,t,...]<0,np.nan, obj_y[n,t,...])
                nobj = np.count_nonzero(~np.isnan(obj_x_cur[:,0]))
                for i in np.arange(nobj):
                    obj_size = np.count_nonzero(~np.isnan(obj_x_cur[i,:]))
                    w_interp = np.full(nz,0.0)
                    for p in np.arange(obj_size):
                        y = obj_y_cur[i,p].astype(int)
                        x = obj_x_cur[i,p].astype(int)
                        z_cur = z[:,y,x]
                        w_cur = w[:,y,x]
                        w_interp = w_interp + np.interp(z_interp,z_cur,w_cur )
                    w_objs[g,j,t,n,i,:] = w_interp/obj_size

w_mean = np.nanmean(w_objs,axis=4)

linecolors=("black","xkcd:grey","xkcd:red","xkcd:blue")
linestyles=("-","-","--","--")
linewidths=(2,2,2,2)
#th = 2
#w_mean_all = np.stack((fobj.OBJ_W[:,th,...],w_mean_filter[:,1,th,...],w_mean_filter[:,0,th,...]))
for th in np.arange(nth):
    w_mean_time = np.nanmean(w_mean[:,:,th,:,:],axis=2) #group,filter level, z
    print(w_mean_time.shape)

    #z_interp = fobj.height
    fig,ax = plt.subplots(1,4,sharey='row')
    fig.set_size_inches(14,9)
    #fig.subplots = plt.subplots_adjust(hspace=0.4, wspace=0.4)
    #fig.suptitle(title, fontsize=24)

    for i in np.arange(ng):
        print(i)
        ax[0].plot(w_mean_time[i,0,:],z_interp,color=linecolors[i],linewidth=linewidths[i],linestyle=linestyles[i],label='{}'.format(group_names[i]))
        #ax[0].plot(np.nanmean(fobj.OBJ_W[i,th,...],axis=0),z_interp,color=linecolors[i],linewidth=linewidths[i],linestyle=linestyles[i],label='{}'.format(group_names[i]))

    ax[0].set_xlim((-2.0,35.0))
    ax[0].set_ylim((0,15000))
    ax[0].set_ylabel('Height AGL (m)')
    #ax[0,0].set_xticks(np.arange(0,60,step=12))
    ax[0].set_xlabel('w orig (m $s^{-1}$)')
    ax[0].legend()

    for i in np.arange(0,ng):
        ax[1].plot(w_mean_time[i,1,:],z_interp,color=linecolors[i],linewidth=linewidths[i],linestyle=linestyles[i])

    ax[1].set_xlim((-2.0,35.0))
    ax[1].set_ylim((0,15000))
    #ax[0,0].set_xticks(np.arange(0,60,step=12))
    ax[1].set_xlabel('w 12x filter (m $s^{-1}$)')

    for i in np.arange(0,ng):
        ax[2].plot(w_mean_time[i,2,:],z_interp,color=linecolors[i],linewidth=linewidths[i],linestyle=linestyles[i])

    ax[2].set_xlim((-2.0,35.0))
    ax[2].set_ylim((0,15000))
    #ax[0,0].set_xticks(np.arange(0,60,step=12))
    ax[2].set_xlabel('w 8x filter (m $s^{-1}$)')

    print(z_interp)
    for i in np.arange(0,ng):
        ax[3].plot(w_mean_time[i,2,2:55]/w_mean_time[i,0,2:55],z_interp[2:55],color=linecolors[i],linewidth=linewidths[i],linestyle=linestyles[i])

    ax[3].set_xlim((0.1,0.75))
    ax[3].set_ylim((0,15000))
    #ax[0,0].set_xticks(np.arange(0,60,step=12))
    ax[3].set_xlabel('w 8x filter % change')

    fig.savefig('%s/%s__%3.1f_Wprofiles_filter.png'%(case_dir,date,fobj.threshes[th]))
