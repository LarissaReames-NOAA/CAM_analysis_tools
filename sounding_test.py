#!/usr/bin/env python
import glob
import os
import numpy as np
from netCDF4 import Dataset
import sharppy.sharptab.profile as profile
import datetime
from test_wof_sounding_plot import plot_wof
import sharppy.sharptab.utils as utils
import sharppy.sharptab.thermo as thermo
import argparse
import time
import multiprocessing
import sys
import warnings
from analsnd import analsnd_mod

warnings.filterwarnings("ignore",message="converting a masked element to nan",append=True)
warnings.filterwarnings("ignore",message="overflow encountered in multiply",append=True)
warnings.simplefilter(action="ignore",category=RuntimeWarning)
warnings.simplefilter(action="ignore",category=UserWarning)

dz = 75
z_top = 20000.0 #m
hgt = np.arange(0.0,z_top, dz)
psfc = 1000.0
tsfc = 300.5 
qvsfc = 14.0
hodo = 1
thermo_snd = 2
nz = len(hgt)
us1 = 12.5

t,q,pi,p,rh,u,v = analsnd_mod.analsnd(dz,hgt,psfc,tsfc,qvsfc,hodo,thermo_snd,us1=us1,nz=nz,zlfc_in=500.0,tlfc_in=303.5)
q = np.array(q[:nz])/1000.
t = np.array(t[:nz])
p = np.array(p[:nz])/100.
t = t * (p/1000.)**(287./1004.)
u = np.array(u[:nz])*1.94384
v = np.array(v[:nz])*1.94384
omega = np.full(v.shape,0.0)

td = thermo.temp_at_mixrat(q*1000., p)
tc = t - 273.15

print(u)
print(v)

start_date = '201905202100'
idateobj = datetime.datetime.strptime(start_date,'%Y%m%d%H%M')
vdateobj = idateobj 
member_profs = np.empty(1, dtype=object)
prof = profile.create_profile(profile='convective', pres=p, hght=hgt, tmpc=tc, \
                    dwpc=td, u=u, v=v, omega=omega, missing=-9999, latitude=35.0, strictQC=False,date=vdateobj )
member_profs[0] = prof
members = {'hght': hgt, 'pres': p, 'tmpc': t, 'dwpc': td, 'u': u, 'v': v, 'member_profs': member_profs}

plot_wof(prof, members, 'test', 0.0, 0.0, idateobj, vdateobj,x_pts=(0,0),y_pts=(0,0))
