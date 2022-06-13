###################################################################################################

import matplotlib
import math
from scipy import *
from scipy import spatial
import pylab as P
import numpy as np
import sys
import os
import netCDF4
from datetime import datetime, timedelta, date, time
from optparse import OptionParser
from news_e_post_cbook import *


def convert(t,indir,outdir,suffix=''):

    #print("TIME : %2d%2d \n" % (fhr,fmin))
    print("suffix = %s" % suffix)
    print("indir = %s" % indir)
    infile = indir + "/fv3_history"+suffix+".nc"
    gridfile = indir + "/grid_spec.nc"
    coordfile = indir + "/atmos_static.nc"

    ##################### Process FV3HISTORY File ########################

    try:                                                 #open 3d file
      fin = netCDF4.Dataset(infile, "r")
      print("Opening %s \n"%infile)
    except:
      print("%s does not exist! \n"% infile)
      sys.exit(1)

    try:                                                 #open grid file
      fgrid = netCDF4.Dataset(gridfile, "r")
      print("Opening %s \n"% gridfile)
    except:
      print("%s does not exist! \n"% gridfile)
      sys.exit(1)

    try:                                                 #open vertical coordinate file
      fcoord = netCDF4.Dataset(coordfile, "r")
      print("Opening %s \n"% coordfile)
    except:
      print("%s does not exist! \n"% coordfile)
      sys.exit(1)  

    ################# Get dimensions and attributes using first file ####################

    ### Get init and valid times ###

    time_var = fin.variables["time"]                             #Read title time string
    start_date=time_var.getncattr("units")
    init_year = start_date[-19:-15]
    init_mon = start_date[-14:-12]
    init_day = start_date[-11:-9]
    init_hr = "%02d"%(int(start_date[-8:-6])+1)
    init_min = start_date[-5:-3]
    start_date_true = init_year + init_mon + init_day + init_hr + init_min

    init_date = init_year + init_mon + init_day               #YYYYMMDD string for output file
    init_time = init_hr + init_min                            #HHMM string for initialization time

    #dt = t * 5
    #print dt
    valid_time_dt = datetime(int(init_year),int(init_mon), int(init_day), int(init_hr), int(init_min)) + timedelta(minutes=float(time_var[t])-60.0)
    valid_time_tt = valid_time_dt.timetuple()
    valid_time = '{:>02d}'.format(valid_time_tt[3]) + '{:>02d}'.format(valid_time_tt[4])              #HHMM string for valid time

    #print valid_time

    valid_yr = valid_time_tt[0]
    valid_mo = valid_time_tt[1]
    valid_dy = valid_time_tt[2]
    valid_hr = valid_time_tt[3]
    valid_min = valid_time_tt[4]
    valid_sc = valid_time_tt[5]

    valid_time_out = "%04d-%02d-%02d_%02d:%02d:%02d" %(valid_yr,valid_mo,valid_dy,valid_hr,valid_min,valid_sc)

    ### Set output path ###
    dt=t
    timestep = str(int(dt))
    if (len(timestep) == 1):
       timestep = '0' + timestep

    outname = "/woffv3_d01_" + valid_time_out + ".nc"         #output file
    output_path = outdir + outname

    ### Get grid/projection info ###

    dx = 3000.0 #fdyn.getncattr("dx")                          #east-west grid spacing (m)
    dy = 3000.0 #fdyn.getncattr("dy")                          #north-south grid spacing (m)
    cen_lat = 0.0 #fdyn.getncattr("cen_lat")              #center of domain latitude (dec deg)
    cen_lon = 0.0 #fdyn.getncattr("cen_lon")              #center of domain longitude (dec deg)
    stand_lon = 0.0 #fdyn.getncattr("cen_lon")             #center lon of Lambert conformal projection
    true_lat_1 = 0.0 #fdyn.getncattr("stdlat1")           #true lat value 1 for Lambert conformal conversion (dec deg)
    true_lat_2 = 0.0 #fdyn.getncattr("stdlat2")           #true lat value 2 for Lambert conformal conversion (dec deg)

    # Constants for PDHS computation
    po = 1.0E5
    To = 300.0
    A = 50.0
    Rd = 287.0

    xlat = fgrid.variables["grid_latt"][:,:]        #latitude (dec deg; Lambert conformal)
    xlon = fgrid.variables["grid_lont"][:,:]        #longitude (dec deg; Lambert conformal)
    hgt = np.squeeze(fin.variables["hgtsfc"][0,:,:])          #terrain height above MSL (m)
    #print "hgt shape = \n" 
    #print np.shape(hgt)
    pdhs = po*np.exp((-To/A)+np.sqrt((To/A)**2-2*hgt*9.81/(A*Rd)))  #dry hydrostatic surface pressure as computed in WRF

    ny = hgt.shape[0]
    nx = hgt.shape[1]
    #ak = np.flip(fcoord.variables["hyam"],0)/1.0E-5
    #bk = np.flip(fcoord.variables["hybm"],0)
    ak = np.flip(fcoord.variables["hyam"],0)/1.0E-5
    bk = np.flip(fcoord.variables["hybm"],0)
    nz = np.size(ak)
    ph = np.zeros((nz,ny,nx))
    p = np.zeros((nz,ny,nx))
    ### Calculate initial and valid time in seconds ###

    init_time_seconds = int(init_hr) * 3600. + int(init_min) * 60.
    valid_time_seconds = int(valid_hr) * 3600. + int(valid_min) * 60.

    if (valid_time_seconds < 43200.):         #Convert values past 0000 UTC, assumes forecast not run past 12 UTC
      valid_time_seconds = valid_time_seconds + 86400.

    ############################### Read WRFOUT variables: ##################################
    #t=0
    #Need to add in 80m winds alongside 100m winds (u100m & v100m) in fv_diagnostics.
    #wspd80= np.hypot(fin.variables["ugrd80m"][t,:,:], fin.variables["vgrd80m"][t,:,:]) 

    w_up = fin.variables["upvvelmax"][t,:,:]

    rain_temp = fin.variables["prateb_ave"][t,:,:]                  #5-min average precip rate (kg/m^2/s)
    rain = rain_temp / 1000.0 * 300.0 * 1000.0                     #convert to mm

    #   soil_moisture[f,:,:] = fin.variables["SMOIS"][0,0,:,:]

    #hail_temp = fin.variables["HAIL_MAXK1"][t,:,:]  

    #hailcast_temp = fdyn.variables["HAILCAST_DHAIL"][t,:,:]
    #hailcast = hailcast_temp * 25.4                     #convert to inches from mm

    u_unstag = np.flip(fin.variables["ugrd"][t,:,:,:],0)                       #expects var dimensions of (nt, nz, ny, nx) with nt = 1
    u = np.zeros((1,nz,ny,nx+1))
    u[0,:,:,0] = u_unstag[:,:,0]
    u[0,:,:,1] = u_unstag[:,:,0]
    for i in np.arange(2,nx+1):
        u[0,:,:,i] = 2*u_unstag[:,:,i-1]-u[0,:,:,i-1]

    v_unstag = np.flip(fin.variables["vgrd"][t,:,:,:],0)
    v = np.zeros((1,nz,ny+1,nx))
    v[0,:,0,:] = v_unstag[:,0,:]
    v[0,:,1,:] = v_unstag[:,0,:]
    for i in np.arange(2,ny+1):
        v[0,:,i,:] = 2*v_unstag[:,i-1,:]-v[0,:,i-1,:]

    w_unstag = np.flip(fin.variables["dzdt"][t,:,:,:],0)
    w = np.zeros((1,nz+1,ny,nx))
    w[0,0,:,:] = w_unstag[0,:,:]
    w[0,1,:,:] = w_unstag[0,:,:]
    for i in np.arange(2,nz+1):
        w[0,i,:,:] = 2*w_unstag[i-1,:,:]-w[0,i-1,:,:]

    ##########ADD IN TO DEFAULT DIAG_TABLE#######
    rel_vort = np.flip(fin.variables["rel_vort"][t,:,:,:],0)

    tk = np.flip(fin.variables["tmp"][t,:,:,:],0)
    print('Min max tk = %f  %f'%(np.min(tk), np.max(tk)))
    qv = np.flip(fin.variables["spfh"][t,:,:,:],0)
    qc = np.flip(fin.variables["clwmr"][t,:,:,:],0)
    qr = np.flip(fin.variables["rwmr"][t,:,:,:],0)
    qi = np.flip(fin.variables["icmr"][t,:,:,:],0)
    qs = np.flip(fin.variables["snmr"][t,:,:,:],0)
    qg = np.flip(fin.variables["grle"][t,:,:,:],0)
    qh = np.flip(fin.variables["hail"][t,:,:,:],0)
    #qh = np.full(qg.shape,0.0)
    qt = qc + qr + qi + qs + qg + qh
    qv = np.where(qv < 0., 0.000001, qv)                      #force qv to be positive definite

    #Creation of "pb" as in wrf from hybrid coordinates
    #for i in arange(nz):
    #  pb[i,...] = ak[i]+bk[i]*pdhs

    #phyd_fv3 = np.flip(fin.variables["pres_hyd"][t,:,:,:],0)

    ##########ADD IN TO DEFAULT DIAG_TABLE#######
    pfull = np.flip(fin.variables["pres"][t,:,:,:],0)
    print(np.shape(pfull))
    print(np.shape(p))
    #print('phyd fv3')
    #print(np.min(phyd_fv3))
    #print(np.max(phyd_fv3))
    #p = pfull-pb   # pseudo-p' from pseudo pb
    th = (100000.00/pfull) ** (287.0/1004.5) * tk - 300.0

    print('Min max th = %f  %f'%(np.min(th), np.max(th)))
    if (np.max(tk) > 350.0) :
      sys.exit()


    #nz = pb.shape[0]
    #phb = np.zeros((nz+1,ny,nx))
    delz_pz = np.zeros((nz+1,ny,nx))
    delz_pz[0,:,:] = hgt
    #if (isnew == 'y') :
    delz_pz[1:,:,:] = -1. * np.flip(fin.variables["delz"][t,:,:,:],0)
    delz = -1. * np.flip(fin.variables["delz"][t,:,:,:],0)
    #else:
    #  delz_pz[1,:,:] = np.flip(fdyn.variables["delz"][t,:,:,:],0)
    #  delz = np.flip(fdyn.variables["delz"][t,:,:,:],0)

    phb = np.cumsum(delz_pz[:,:,:], axis=0)*9.81
    #print(phb[:,100,100])
    ph = np.zeros((nz+1,ny,nx))
    z, dz = calc_height(ph,phb)
    t_v = calc_thv(tk, qv, qt)
    print(np.min(t_v))
    print(np.max(t_v))
    print(np.min(dz))
    print(np.max(dz))

    #phb = np.cumsum(delz[:,:,:], axis=0) * 9.81
    #ph = np.zeros((nz,ny,nx))
    #print('pfull')
    #print(pfull[nz-1,177,267])
    #print(np.max(pfull[nz-1,...]))
    #phyd = calc_phyd(np.squeeze(t_v),np.squeeze(z),pfull[nz-1,:,:])
    #print('phyd calc')
    #print(phyd[nz-1,177,267])
    #print(np.max(phyd[nz-1,...]))
    dbz = np.flip(fin.variables["refl_10cm"][t,:,:,:],0)
    #wz = fin.variables["REL_VORT"][t,:,:,:],0)

    u10 = fin.variables["ugrd10m"][t,:,:]
    v10 =fin.variables["vgrd10m"][t,:,:]

    q2 = fin.variables["spfh2m"][t,:,:]
    t2 = fin.variables["tmp2m"][t,:,:]
    psfc = fin.variables["pressfc"][t,:,:]
    th2 = (100000.00/psfc) ** (287/1004.5) * t2

    #uh25 = fin.variables["uhmax25"][t,:,:]
    #print(np.max(uh25))
    #uh02 = fin.variables["uhmax03"][t,:,:]
    #wz02 = fin.variables["maxvort02"][t,:,:]
    uh02 = calc_uh(w_unstag,rel_vort,np.cumsum(delz[:,:,:], axis=0),delz,0.,2000.)
    uh25 = calc_uh(w_unstag,rel_vort,np.cumsum(delz[:,:,:], axis=0),delz,2000.,5000.)
    #print(np.max(rel_vort))
    #print(np.max(uh25))
    wz02 = calc_avg_vort(rel_vort,np.cumsum(delz[:,:,:], axis=0),delz,0.,2000.)


    swd = fin.variables["dswrf_ave"][t,:,:]
    iwp = np.zeros((ny,nx))
    lwp = np.zeros((ny,nx))

    ### Close WRFOUT file ###

    fin.close()
    del fin

    fgrid.close()
    del fgrid

    fcoord.close()
    del fcoord

    ### Create file and dimensions: ###

    try:
       fout = netCDF4.Dataset(output_path, "w")
    except:
       print( "Could not create %s!\n" % output_path)

    fout.createDimension('time', None)
    fout.createDimension('west_east', nx)
    fout.createDimension('south_north', ny)
    fout.createDimension('bottom_top', nz)
    fout.createDimension('west_east_stag', nx+1)
    fout.createDimension('south_north_stag', ny+1)
    fout.createDimension('bottom_top_stag', nz+1)

    ### Set Attributes: ###

    setattr(fout,'DX',dx)
    setattr(fout,'DY',dy)
    setattr(fout,'CEN_LAT',cen_lat)
    setattr(fout,'CEN_LON',cen_lon)
    setattr(fout,'STAND_LON',stand_lon)
    setattr(fout,'TRUELAT1',true_lat_1)
    setattr(fout,'TRUELAT2',true_lat_2)
    setattr(fout,'PROJECTION','Lambert Conformal')
    setattr(fout,'START_DATE',start_date_true)
    setattr(fout,'INIT_TIME_SECONDS',init_time_seconds)
    setattr(fout,'VALID_TIME_SECONDS',valid_time_seconds)
    setattr(fout,'FORECAST_TIME_STEP',dt)
    setattr(fout,'MAP_PROJ_CHAR',1)
    setattr(fout,'I_PARENT_START',1)
    setattr(fout,'J_PARENT_START',1)
    setattr(fout,'SOUTH-NORTH_PATCH_END_UNSTAG',ny)
    setattr(fout,'WEST-EAST_PATCH_END_UNSTAG',nx)
    setattr(fout,'BOTTOM-TOP_PATCH_END_UNSTAG',nz)
    setattr(fout,'SOUTH-NORTH_GRID_DIMENSION',nx+1)
    setattr(fout,'WEST-EAST_GRID_DIMENSION',ny+1)

    print( 'Create variables')

    xlat1 = fout.createVariable('XLAT', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    xlat1.long_name = "Latitude"
    xlat1.units = "degrees North"

    xlon1 = fout.createVariable('XLON', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    xlon1.long_name = "Longitude"
    xlon1.units = "degrees West"

    hgt1 = fout.createVariable('HGT', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    hgt1.long_name = "Height AGL"
    hgt1.units = "m"

    ppsfc = fout.createVariable('PSFC', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    ppsfc.long_name = "Surface Pressure"
    ppsfc.units = "Pa"

    uu10 = fout.createVariable('U10', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    uu10.long_name = "10-m Zonal Wind"
    uu10.units = "m/s"

    vv10 = fout.createVariable('V10', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    vv10.long_name = "10-m Meridional Wind"
    vv10.units = "m/s"

    tt2 = fout.createVariable('T2', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    tt2.long_name = "2-m Temperature"
    tt2.units = "K"

    tth2 = fout.createVariable('TH2', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    tth2.long_name = "2-m Potential Temperature"
    tth2.units = "K"

    qq2 = fout.createVariable('Q2', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    qq2.long_name = "2-m Water Vapor Mixing Ratio"
    qq2.units = "kg/kg"

    wwspd80 = fout.createVariable('WSPD80', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    wwspd80.long_name = "Maximum 80-m Wind Speed"
    wwspd80.units = "m/s"

    wwup = fout.createVariable('W_UP_MAX', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    wwup.long_name = "Maximum Upwards Vertical Wind Speed"
    wwup.units = "m/s"

    uuh02 = fout.createVariable('UH_02_MAX', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    uuh02.long_name = "Maximum 0-3 km Updraft Helicity"
    uuh02.units = "m^2 s^-2"

    uuh25 = fout.createVariable('UH_25_MAX', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    uuh25.long_name = "Maximum 2-5 km Updraft Helicity"
    uuh25.units = "m^2 s^-2"

    wwz02 = fout.createVariable('WZ_02_MAX', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    wwz02.long_name = "Maximum 0-2 km Vertical Vorticity"
    wwz02.units = "s^-1"

    pprcp = fout.createVariable('PREC_ACC_NC', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    pprcp.long_name = "5-minutes Accumulated Precipitation"
    pprcp.units = "mm"

    hhailmax = fout.createVariable('HAIL_MAXK1', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    hhailmax.long_name = "Max Hail Diameter K=1 (not used)"
    hhailmax.units = "m"

    hhailcastmax = fout.createVariable('HAILCAST_DIAM_MAX', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    hhailcastmax.long_name = "Hailcast Max Hail Diameter"
    hhailcastmax.units = "mm"

    sswdown = fout.createVariable('SWDOWN', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    sswdown.long_name = "Downward Short Wave Flux at Ground Surface"
    sswdown.units = "W/m^2"

    iiwpath = fout.createVariable('IWP', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    iiwpath.long_name = "Ice Water Path (not used)"
    iiwpath.units = "g/m^2"

    llwpath = fout.createVariable('LWP', 'f4', ('time','south_north','west_east'),fill_value=-1e20)
    llwpath.long_name = "Liquid Water Path (not used)"
    llwpath.units = "g/m^2"

    rradar = fout.createVariable('REFL_10CM', 'f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    rradar.long_name = "Simulated Radar Reflectivity"
    rradar.units = "dBz"

    uu = fout.createVariable('U', 'f4', ('time','bottom_top','south_north','west_east_stag'),fill_value=-1e20)
    uu.long_name = "Zonal Wind Speed"
    uu.units = "m/s"

    vv = fout.createVariable('V', 'f4', ('time','bottom_top','south_north_stag','west_east'),fill_value=-1e20)
    vv.long_name = "Meridional Wind Speed"
    vv.units = "m/s"

    ww = fout.createVariable('W', 'f4', ('time','bottom_top_stag','south_north','west_east'),fill_value=-1e20)
    ww.long_name = "Vertical Velocity"
    ww.units = "m/s"

    tt = fout.createVariable('T', 'f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    tt.long_name = "Perturbation Potential Temperature"
    tt.units = "K"

    pp = fout.createVariable('P','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    pp.long_name = "Full Pressure" 
    pp.units = "Pa"

    ppb = fout.createVariable('PB','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    ppb.long_name = "WRF-like Base State Pressure" 
    ppb.units = "Pa"

    ddz = fout.createVariable('DELZ','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    ddz.long_name = "Centered Level Thickness" 
    ddz.units = "m"

    pphb = fout.createVariable('PHB','f4', ('time','bottom_top_stag','south_north','west_east'),fill_value=-1e20)
    pphb.long_name = "Base State Potential Height" 
    pphb.units = "gpm"

    pph = fout.createVariable('PH','f4', ('time','bottom_top_stag','south_north','west_east'),fill_value=-1e20)
    pph.long_name = "Pertubation Potential Height (=0)" 
    pph.units = "gpm"

    rrv = fout.createVariable('REL_VORT','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    rrv.log_name = "Relative Vorticity" 
    rrv.units = "1/s"

    qqv = fout.createVariable('QVAPOR','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    qqv.log_name = "Water Vapor Mixing Ratio" 
    qqv.units = "kg/m^3"

    qqc = fout.createVariable('QCLOUD','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    qqc.log_name = "Cloud Water Mixing Ratio" 
    qqc.units = "kg/m^3"

    qqr = fout.createVariable('QRAIN','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    qqr.log_name = "Rain Water Mixing Ratio" 
    qqr.units = "kg/m^3"

    qqi = fout.createVariable('QICE','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    qqi.log_name = "Ice Mixing Ratio" 
    qqi.units = "kg/m^3"

    qqs = fout.createVariable('QSNOW','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    qqs.log_name = "Snow Mixing Ratio" 
    qqs.units = "kg/m^3"

    qqg = fout.createVariable('QGRAUP','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    qqg.log_name = "Graupel Mixing Ratio" 
    qqg.units = "kg/m^3"

    qqh = fout.createVariable('QHAIL','f4', ('time','bottom_top','south_north','west_east'),fill_value=-1e20)
    qqh.log_name = "Hail Mixing Ratio"
    qqh.units = "kg/m^3"

    print( 'Write variables')

    fout.variables['XLAT'][0,:] = xlat
    fout.variables['XLON'][0,:]= xlon
    fout.variables['HGT'][0,:]= hgt
    fout.variables['PSFC'][0,:] = psfc
    fout.variables['U10'][0,:] = u10
    fout.variables['V10'][0,:] = v10
    fout.variables['T2'][0,:] = t2
    fout.variables['TH2'][0,:] = th2
    fout.variables['Q2'][0,:] = q2
    #fout.variables['WSPD80'][0,:] = wspd80
    fout.variables['W_UP_MAX'][0,:]= w_up
    fout.variables['UH_25_MAX'][0,:]= uh25
    fout.variables['UH_02_MAX'][0,:]= uh02
    fout.variables['WZ_02_MAX'][0,:]= wz02
    fout.variables['REFL_10CM'][0,:] = dbz
    fout.variables['PREC_ACC_NC'][0,:] = rain
    fout.variables['HAIL_MAXK1'][0,:] = hhailmax._FillValue
    fout.variables['HAILCAST_DIAM_MAX'][0,:] = hhailcastmax._FillValue
    fout.variables['SWDOWN'][0,:] = swd
    fout.variables['IWP'][0,:] = iwp
    fout.variables['LWP'][0,:] = lwp

    fout.variables['U'][0,:] = u
    fout.variables['V'][0,:] = v
    fout.variables['W'][0,:] = w
    fout.variables['PB'][0,:] = pfull
    fout.variables['P'][0,:] = p
    fout.variables['DELZ'][0,:] = delz
    fout.variables['PH'][0,:] = ph
    fout.variables['PHB'][0,:] = phb
    fout.variables['T'][0,:] = th
    fout.variables['REL_VORT'][0,:] = rel_vort
    fout.variables['QVAPOR'][0,:] = qv
    fout.variables['QCLOUD'][0,:] = qc
    fout.variables['QRAIN'][0,:] = qr
    fout.variables['QICE'][0,:]= qi
    fout.variables['QSNOW'][0,:] = qs
    fout.variables['QGRAUP'][0,:] = qg
    fout.variables['QHAIL'][0,:] = qh

    ### Close output file ###

    fout.close()
    del fout
