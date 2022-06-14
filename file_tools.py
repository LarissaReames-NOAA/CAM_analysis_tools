import numpy as np
import netCDF4 as ncdf
import matplotlib.pyplot as plt
import xarray as xr
import glob as glob
import os as os
import sys as sys
import pygrib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from gridrad import *
import xesmf as xe
from cmpref import cmpref_mod as cmpref
from functools import partial
import time
import multiprocessing
from contextlib import contextmanager
from zipfile import *
import io
import dask.array as da

import warnings
warnings.filterwarnings("ignore")

# def fxn():
#     warnings.warn("deprecated", DeprecationWarning)

# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     fxn()
    
_nthreads = 2

_Rgas       = 287.04
_gravity    = 9.806

# These are 45 vertical levels that the FV3 puts out - use them here to map ARW to that grid for comparison

plevels = np.asarray([100000.,  97500.,  95000.,  92500.,  90000.,  87500.,  85000.,  82500.,
                       80000.,  77500.,  75000.,  72500.,  70000.,  67500.,  65000.,  62500.,
                       60000.,  57500.,  55000.,  52500.,  50000.,  47500.,  45000.,  42500.,
                       40000.,  37500.,  35000.,  32500.,  30000.,  27500.,  25000.,  22500.,
                       20000.,  17500.,  15000.,  12500.,  10000.,   7000.,   5000.,   3000.,
                        2000.,   1000.,    700.,    500.,    200.])

nz_new = plevels.shape[0]

#--------------------------------------------------------------------------------------------------
# Special thanks to Scott Ellis of DOE for sharing codes for reading grib2

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()
'''
def gridrad_loop(i, files, input_grid, target_grid, regrid):
    file = files[i]
    print(file)
    f = read_file(file)
    f = remove_clutter(f)               
    ds_conus, _ = gridrad_2_xarray(f)

    cref_in = np.asarray(ds_conus['CREF']).transpose()
    ref_in = np.asarray(ds_conus['REF']).transpose()
    hgt = np.asarray(ds_conus['hgt'])

    input_field = esmf.Field(input_grid, name='input grid Composite Reflectivity')
    input_field.data[...] = cref_in

    target_field = esmf.Field(target_grid, name='target grid Composite Reflectivity')
    target_file = regrid(input_field,target_field)
    cref_out = target_field.data[...].transpose()

    input_field = esmf.Field(input_grid, name = 'input grid 3D Reflectivity', ndbounds=[ref_in.shape[2]])
    target_field = esmf.Field(target_grid, name = 'input grid 3D Reflectivity', ndbounds=[ref_in.shape[2]])

    input_field.data[...] = ref_in

    target_field = regrid(input_field,target_field)

    ref_out = target_field.data[...].transpose()

    ref_out = np.where(ref_out < 0.01, np.nan, ref_out)
    cref_out = np.where(cref_out < 0.01, np.nan, cref_out)

    path, base = os.path.split(file)
    path = '/'.join(path.split("/")[:-1])
    dt = base.split('_')[-1]
    date     = dt.split('T')[0]
    time     = dt.split('T')[-1]

    outfilename = "%s/%s_%s_gridrad_econus.nc"%(path,date,time[0:4])


    new = xr.DataArray( ref_out, dims=['nz','ny','nx'], 
                     coords={"height": (["nz"], hgt),
                             "y" : (["ny"], np.arange(ny).astype(float)),
                             "x" : (["nx"], np.arange(nx).astype(float))} )  
    ds_conus = new.to_dataset(name = 'REFL_10CM')

    new = xr.DataArray( cref_out, dims=['ny','nx'], 
                     coords={"y" : (["ny"], np.arange(ny).astype(float)),
                             "x" : (["nx"], np.arange(nx).astype(float))} ) 
    ds_conus['CREF'] = new

    new = xr.DataArray( lat_out.transpose(), dims=['ny','nx'], 
                     coords={"y" : (["ny"], np.arange(ny).astype(float)),
                             "x" : (["nx"], np.arange(nx).astype(float))} ) 
    ds_conus['XLAT'] = new

    new = xr.DataArray( lon_out.transpose(), dims=['ny','nx'], 
                     coords={"y" : (["ny"], np.arange(ny).astype(float)),
                             "x" : (["nx"], np.arange(nx).astype(float))} ) 
    ds_conus['XLONG'] = new

     # Add attributes

    ds_conus.attrs['date']       = date
    ds_conus.attrs['time']       = time
    ds_conus.attrs['gridType']   = 'econus'
    ds_conus.attrs['DateTime']   = datetime.now().strftime("%Y%m%d_%H:%M:%S")
    ds_conus.attrs['TimeStamp']  = datetime.timestamp(datetime.now())

    ds_conus.to_netcdf(outfilename, mode='w') 

def regrid_gridrad(date,fins,ftarget):
    print(fins)
    
    f = read_file(fins[0])           
    lat_in_2d = np.asarray(f['y']['values'])
    lon_in_2d = np.asarray(f['x']['values'])
                  
    lat_in = np.tile(lat_in_2d, (np.shape(lon_in_2d)[0],1))
    lon_in = np.tile(lon_in_2d, (np.shape(lat_in_2d)[0],1)).transpose()

    temp = np.full(lat_in.shape,np.nan)
    
    grid = xr.open_dataset(ftarget,engine='netcdf4')
    lat_out = np.asarray(grid.variables['lats']).transpose()
    lon_out = np.asarray(grid.variables['lons']).transpose()
    nx,ny = lat_out.shape
    
    _ , input_grid, target_grid, regrid = esmf_regrid(lat_in,lon_in,
                                                            lat_out,lon_out,temp,
                                                            return_regrid=True)
    
    with poolcontext(processes=4) as pool:
        t0 = time.time()
        results = pool.map(partial(gridrad_loop, files = fins, input_grid=input_grid, target_grid=target_grid, regrid=regrid ), np.arange(len(fins)))
        t1 = time.time()
'''      
def regrid_mrms(t,date,hours,fins,ftarget,outdir,weights_file,mrms_vars):
    
        
    #Extract the correct zip files
    file = fins[t]
    tmpdir = "%s/tmp"%outdir
    cmd = "mkdir -p %s"%tmpdir
    with ZipFile(file) as z:
        zip_vars = z.namelist()
        zip_hours = [s for s in zip_vars if "00.zip" in s]
        for hr in zip_hours:
            if not any('Rotation' in var for var in mrms_vars ):
                outfname = "%s/%s_MRMS_ECONUS.nc"%(outdir,os.path.splitext(hr)[0])
            else:
                outfname = "%s/%s_ROTMRMS_ECONUS.nc"%(outdir,os.path.splitext(hr)[0])

            hour = os.path.splitext(hr)[0][-4:]
            print(hr)
            if not any (hour in h for h in hours): continue
            
            with z.open(hr) as z2:
                z2_filedata = io.BytesIO(z2.read())
                with ZipFile(z2_filedata) as nested_zip:
                    zip_vars = nested_zip.namelist()
                    matching = []
                    
                    for var in mrms_vars:
                        var_file = [s for s in zip_vars if var in s and not any(x in s for x in
                                               ("conusPlus","alaska","guam","hawaii","carib"))][0]
                        grb_file = '%s/%s'%(tmpdir,var_file)
                        
                        if not os.path.exists(grb_file):
                            nested_zip.extract(var_file,path=tmpdir)
                            
                        ncf = '%s.nc'%(os.path.splitext(grb_file)[0])
                        print(ncf)
                        if not os.path.exists(ncf):
                            cmd = "ncl_convert2nc %s"%(grb_file)
                            os.system(cmd)
                            path,file = os.path.split(ncf)
                            cmd = "mv %s %s"%(file,path)
                            os.system(cmd)

                        matching.append(ncf)
                    
                    print("reading input data into dataset")
                    print(matching)
                    isds = False
                    for i,fname in enumerate(matching):
                        if not isds:
                            tmp = xr.open_dataset(fname)
                            ds_in = tmp.rename(name_dict={list(tmp.keys())[0]:mrms_vars[i]})
                            isds = True
                        else:
                            tmp = xr.open_dataset(fname)
                            ds_in = xr.merge([ds_in,tmp.rename(name_dict={list(tmp.keys())[0]:mrms_vars[i]})],join='override')
                                     
                    print(ds_in)                
                    print("Getting output lat/lons")
                    
                    grid = xr.open_dataset(ftarget,engine='netcdf4')
                    lat_out= grid.variables['lats']
                    lon_out = grid.variables['lons']
                    
                    ds_out = xr.Dataset(coords=dict(lats = (["south_north","west_east"], lat_out,{'units':'degrees_north'}),
                                              lons = (["south_north","west_east"], lon_out, {'units':'degrees_east'})))
                    
                    print("creating regriddiner")   
                    regridder = xe.Regridder(ds_in, ds_out, "bilinear",
                                             filename=weights_file, 
                                             reuse_weights=True)
                    print("create dask graph")
                    ds_out = regridder(ds_in)
                    print("apply regridder")
                    ds_out.compute()
                    
                    print("Creating dataset and writing to file")

                    ds_out.attrs['date']       = date
                    ds_out.attrs['time']       = os.path.splitext(hr)[0][-4:]
                    ds_out.attrs['gridType']   = 'HRRR_ECONUS'
                    ds_out.attrs['DateTime']   = datetime.now().strftime("%Y%m%d_%H:%M:%S")
                    ds_out.attrs['TimeStamp']  = datetime.timestamp(datetime.now())

                    ds_out.to_netcdf(outfname, mode='w') 
                    
def regrid_rrfs(f,date,hours,fins,ftarget,weights_file,outdir,grib2_vars=None):    
    if not grib2_vars:
        grib2_vars = {     #  keys                 / No. of Dims /  Type   / bottomLevel / paramCategory / paramNumber
               'TEMP':     [{'shortName':'t','typeOfLevel':'hybrid'}],
               'OMEGA':    [{'shortName':'wz','typeOfLevel':'hybrid'}],
               'U':        [{'shortName':'u','typeOfLevel':'hybrid'}],
               'V':        [{'shortName':'v','typeOfLevel':'hybrid'}],
               'GPH':      [{'shortName':'gh','typeOfLevel':'hybrid'}],
               'UH':       [{'typeOfLevel': 'heightAboveGroundLayer', 'level':2000, 'shortName':'uphl'}],
               'CREF':     [{'typeOfLevel': 'atmosphere', 'shortName':'refc'}],
               'HGT':      [{'shortName':'orog','typeOfLevel':'surface'}],
               'REFL_10CM':[{'shortName':'refd','typeOfLevel':'hybrid'}],
               }
    #Extract the correct zip files
    file = fins[f]
    tmpdir = "%s/tmp"%outdir
    cmd = "mkdir -p %s"%tmpdir
    
    with ZipFile(file) as z:
        zip_files = z.namelist()
        for file in zip_file:
            if not any (hour in file for h in hours): continue
            
            print("reading input data into dataset")
            for k, key in enumerate(grib2_vars):
                if k == 0:
                    ds_in = xr.open_dataset(file,filter_by_keys=grib_vars[key][0])
                else:
                    ds_in = xr.merge(ds_in,xr.open_dataset(file,filter_by_keys=grib_vars[key][0]))
            
            print("Getting output lat/lons")

            grid = xr.open_dataset(ftarget,engine='netcdf4')
            lat_out= grid.variables['lats']
            lon_out = grid.variables['lons']

            ds_out = xr.Dataset(coords=dict(lats = (["south_north","west_east"], lat_out,{'units':'degrees_north'}),
                                      lons = (["south_north","west_east"], lon_out, {'units':'degrees_east'})))

            print("creating regriddiner")   
            regridder = xe.Regridder(ds_in, ds_out, "bilinear",
                                     filename=weights_file, 
                                     reuse_weights=True)
            print("create dask graph")
            ds_out = regridder(ds_in)
            print("apply regridder")
            ds_out.compute()
       
            ds_out = ds_out.reindex(hybrid=list(reversed(ds.hybrid)))

            print("Creating dataset and writing to file")

            ds_out.attrs['date']       = date
            ds_out.attrs['time']       = os.path.splitext(hr)[0][-4:]
            ds_out.attrs['gridType']   = 'HRRR_ECONUS'
            ds_out.attrs['DateTime']   = datetime.now().strftime("%Y%m%d_%H:%M:%S")
            ds_out.attrs['TimeStamp']  = datetime.timestamp(datetime.now())

            ds_out.to_netcdf(outfname, mode='w') 

def esmf_regrid(lat_in, lon_in, lat_out, lon_out, data, data_name='', return_regrid=False):

    
    input_grid = esmf.Grid(max_index = np.asarray(lat_in.shape), 
                          coord_sys=esmf.CoordSys.SPH_DEG,
                          staggerloc=esmf.StaggerLoc.CENTER)
    
    slat = input_grid.get_coords(1)
    slon = input_grid.get_coords(0)
    slat[...] = lat_in
    slon[...] = lon_in


    input_field = esmf.Field(input_grid, name='input grid %s'%data_name)
    input_field.data[...] = data
    
    target_grid = esmf.Grid(max_index = np.asarray(lat_out.shape), 
                          coord_sys=esmf.CoordSys.SPH_DEG,
                          staggerloc=esmf.StaggerLoc.CENTER)
    
    tlat = target_grid.get_coords(1)
    tlon = target_grid.get_coords(0)
    tlat[...] = lat_out
    tlon[...] = lon_out
    
    target_field = esmf.Field(target_grid, name='target grid %s'%data_name)
    
    regrid = esmf.Regrid(input_field, target_field, regrid_method=esmf.RegridMethod.BILINEAR,
                     unmapped_action=esmf.UnmappedAction.IGNORE)
    target_field = regrid(input_field,target_field)
        
    input_field.destroy()
    
    if return_regrid :
        return target_field.data[...], input_grid, target_grid, regrid
    else :
        return target_field[...]

def grbFile_attr(grb_file):

    dataloc = np.array(grb_file[1].latlons())

    return np.float32(dataloc[0]), np.float32(dataloc[1])

def grbVar_to_slice(grb_obj, type=None):

    """Takes a single grb object for a variable returns a 2D plane"""

    return {'data' : np.float32(grb_obj[0].values), 'units' : grb_obj[0]['units'],
            'date' : grb_obj[0].date, 'fcstStart' : grb_obj[0].time, 'fcstTime' : grb_obj[0].step}

def grbVar_to_cube(grb_obj, type='isobaricInhPa'):

    """Takes a single grb object for a variable containing multiple
    levels. Can sort on type. Compiles to a cube"""

    all_levels = np.array([grb_element['level'] for grb_element in grb_obj])
    types      = np.array([grb_element['typeOfLevel'] for grb_element in grb_obj])

    if type != None:
        levels = []
        for n, its_type in enumerate(types):
            if type == types[n]:
                levels.append(all_levels[n])
        levels = np.asarray(levels)
    else:
        levels = all_levels

    n_levels   = len(levels)
    
    indexes    = np.argsort(levels)[::-1] # highest pressure first
    cube       = np.zeros([n_levels, grb_obj[0].values.shape[0], grb_obj[1].values.shape[1]])

    for k in range(n_levels):
        cube[k,:,:] = grb_obj[indexes[k]].values
        #print("k %d ind %d Level %d obj_level %d:    Min %4.1f      Max %4.1f"%(k, indexes[k], levels[indexes[k]], grb_obj[indexes[k]].level, np.min(cube[k,:,:]), np.max(cube[k,:,:])))

    return {'data' : np.float32(cube), 'units' : grb_obj[0]['units'], 'levels' : levels[indexes],
            'date' : grb_obj[0].date, 'fcstStart' : grb_obj[0].time, 'fcstTime' : grb_obj[0].step}

#--------------------------------------------------------------------------------------------------
# Plotting strings for filtered field labels
    
def title_string(file, level, label, wmax, wmin, eps=None):
    if eps:
        return ("%s at level=%2.2i with EPS=%5.1f \n %s-max: %3.1f       %s-min: %4.2f" % (file, level, label, eps, label, wmax, label, wmin))
    else:
        return ("%s at level=%2.2i \n %s-Max: %3.1f     %s-Min: %4.2f" % (file, level, label, wmax, label, wmin))

#--------------------------------------------------------------------------------------------------
# Quick W plot to figure out where stuff is

def quickplotgrib(file, klevel= 20, cmap = 'turbo', ax=None, filetype='hrrr', \
                                   newlat=[25,50], 
                                   newlon=[-130,-65]):
    """
        Meant to be a quick look at a horizontal field from a grib file, using W.
        Defaults should be good enough, but its possible for some cases to have to changed these
        to look at the whole grid.
        
        Input:
        ------
        
        file:  name of grib file
        
        Options:
        --------
        
        klevel:    the horizontal level to plot.  Values of 20-30 should capture updraft
        cmap:      whatever you like, use a diverging colormap.  Values are automatically scaled.
        ax:        dont mess with this, unless your ax setup is done right (like below)
        filetype:  either HRRR or RRFS grib, the VV variable is hard coded to pick up what you need.
        newlat:    tuple of two lats for zoom, or to test the values for a region grid. Set to "None" to see whole grid.
        newlon:    tuple of two lons for zoom, or to test the values for a region grid. Set to "None" to see whole grid.
        
        LJW December 2021
    """
    
    # open file

    grb_file = pygrib.open(file)

    # Get lat lons

    lats, lons = grbFile_attr(grb_file)
    
    if filetype == 'hrrr':
        grb_var = grb_file.select(name='Vertical velocity')
        cube = grbVar_to_cube(grb_var, type='hybrid')['data']
    else:
        grb_var = grb_file.select(name='Geometric vertical velocity')
        cube = grbVar_to_cube(grb_var)['data']
        
    glat_min = lats.min()
    glat_max = lats.max()
    glon_min = lons.min()
    glon_max = lons.max()

    print(f'\nGrib File Lat Min: %4.1f  Lat Max:  %4.1f' % (glat_min, glat_max))
    print(f'\nGrib File Lon Min: %4.1f  Lon Max:  %4.1f' % (glon_min, glon_max))
    
    print(f'\nGrib File W Min: %4.1f  W Max:  %4.1f\n' %(-cube[klevel].max(), cube[klevel].min()))
    
    if ax == None:
        
        proj = ccrs.LambertConformal(central_latitude = 30, 
                                     central_longitude = 265., 
                                     standard_parallels = (10,10))
        
        fig = plt.figure(figsize=(20, 20))

        ax = plt.axes(projection = proj)

        ax.set_global()
        ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.STATES, linestyle=':')
        
        if newlat == None:
            lat_min, lat_max = glat_min, glat_max
        else:
            lat_min, lat_max = newlat
            
        if newlon == None:
            lon_min, lon_max = glon_min, glon_max
        else:
            lon_min, lon_max = newlon
            
        print(f'\nPlot Lat Min: %4.1f  Lat Max:  %4.1f' % (lat_min, lat_max))
        print(f'\nPlot Lon Min: %4.1f  Lon Max:  %4.1f' % (lon_min, lon_max))

        ax.set_extent([lon_min, lon_max, lat_min, lat_max])

# # Add variety of features
#         # ax.add_feature(cfeature.LAND)
#         # ax.add_feature(cfeature.OCEAN)
#         # ax.add_feature(cfeature.COASTLINE)

# # Can also supply matplotlib kwargs
#        ax.add_feature(cfeature.LAKES, alpha=0.5)    
# if klevel < 10:
#     vmin = -5.
#     vmax = 10.
#     clevels = np.linspace(vmin, vmax, 16)
# else:
#     vmin = -10.
#     vmax = 20.
#     clevels = np.linspace(vmin, vmax, 16)
                
    title = title_string(os.path.basename(file), klevel, 'W', cube[klevel].max(), cube[klevel].min())
    
    ax.pcolormesh(lons, lats, -cube[klevel], cmap=cmap, transform=ccrs.PlateCarree())
                                    
    ax.set_title(title,fontsize=20)
    
    return 

#--------------------------------------------------------------------------------------------------
# Interp from 3D pressure to 1D pressure (convert from hybrid to constant p-levels)

def interp3d_np(data, p3d, p1d, nthreads = _nthreads):
    
    dinterp = np.zeros((len(p1d),data.shape[1],data.shape[2]),dtype=np.float32)

    if nthreads < 0:  # turning this off for now.
        def worker(i):
            print("running %d %s" % (i, data.shape))
            for j in np.arange(data.shape[1]):
                  dinterp[:,j,i] = np.interp(p1d[::-1], p3d[:,j,i], data[:,j,i])

        pool = mp.Pool(nthreads)
        for i in np.arange(data.shape[2]):
            pool.apply_async(worker, args = (i, ))
        pool.close()
        pool.join()
        
        return dinterp[::-1,:,:]
    
    else:        
        for i in np.arange(data.shape[2]):
            for j in np.arange(data.shape[1]):
                dinterp[:,j,i] = np.interp(p1d[::-1], p3d[:,j,i], data[:,j,i])
        
        return dinterp[::-1,:,:]

#--------------------------------------------------------------------------------------------------
# Choose a section of the grid based on lat/lon corners - excludes the rest of grid from xarray

def extract_subregion(xr_obj, sw_corner=None, ne_corner=None, drop=True):
    
    if (sw_corner and len(sw_corner) > 1) and (ne_corner and len(ne_corner) > 1):
        lat_min = min(sw_corner[0], ne_corner[0])
        lat_max = max(sw_corner[0], ne_corner[0])
        lon_min = min(sw_corner[1], ne_corner[1])
        lon_max = max(sw_corner[1], ne_corner[1])
        
        print(f'Creating a sub-region of grid: {lat_min:.2f}, {lon_min:.2f}, {lat_max:.2f}, {lon_max:5.2f}','\n')
        
        xr_obj.attrs['gridType'] = 'econus'
        lats = xr_obj.lats.values
        lons = xr_obj.lons.values
        print(xr_obj)
        return xr_obj.where(xr.DataArray(np.logical_and(np.logical_and(lats>sw_corner[0],lats<ne_corner[0]), 
                            np.logical_and(lons>sw_corner[1],lons<ne_corner[1])),dims=['south_north','west_east']),drop=True)
    else:
        print(f"No grid information supplied - returning original grid!\n")
        return xr_obj

#--------------------------------------------------------------------------------------------------

def hrrr_grib_read_variable(file, sw_corner=None, ne_corner=None, var_list=[''], interpP=True, writeout=True, prefix=None) :
    
    # Special thanks to Scott Ellis of DOE for sharing codes for reading grib2
    
    default = {            #  Grib2 name                 / No. of Dims /  Type   / bottomLevel / paramCategory / paramNumber
               'TEMP':     ['Temperature',                           3, 'hybrid'],
               'OMEGA':    ['Vertical velocity',                     3, 'hybrid'],
               'U':        ['U component of wind',                   3, 'hybrid'],
               }
 #              'V':        ['V component of wind',                   3, 'hybrid'],
 #              'GPH':      ['Geopotential Height',                   3,  'hybrid'],
 #              'UH':       ['unknown',                               2,  'hybrid', 2000,         7,              199],
 #              'CREF':     ['Maximum/Composite radar reflectivity',  2,  'atmosphere',           16,             196],
 #              'HGT':      ['Orography',                             2,  'surface'],
 #              'REFL_10CM':['Reflectivity',                          3, 'hybrid']
 #              }

    if var_list != ['']:
        variables = {k: var_list[k] for k in var_list.keys() & set(var_list)}  # yea, I stole this....
    else:
        variables = default

    if prefix == None:
        prefix = 'hrrr'

    print(f'-'*120,'\n')
    print(f'HRRR Grib READ: Extracting variables from grib file: {file}','\n')
    
    level_type_grib = 105
    level_type_string = 'hybrid'

    # open file

    grb_file = pygrib.open(file)

    # Get lat lons

    lats, lons = grbFile_attr(grb_file)
    lats = xr.DataArray(lats,dims=['south_north','west_east'],attrs=dict(units='degrees_north'))
    lons = xr.DataArray(lons,dims=['south_north','west_east'],attrs=dict(units='degrees_east'))
    
    if (np.amax(lons) > 180.0): lons = lons - 360.0

    pres = None

    
    grb_var = grb_file.select(name='Pressure', typeOfLevel=level_type_string)
    levels = np.array([grb_element['level'] for grb_element in grb_var])
    indexes    = np.argsort(levels)
    levels = levels[indexes]
    cube = grbVar_to_cube(grb_var, type='hybrid')
    p3d = cube['data']
    print(f'InterpP is True, Read 3D pressure field from GRIB file\n')
    print(f'P-max:  %5.2f  P-min:  %5.2f\n' % (p3d.max(), p3d.min()))

    for n, key in enumerate(variables):

        print('Reading my variable: ',key, 'from GRIB file variable: %s\n' % (variables[key][0]))
        if key == 'REFL_10CM' :
            cmr = grbVar_to_cube(grb_file.select(name='Cloud mixing ratio'), type=variables[key][2])['data']
            rmr = grbVar_to_cube(grb_file.select(name='Rain mixing ratio'), type=variables[key][2])['data']
            smr = grbVar_to_cube(grb_file.select(name='Snow mixing ratio'), type=variables[key][2])['data']
            grl = grbVar_to_cube(grb_file.select(name='Graupel (snow pellets)'), type=variables[key][2])['data']
            qv = grbVar_to_cube(grb_file.select(name='Specific humidity'), type=variables[key][2])['data']
            t = grbVar_to_cube(grb_file.select(name='Temperature'), type=variables[key][2])['data']
            p = grbVar_to_cube(grb_file.select(name='Pressure'), type=variables[key][2])
            nlev = np.shape(p['levels'])
            nz,ny,nx = np.shape(t)
            grb_var = cmpref.calcrefl10cm(qv, cmr, rmr, smr, grl, t, p['data'], nlev, nx, ny)
            cube = {'data' : np.float32(grb_var), 'units' : 'dBZ', 'levels' : p['levels'],
            'date' : p['date'], 'fcstStart' : p['fcstStart'], 'fcstTime' : p['fcstTime']}
        elif type(variables[key][0]) == type('1'):
            if len(variables[key]) == 3:
                grb_var = grb_file.select(name=variables[key][0],typeOfLevel=variables[key][2])
            elif len(variables[key]) == 6:
                grb_var = grb_file.select(bottomLevel=variables[key][3],parameterCategory=variables[key][4],parameterNumber=variables[key][5])
            else :
                grb_var = grb_file.select(parameterCategory=variables[key][3],parameterNumber=variables[key][4])
        else:
            grb_var = [grb_file.message(variables[key][0])]

        if variables[key][1] == 3:
    
            if (key != 'REFL_10CM') : cube = grbVar_to_cube(grb_var, type=variables[key][2])
            
            if interpP:
                   
                cubeI = interp3d_np(cube['data'], p3d, plevels)
                plevels = xr.DataArray(plevels, dims=['bottom_top'], units='Pa')
                    
                new = xr.DataArray( cubeI, dims = ['bottom_top','south_north','west_east'], 
                                   coords=dict(lats = (["south_north","west_east"], lats, {'units':'degrees_north'}),
                                               lons = (["south_north","west_east"], lons, {'units':'degrees_east'}),
                                               pres = (["bottom_top"], plevels)))
                
            else:

                new = xr.DataArray( cube['data'], dims = ['bottom_top','south_north','west_east'],
                                   coords=dict(lats = (["south_north","west_east"], lats,{'units':'degrees_north'}),
                                               lons = (["south_north","west_east"], lons, {'units':'degrees_east'}),
                                               hybrid = (["bottom_top"], cube['levels'])))
        if variables[key][1] == 2:

            cube = grbVar_to_slice(grb_var, type=variables[key][2])

            new = xr.DataArray( cube['data'], dims=['south_north','west_east'], 
                                             coords=dict(lats = (["south_north","west_east"], lats,{'units':'degrees_north'}),
                                              lons = (["south_north","west_east"], lons, {'units':'degrees_east'}))) 

        if n == 0:
            ds_conus = new.to_dataset(name=key)
        else:         

            ds_conus[key] = new

        del(new)
        
        date      = cube['date'] 
        fcstStart = cube['fcstStart']
        fcstHour  = cube['fcstTime']
        
    # Add attributes
    
    ds_conus.attrs['date']       = date
    ds_conus.attrs['fcstStart']  = fcstStart
    ds_conus.attrs['fcstHour']   = fcstHour
    ds_conus.attrs['gridPrefix'] = prefix
    ds_conus.attrs['gridType']   = 'conus'
    ds_conus.attrs['DateTime']   = datetime.now().strftime("%Y%m%d_%H:%M:%S")
    ds_conus.attrs['TimeStamp']  = datetime.timestamp(datetime.now())

    # clean up grib file
    
    grb_file.close()
            
    # Convert omega --> w
    
    if (interpP) :
        pp3 = np.broadcast_to(plevels, (p3d.shape[2],p3d.shape[1],len(plevels))).transpose()
    else:
        pp3 = p3d
        
        
    print(np.amax(ds_conus['OMEGA'].data[0,:,:]))
    print(np.amax(pp3[0,:,:]))
    print(np.amax(ds_conus['TEMP'].data[0,:,:]))
          
    w_new = -ds_conus['OMEGA'].data / ( (_gravity * pp3 ) / (_Rgas * ds_conus['TEMP'].data) )
   #coords=dict(
#...         lon=(["x", "y"], lon),
#...         lat=(["x", "y"], lat),
    if interpP : 
        ds_conus['W'] = xr.DataArray( w_new, dims = ['bottom_top','south_north','west_east'], 
                                      coords=dict(lats = (["south_north","west_east"], lats, {'units':'degrees_north'}),
                                              lons = (["south_north","west_east"], lons, {'units':'degrees_east'}),
                                              pres = (["bottom_top"], plevels)))
        ds_conus['P'] = xr.DataArray( pp3, dims = ['bottom_top','south_north','west_east'], 
                                      coords=dict(lats = (["south_north","west_east"], lats, {'units':'degrees_north'}),
                                              lons = (["south_north","west_east"], lons, {'units':'degrees_east'}),
                                              pres = (["bottom_top"], plevels)))  
    else:

        ds_conus['W'] = xr.DataArray( w_new, dims = ['bottom_top','south_north','west_east'], 
                                  coords=dict(lats = (["south_north","west_east"], lats, {'units':'degrees_north'}),
                                              lons = (["south_north","west_east"], lons, {'units':'degrees_east'}),
                                              hybrid = (["bottom_top"], cube['levels'])))
        ds_conus['P'] = xr.DataArray( pp3, dims = ['bottom_top','south_north','west_east'], 
                                  coords=dict(lats = (["south_north","west_east"], lats, {'units':'degrees_north'}),
                                              lons = (["south_north","west_east"], lons, {'units':'degrees_east'}),
                                              hybrid = (["bottom_top"], cube['levels'])))
    # extract region
    ds_conus = extract_subregion(ds_conus, sw_corner=sw_corner, ne_corner=ne_corner)
    print(ds_conus)
    # add some file and directory attributes
    
    dir, base = os.path.split(file)
    outfilename = os.path.join(dir, '%s_%8.8i%04d_F%2.2i.nc' % (ds_conus.attrs['gridType'], date, fcstStart, fcstHour))
        
    ds_conus.attrs['srcdir']   = dir
    ds_conus.attrs['filename'] = os.path.basename(outfilename)
 
    if writeout:
        os.system("export HDF5_USE_FILE_LOCKING=FALSE")
        ds_conus.to_netcdf(outfilename, mode='w', engine='netcdf4')  
        print(f'Successfully wrote new data to file:: {outfilename}','\n')
        return ds_conus, outfilename 
    
    else:
        
        return ds_conus, outfilename
    
#--------------------------------------------------------------------------------------------------

def fv3_grib_read_variable(file, sw_corner=None, ne_corner=None, var_list=[''], writeout=True, prefix=None):
    
    # Special thanks to Scott Ellis of DOE for sharing codes for reading grib2
    
    default = {             #  Grib2 name                 / No. of Dims /  Type
               'TEMP':     ['Temperature',                 3, 'isobaricInhPa'],
               'HGT':      ['Orography',                   2, 'surface'],
               'GPH':      ['Geopotential Height',         3, 'isobaricInhPa'],              
               'W':        ['Geometric vertical velocity', 3, 'isobaricInhPa'],
               'U':        ['U component of wind',         3, 'isobaricInhPa'],
               'V':        ['V component of wind',         3, 'isobaricInhPa'],
               'UH':       ['unknown',                               2,  'isobaricInhPa', 2000,         7,              199],
               'CREF':     ['Maximum/Composite radar reflectivity',  2, 'atmosphereSingleLayer'      ],
               'REFL_10CM':['Reflectivity',                          3, 'isobaricInhPa']
               }

    if var_list != ['']:
        variables = {k: default[k] for k in default.keys() & set(var_list)}  # yea, I stole this....
    else:
        variables = default

    if prefix == None:
        prefix = 'rrfs'

    print(f'-'*120,'\n')
    print(f'RRFS Grib READ: Extracting variables from grib file: {file}','\n')

    # open file

    grb_file = pygrib.open(file)

    # Get lat lons

    lats, lons = grbFile_attr(grb_file)

    if np.amax(lons) > 180.0: lons = lons - 360.0

    for n, key in enumerate(variables):

        print('Reading my variable: ', key, 'from GRIB variable: %s\n' % (variables[key][0]))
        if key == 'REFL_10CM' :
            cmr = grbVar_to_cube(grb_file.select(name='Cloud mixing ratio'), type=variables[key][2])['data'][::-1,:,:]
            cmr = np.where(cmr == 9999.0, 0.0, cmr)
            
            rmr = grbVar_to_cube(grb_file.select(name='Rain mixing ratio'), type=variables[key][2])['data'][::-1,:,:]
            rmr = np.where(rmr == 9999.0, 0.0, rmr)
            
            smr = grbVar_to_cube(grb_file.select(name='Snow mixing ratio'), type=variables[key][2])['data'][::-1,:,:]
            smr = np.where(smr == 9999.0, 0.0, smr)
            
            grl = grbVar_to_cube(grb_file.select(name='Graupel (snow pellets)'), type=variables[key][2])['data'][::-1,:,:]
            grl = np.where(grl == 9999.0, 0.0, grl)
            
            qv = grbVar_to_cube(grb_file.select(name='Specific humidity'), type=variables[key][2])['data'][::-1,:,:]
            qv = np.where(qv == 9999.0, 1.E-8, qv)
            
            t = grbVar_to_cube(grb_file.select(name='Temperature'), type=variables[key][2])
            tdata = t['data']
            tdata = np.where(tdata == 9999.0, 200.0, tdata)
            
            nz,ny,nx = np.shape(qv)
            p = np.broadcast_to(t['levels'][::-1], (nx,ny,nz)).transpose()
            if (np.nanmax(p) < 100000.0): p = p * 100.0
            grb_var = cmpref.calcrefl10cm(qv, cmr, rmr, smr, grl, tdata[::-1,:,:], p, nz, nx, ny)
            print('after cmpref')
            cube = {'data' : np.float32(grb_var[::-1,:,:]), 'units' : 'dBZ', 'levels' : t['levels'],
            'date' : t['date'], 'fcstStart' : t['fcstStart'], 'fcstTime' : t['fcstTime']}
        elif type(variables[key][0]) == type('1'):
            if len(variables[key]) == 3:
                grb_var = grb_file.select(name=variables[key][0],typeOfLevel=variables[key][2])
            else:
                grb_var = grb_file.select(bottomLevel=variables[key][3],parameterCategory=variables[key][4],parameterNumber=variables[key][5])
        else:
            grb_var = [grb_file.message(variables[key][0])]
        
        if variables[key][1] == 3:
            if key != 'REFL_10CM' : cube = grbVar_to_cube(grb_var, type='isobaricInhPa')
            pres = cube['levels']
            new = xr.DataArray( cube['data'], dims=['bottom_top','south_north','west_east'], 
                               coords={'pres': (['bottom_top'], pres),
                                       "lons": (['south_north','west_east'], lons),
                                       "lats": (['south_north','west_east'], lats)} )      
        if variables[key][1] == 2:
            cube = grbVar_to_slice(grb_var)
            new = xr.DataArray( cube['data'], dims=['south_north','west_east'], 
                               coords={"lons": (['south_north','west_east'], lons),
                                       "lats": (['south_north','west_east'], lats)} )      

        if n == 0:
            
            ds_conus = new.to_dataset(name = key)
            
        else:
            
            ds_conus[key] = new
            
        del(new)
        
        date      = cube['date'] 
        fcstStart = cube['fcstStart']
        fcstHour  = cube['fcstTime']

    # clean up grib file
    
    grb_file.close()
    
    # Useful info for global attributes
         
    ds_conus.attrs['date']       = date
    ds_conus.attrs['fcstStart']  = fcstStart
    ds_conus.attrs['fcstHour']   = fcstHour
    ds_conus.attrs['gridPrefix'] = prefix
    ds_conus.attrs['gridType']   = 'conus'
    ds_conus.attrs['DateTime']   = datetime.now().strftime("%Y%m%d_%H:%M:%S")
    ds_conus.attrs['TimeStamp']  = datetime.timestamp(datetime.now())

    # extract region
    
    ds_conus = extract_subregion(ds_conus, sw_corner=sw_corner, ne_corner=ne_corner)

    # add some file and directory attributes
    
    dir, base   = os.path.split(file)
    outfilename = os.path.join(dir, '%s_%8.8i%04d_F%2.2i.nc' % (ds_conus.attrs['gridType'], date, fcstStart, fcstHour))
        
    ds_conus.attrs['srcdir']   = dir
    ds_conus.attrs['filename'] = os.path.basename(outfilename)
 
    if writeout:
       
        ds_conus.to_netcdf(outfilename, mode='w')  
        print(f'Successfully wrote new data to file:: {outfilename}','\n')
        return ds_conus, outfilename 
    
    else:
        
        return ds_conus, outfilename

#--------------------------------------------------------------------------------------------------


