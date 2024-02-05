#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import sys
import matplotlib
import matplotlib.pyplot as plt
matplotlib.interactive('t')

def PeriodicBoundaryConditions(lon,lat,data,atmos=True,clon=0.,clat=0.):
    """                                                                                                                                                      
    Close the data set by setting additional values at the end of the arrays                                                                                 
    """
    last = 360.

    # Two possibilities: the mask is an array or it is just equal to False                                                                                   
    if type(data.mask)==np.bool_:
        bc_lon = np.zeros(len(lon)+1)
        bc_lat = np.zeros(lat.size)
        bc_data = np.zeros((data.shape[0],len(lon)+1))
        bc_lon[:-1] = lon[:]-clon
        bc_lon[-1] = last-clon
        bc_lat = lat-clat
        bc_data[:,:-1] = data.data[:,:]
        bc_data[:,-1] = data.data[:,0]

        return bc_lon,bc_lat,np.ma.array(bc_data,mask=False)
    else:
        bc_lon = np.zeros(len(lon)+1)
        bc_data = np.zeros((data.shape[0],len(lon)+1))
        bc_lat = np.zeros(lat.size)
        bc_mask = np.zeros((data.shape[0],len(lon)+1))
        bc_lon[:-1] = lon[:]-clon
        bc_lon[-1] = last-clon
        bc_lat = lat-clat
        bc_data[:,:-1] = data.data[:,:]
        bc_mask[:,:-1] = data.mask[:,:]
        bc_data[:,-1] = data.data[:,0]
        bc_mask[:,-1] = data.mask[:,0]

        return bc_lon,bc_lat,np.ma.array(bc_data,mask=bc_mask.astype(bool))

def Plot(lon,lat,data,i=0,clon=0.,clat=0.,vmin=None,vmax=None,cont=10,fignumber=1,cmap=plt.get_cmap('YlOrRd')):
    """
    Plot using options
    """
    plt.figure(fignumber, figsize=(16, 8))
    plt.clf()
    
    # projection
    projectionfun = ccrs.PlateCarree
    projection = projectionfun(clon)

    ax = plt.axes(projection=projection)
    ax.set_global()
    ax.coastlines()
    ax.set_title("{:.2f}".format(2015 + (i-3780)/12))
    bodr = cartopy.feature.NaturalEarthFeature(category='cultural',
        name='admin_0_boundary_lines_land', scale='50m', facecolor='none', alpha=0.7)
    ax.add_feature(bodr, linestyle='--', edgecolor='k', alpha=1)

    # colorbar
    if vmin is None:
        vmin = np.amin(data)
    if vmax is None:
        vmax = np.amax(data)
    levels = np.linspace(start=vmin, stop=vmax, num=cont)

    # plot
    # load data
    tmp_atmos = data.copy()
    tmp_atmos[tmp_atmos<vmin] = vmin-0.001
    tmp_atmos[tmp_atmos>vmax] = vmax

    # contourf
    im = ax.contourf(lon-clon, lat-clat, tmp_atmos, levels=levels,
                     transform=ccrs.PlateCarree(clon),
                     cmap=cmap)

    # post process image
    cb = plt.colorbar(im)
    cb.ax.set_ylabel('Temperature above threshold')
    
    fname = "mora/mora_{:d}.png".format(i)
    plt.savefig(fname, bbox_inches="tight")
    plt.close("all")

def DS(filename='atminst001700.nc'):
    return nc.Dataset(filename)


def LonLatDat(ds,itimeseries):
    
    KtoC = 273.15
    # # for the temperature used in Mora's maps
    # # get Earth's absolute temperature from 1961 to 1990: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/jgrd.50359 
    # data_abso_mean_1961_1990 = 0.5*(13.7 + 14.0)
    # # get Earth's relative temperature from 1961 to 1990 with base 1951-1980: https://data.giss.nasa.gov/gistemp/
    # data_rel_mean_1961_1990_1951_1980 = 0.01
    # # get Earth's relative temperature from 1880 to 1900 with base 1951-1980: https://data.giss.nasa.gov/gistemp/
    # data_rel_mean_1880_1900_1951_1980 = -0.21
    # # compute Earth's relative temperature from 1961 to 1990 with base 1880-1900
    # rel_mean_1961_1990_1880_1900 = data_rel_mean_1961_1990_1951_1980 - data_rel_mean_1880_1900_1951_1980
    # # compute Earth's absolute temperature from 1880-1900:
    # abso_mean_1880_1900 = data_abso_mean_1961_1990 - rel_mean_1961_1990_1880_1900
    # moy_preind_Mora = abso_mean_1880_1900

    # moy_preind_iloveclim_SSP = 17.35 # Computed from a Business as Usual scenario on the period 1850-1900
    MORAT = np.array([26.5, 27. , 27.5, 28. , 28.5, 29. , 29.5, 30. , 30.5, 31. , 31.5,
                      32. , 32.5, 33. , 33.5, 34. , 34.5, 35. , 35.5, 36. , 36.5, 37. ,
                      37.5, 38. , 38.5, 39. , 39.5, 40. , 40.5, 41. , 41.5, 42. , 42.5,
                      43. , 43.5, 44. , 44.5, 45. , 45.5, 46. , 46.5, 47. , 47.5, 48. ,
                      48.5, 50.0])
    MORAH = np.array([100.109,  89.015,  82.441,  76.078,  71.895,  67.712,  63.529,
                      59.493,  57.152,  54.812,  52.471,  50.131,  47.79 ,  45.449,
                      43.109,  40.807,  39.081,  37.355,  35.628,  33.902,  32.175,
                      30.449,  28.723,  26.996,  25.291,  24.016,  22.741,  21.466,
                      20.191,  18.916,  17.641,  16.366,  15.091,  13.817,  12.55 ,
                      11.614,  10.678,   9.742,   8.805,   7.869,   6.933,   5.997,
                      5.06 ,   4.124,   3.188, 0.000])


    Threshold_T = MORAT
    Threshold_H = MORAH

    MoraDeadlyThreshold = interp1d(Threshold_H,Threshold_T)

    itime = 0

    lon = ds['lon'][:]
    lat = ds['lat'][:]
    t2m = 0.
    rh = 0.
    i = 0
    for itime in itimeseries:
        i=i+1
        #print(i,len(itimeseries))
        t2m += ds['t2m'][:][itime,:,:]-KtoC
        rh += ds['r'][:][itime,:,:]*100.

    dT = np.load('Deltas.npz')['DeltaT']
    drh = np.load('Deltas.npz')['RHFactor']

    t2m = t2m/float(len(itimeseries)) + dT
    rh = rh/float(len(itimeseries))*drh

    rh[rh>100] = 100.
    
    # t2m -= moy_preind_iloveclim_SSP
    data = t2m.data - MoraDeadlyThreshold(rh.data)

    rawdata = np.ma.masked_array(data,mask=t2m.mask)

    lon,lat,rawdata = PeriodicBoundaryConditions(lon,lat,rawdata)

    return lon,lat,rawdata
    
def Series(sy,ny,sd,nd):
    s = []
    for i in range(sy,sy+ny):
        for j in range(sd,sd+nd):
            s.append(j+360*i)
    return s
