#!/usr/bin/env python3

import numpy as np
import script_new as script
import netCDF4 as nc
import matplotlib.pyplot as plt

ds = nc.Dataset('atminst001700.nc')

#print('Start south 2015')

#its_south_2015 = script.Series(310,10,-15,30)
#lon,lat,data_south_2015 = script.LonLatDat(ds,its_south_2015)
#np.savez('MoraK_10years_averaged_SouthernSummer_2015_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_south_2015)

#print('Start north 2015')

#its_north_2015 = script.Series(310,10,165,30)
#lon,lat,data_north_2015 = script.LonLatDat(ds,its_north_2015)
#np.savez('MoraK_10years_averaged_NorthernSummer_2015_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_north_2015)

data = np.load('MoraK_10years_averaged_SouthernSummer_2015_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'January 2015',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_SouthernSummer_2015_ERA5Corr_TandH.pdf')

data = np.load('MoraK_10years_averaged_NorthernSummer_2015_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'July 2015',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_NorthernSummer_2015_ERA5Corr_TandH.pdf')

print('Start south 2050')

its_south_2050 = script.Series(345,10,-15,30)
lon,lat,data_south_2050 = script.LonLatDat(ds,its_south_2050)
np.savez('MoraK_10years_averaged_SouthernSummer_2050_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_south_2050)

print('Start north 2050')

its_north_2050 = script.Series(345,10,165,30)
lon,lat,data_north_2050 = script.LonLatDat(ds,its_north_2050)
np.savez('MoraK_10years_averaged_NorthernSummer_2050_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_north_2050)

# 2050

data = np.load('MoraK_10years_averaged_SouthernSummer_2050_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'January 2050',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_SouthernSummer_2050_ERA5Corr_TandH.pdf')

data = np.load('MoraK_10years_averaged_NorthernSummer_2050_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'July 2050',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_NorthernSummer_2050_ERA5Corr_TandH.pdf')

print('Start south 2068')

its_south_2068 = script.Series(363,10,-15,30)
lon,lat,data_south_2068 = script.LonLatDat(ds,its_south_2068)
np.savez('MoraK_10years_averaged_SouthernSummer_2068_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_south_2068)

print('Start north 2068')

its_north_2068 = script.Series(363,10,165,30)
lon,lat,data_north_2068 = script.LonLatDat(ds,its_north_2068)
np.savez('MoraK_10years_averaged_NorthernSummer_2068_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_north_2068)

# 2068

data = np.load('MoraK_10years_averaged_SouthernSummer_2068_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'January 2068',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_SouthernSummer_2068_ERA5Corr_TandH.pdf')

data = np.load('MoraK_10years_averaged_NorthernSummer_2068_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'July 2068',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_NorthernSummer_2068_ERA5Corr_TandH.pdf')

print('Start south 2100')

its_south_2100 = script.Series(395,10,-15,30)
lon,lat,data_south_2100 = script.LonLatDat(ds,its_south_2100)
np.savez('MoraK_10years_averaged_SouthernSummer_2100_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_south_2100)

print('Start north 2100')

its_north_2100 = script.Series(395,10,165,30)
lon,lat,data_north_2100 = script.LonLatDat(ds,its_north_2100)
np.savez('MoraK_10years_averaged_NorthernSummer_2100_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_north_2100)

# 2100

data = np.load('MoraK_10years_averaged_SouthernSummer_2100_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'January 2100',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_SouthernSummer_2100_ERA5Corr_TandH.pdf')

data = np.load('MoraK_10years_averaged_NorthernSummer_2100_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'July 2100',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_NorthernSummer_2100_ERA5Corr_TandH.pdf')

print('Start south 2134')

its_south_2134 = script.Series(429,10,-15,30)
lon,lat,data_south_2134 = script.LonLatDat(ds,its_south_2134)
np.savez('MoraK_10years_averaged_SouthernSummer_2134_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_south_2134)

print('Start north 2134')

its_north_2134 = script.Series(429,10,165,30)
lon,lat,data_north_2134 = script.LonLatDat(ds,its_north_2134)
np.savez('MoraK_10years_averaged_NorthernSummer_2134_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_north_2134)

# 2134

data = np.load('MoraK_10years_averaged_SouthernSummer_2134_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'January 2134',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_SouthernSummer_2134_ERA5Corr_TandH.pdf')

data = np.load('MoraK_10years_averaged_NorthernSummer_2134_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'July 2134',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_NorthernSummer_2134_ERA5Corr_TandH.pdf')

print('Start south 2150')

its_south_2150 = script.Series(445,10,-15,30)
lon,lat,data_south_2150 = script.LonLatDat(ds,its_south_2150)
np.savez('MoraK_10years_averaged_SouthernSummer_2150_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_south_2150)

print('Start north 2150')

its_north_2150 = script.Series(445,10,165,30)
lon,lat,data_north_2150 = script.LonLatDat(ds,its_north_2150)
np.savez('MoraK_10years_averaged_NorthernSummer_2150_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_north_2150)

# 2150

data = np.load('MoraK_10years_averaged_SouthernSummer_2150_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'January 2150',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_SouthernSummer_2150_ERA5Corr_TandH.pdf')

data = np.load('MoraK_10years_averaged_NorthernSummer_2150_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'July 2150',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_NorthernSummer_2150_ERA5Corr_TandH.pdf')

print('Start south 2195')

its_south_2195 = script.Series(490,10,-15,30)
lon,lat,data_south_2195 = script.LonLatDat(ds,its_south_2195)
np.savez('MoraK_10years_averaged_SouthernSummer_2195_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_south_2195)

print('Start north 2195')

its_north_2195 = script.Series(490,10,165,30)
lon,lat,data_north_2195 = script.LonLatDat(ds,its_north_2195)
np.savez('MoraK_10years_averaged_NorthernSummer_2195_ERA5Corr_TandH.npz',lon=lon,lat=lat,data=data_north_2195)

# 2195

data = np.load('MoraK_10years_averaged_SouthernSummer_2195_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'January 2195',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_SouthernSummer_2195_ERA5Corr_TandH.pdf')

data = np.load('MoraK_10years_averaged_NorthernSummer_2195_ERA5Corr_TandH.npz')
lon = data['lon']
lat = data['lat']
MoraData = data['data']
plt.figure(1,figsize=(16,6))
script.Plot(lon,lat,MoraData,'July 2195',fignumber=1,cont=20,cmap=plt.get_cmap('YlOrRd'),vmin=0.,vmax=5.)
plt.savefig('/home/timothee/iLOVECLIM-GEMMES/data/Until2200_Dhigh_Pcnull_techEhigh_Instdata_ecbilt_clio/output001700/atmos/AboveMora/T2mAboveMora/WithERA5Correction/MoraK_10years_averaged_NorthernSummer_2195_ERA5Corr_TandH.pdf')


