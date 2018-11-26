print 'Running'

from datetime import datetime
start = datetime.now()

import numpy as np 
import netCDF4 as nc
import scipy.io

# read bottom index and bottom drag coefficient from one arbitrary file 
grid_o65 = nc.Dataset('/g/data1/v45/mom/mom01v3/output065/ocean_grid.nc','r') 

x = grid_o65.variables['xu_ocean'][:]
y = grid_o65.variables['yu_ocean'][:]

# find the latitude of the Southern Ocean --- 40s~65S
idxy = np.where(np.logical_and(y>=-65,y<=-40))
lat = y[idxy]                         

# find the longitude of the Southern Ocean --- 180W~180E 
idxx = np.where(np.logical_and(x>=-280,x<=80))
lon = x[idxx]

bot = grid_o65.variables['kmu'][idxy[0][:],idxx[0][:]]
bot = np.array(bot,dtype='i')
bot = np.maximum(bot,0)

# ocean depths on u-cells
hu = grid_o65.variables['hu'][idxy[0][:],idxx[0][:]]

rho0 = 1035
cd = grid_o65.variables['drag_coeff'][0,idxy[0][:],idxx[0][:]] # cd has only one t-dimension for one-month
cd = cd.filled(np.nan)
#---

ny = np.size(idxy[0]) # ny = 426
nx = np.size(idxx[0]) # nx = 3600

ub = np.zeros((1,ny,nx))      # bottom velocity - u (the first dimension denotes time)
vb = np.zeros((1,ny,nx)) 
ubt = np.full((1,ny,nx),np.nan)
vbt = np.full((1,ny,nx),np.nan)
bblx = np.zeros((1,ny,nx))
bbly = np.zeros((1,ny,nx))
bblxt = np.full((1,ny,nx),np.nan)
bblyt = np.full((1,ny,nx),np.nan)
MKE = np.full((ny,nx),np.nan)
mTKE = np.full((ny,nx),np.nan)
mEKE = np.full((ny,nx),np.nan)
Mdisp = np.full((ny,nx),np.nan)
mdisp = np.full((ny,nx),np.nan)
mEdisp = np.full((ny,nx),np.nan)
Cdisp = np.full((ny,nx),np.nan)

# get bottom velocity
for opt in xrange(80,92,1):
#for opt in xrange(65,77,1):
#for opt in xrange(65,66,1):
  print '--------'
  print opt
# print '/g/data1/v45/mom/mom01v3/output0%d/u_snap.nc' % opt
  dir_opt_u = '/g/data1/v45/mom/mom01v3/output0%d/u_snap.nc' % opt
  dir_opt_v = '/g/data1/v45/mom/mom01v3/output0%d/v_snap.nc' % opt

  data_u = nc.Dataset(dir_opt_u,'r')
  data_v = nc.Dataset(dir_opt_v,'r')
  
  t = data_u.variables['time'][:]
  print np.shape(t)

# for k in xrange(np.size(t)):
  for k in xrange(0,np.size(t),20):
# for k in xrange(3):
    print k
    u = data_u.variables['u'][k,:,idxy[0][:],idxx[0][:]]
    v = data_v.variables['v'][k,:,idxy[0][:],idxx[0][:]]

    # change the masked array into ndarray, fill all the masked array with np.nan
    u = u.filled(np.nan)
    v = v.filled(np.nan)

    ub_choose = np.array([u[bot[I]-1][I] for I in np.ndindex(bot.shape)])
    ubt[0,:,:] = np.reshape(ub_choose,(ny,nx))
    vb_choose = np.array([v[bot[I]-1][I] for I in np.ndindex(bot.shape)])
    vbt[0,:,:] = np.reshape(vb_choose,(ny,nx))
    # bot[I]: numbers of level
    # bot[I]-1: the index of bottom level in Python
    # u[bot[I]-1]: go to bottom level in u
    # u[bot[I]-1][I]: find the bottom value indicated by I in the bottom level

    ub = np.append(ub,ubt,0) 
    vb = np.append(vb,vbt,0)  

print 'Finish finding bottom velocity'
#---


# decompose bottom velocity
# the first dimension of ub and vb are one more than the time steps
# get rid of the first zero values
ub = ub[1:,:,:]  
vb = vb[1:,:,:]
um = np.mean(ub,axis=0)
vm = np.mean(vb,axis=0)
up = ub - um
vp = vb - vm

nt = np.size(ub,0)
Ub = np.full((nt,ny,nx),np.nan)
TKE = np.full((nt,ny,nx),np.nan)
EKE = np.full((nt,ny,nx),np.nan)
disp = np.full((nt,ny,nx),np.nan)
Edisp = np.full((nt,ny,nx),np.nan)

Ub = (ub**2 + vb**2)**0.5
MKE = 0.5*(um**2 + vm **2)
TKE = 0.5*(ub**2 + vb **2)
EKE = 0.5*(up**2 + vp **2)

disp = -rho0*cd*(ub**2 + vb**2)**1.5
Cdisp = -rho0*cd*(um**2 + vm**2)**1.5

mTKE = np.mean(TKE,axis=0)
mEKE = np.mean(EKE,axis=0)
mdisp = np.mean(disp,axis=0)

for l in xrange(nt):  
  bblxt[0,:,:] = -rho0*cd*Ub[l,:,:]*ub[l,:,:]
  bblyt[0,:,:] = -rho0*cd*Ub[l,:,:]*vb[l,:,:]
  bblx = np.append(bblx,bblxt,0)
  bbly = np.append(bbly,bblyt,0)  
  
bblx = bblx[1:,:,:]
bbly = bbly[1:,:,:]
mbblx = np.mean(bblx,axis=0)
mbbly = np.mean(bbly,axis=0)
bblxp = bblx - mbblx
bblyp = bbly - mbbly

Mdisp = mbblx*um + mbbly*vm
Edisp = bblxp*up + bblyp*vp
mEdisp = np.mean(Edisp,axis=0)

print 'Finish decomposing and KE computation'

time_diff = datetime.now()-start

# save results in a .nc file
#fileobj = nc.Dataset('/short/v45/lxy581/pyfiles/decomp/decomp_KE_disp_SO_10_160804.nc','w')
fileobj = nc.Dataset('/short/v45/lxy581/pyfiles/decomp/decomp_KE_disp_SO_10_171117_1yr.nc','w')

fileobj.createDimension('lat',len(lat)) # ny = len(lat)
fileobj.createDimension('lon',len(lon)) # nx = len(lon)
fileobj.createDimension('time',nt) # nt = len(time)

lat_var = fileobj.createVariable('lat','f',('lat',))
lon_var = fileobj.createVariable('lon','f',('lon',))
ub_var = fileobj.createVariable('ub','f',('time','lat','lon'),fill_value=-1e+20)
vb_var = fileobj.createVariable('vb','f',('time','lat','lon'),fill_value=-1e+20)
um_var = fileobj.createVariable('um','f',('lat','lon'),fill_value=-1e+20)
vm_var = fileobj.createVariable('vm','f',('lat','lon'),fill_value=-1e+20)
MKE_var = fileobj.createVariable('MKE','f',('lat','lon'),fill_value=-1e+20)
mTKE_var = fileobj.createVariable('mTKE','f',('lat','lon'),fill_value=-1e+20)
mEKE_var = fileobj.createVariable('mEKE','f',('lat','lon'),fill_value=-1e+20)
mdisp_var = fileobj.createVariable('mdisp','f',('lat','lon'),fill_value=-1e+20)
mEdisp_var = fileobj.createVariable('mEdisp','f',('lat','lon'),fill_value=-1e+20)
Mdisp_var = fileobj.createVariable('Mdisp','f',('lat','lon'),fill_value=-1e+20)
Cdisp_var = fileobj.createVariable('Cdisp','f',('lat','lon'),fill_value=-1e+20)

hu_var = fileobj.createVariable('hu','f',('lat','lon'),fill_value=-1e+20)
mbblx_var = fileobj.createVariable('mbblx','f',('lat','lon'),fill_value=-1e+20)
mbbly_var = fileobj.createVariable('mbbly','f',('lat','lon'),fill_value=-1e+20)

lat_var[:] = lat[:]
lon_var[:] = lon[:]
ub_var[:,:,:] = ub[:,:,:]
vb_var[:,:,:] = vb[:,:,:]
um_var[:,:]=um[:,:]
vm_var[:,:]=vm[:,:]
MKE_var[:,:] = MKE[:,:]
mTKE_var[:,:] = mTKE[:,:]
mEKE_var[:,:] = mEKE[:,:]
mdisp_var[:,:] = mdisp[:,:]
mEdisp_var[:,:] = mEdisp[:,:]
Mdisp_var[:,:] = Mdisp[:,:]
Cdisp_var[:,:] = Cdisp[:,:]

hu_var[:,:] = hu[:,:]
mbblx_var[:,:] = mbblx[:,:]
mbbly_var[:,:] = mbbly[:,:]

lat_var.units='degree N'
lon_var.units='degree E' 
ub_var.units='m/s'
vb_var.units='m/s'
um_var.units='m/s'
vm_var.units='m/s'
MKE_var.units='m^2/s^2'
mTKE_var.units='m^2/s^2'
mEKE_var.units='m^2/s^2'
mdisp_var.units='W/m^2'
mEdisp_var.units='W/m^2'
Mdisp_var.units='W/m^2'
Cdisp_var.units='W/m^2'

hu_var.units='m'
mbblx_var.units='N/m^2'
mbbly_var.units='N/m^2'

lat_var.long_name='Latitude'
lon_var.long_name='Longitude'
ub_var.long_name='zonal velocity with time dimension'
vb_var.long_name='meridional velocity with time dimension'
um_var.long_name='time-mean zonal velocity'
vm_var.long_name='time_mean meridional velocity'
MKE_var.long_name='Kinetic energy of bottom mean flow'
mTKE_var.long_name='time-mean Kinetic Energy of bottom total flow'
mEKE_var.long_name='time-mean Kinetic Energy of bottom eddy flow'
mdisp_var.long_name='time-mean dissipation rate due to bottom boundary layer drag - total'
mEdisp_var.long_name='time-mean dissipation rate due to bottom boundary layer drag - eddy'
Mdisp_var.long_name='time-mean dissipation rate due to bottom boundary layer drag - mean'
Cdisp_var.long_name='contribution from mean flow to dissipation to mean flow'

hu_var.long_name='ocean depth'
mbblx_var.long_name='time-mean tao bottom x'
mbbly_var.long_name='time-mean tao bottom y'

print 'Finished! ^.^'

# changes made on 21 July, 2016:
# (1) add density into all KE computations;
# (2) save hu(ocean depth, topo) in the region;
# (3) save mbblx and mbbly;

# changes made on 3 Aug, 2016:
# (1) delete density from all KE calculations;
# (2) save another quantity - Cdisp for future percentage calculation.
#     1-Cdisp_integrated/Mdisp_integrated = eddy contribution percentage.

# Rerun
# because of the stupid mistake about disp and Cdisp
