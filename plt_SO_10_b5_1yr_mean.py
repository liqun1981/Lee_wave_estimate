import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap
import cmocean

plt.close('all')

# TBBL
data = nc.Dataset('/short/v45/lxy581/pyfiles/decomp/decomp_KE_disp_SO_10_171117_1yr.nc','r')
data_g = nc.Dataset('/g/data1/v45/mom/mom01v3/output065/ocean_grid.nc','r')

# NF2011
nf_fr4 = scipy.io.loadmat('/short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_NF11_180129_fr4.mat')

# G2010
g_fr4 = scipy.io.loadmat('/short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_g2010_180129_fr4.mat')

# GA2010
ga_fr4 = scipy.io.loadmat('/short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_ga2010_180129_fr4.mat')

# read - grid
lon = nf_fr4['lon']
lat = nf_fr4['lat']
[LON,LAT] = np.meshgrid(lon,lat)

# read - energy loss from the mean flow
tb_mean = data.variables['Mdisp'][:,:] # TBBL
nf_mean = nf_fr4['Mdisp_nl']
g_mean  = g_fr4['Mdisp_nl']
ga_mean = ga_fr4['Mdisp_nl']

# transpose
nf_mean = np.transpose(nf_mean)
g_mean  = np.transpose(g_mean)
ga_mean = np.transpose(ga_mean)

nf_mean[nf_mean==0]=np.nan
g_mean[g_mean==0]=np.nan
ga_mean[ga_mean==0]=np.nan

# log10
tb_log = np.log10(-tb_mean)+3
nf_log = np.log10(-nf_mean)+3
g_log = np.log10(g_mean)+3
ga_log = np.log10(ga_mean)+3

# area of grid cells
data_u = nc.Dataset('/g/data1/v45/mom/mom01v3/output065/u_snap.nc','r')
data_g = nc.Dataset('/g/data1/v45/mom/mom01v3/output065/ocean_grid.nc','r')

x = data_u.variables['xu_ocean'][:]
y = data_u.variables['yu_ocean'][:]

# find the latitude of the SO --- 40s~65S
idxy = np.where(np.logical_and(y>=-65,y<=-40))                        

# find the longitude of the SO --- 180W~180E 
idxx = np.where(np.logical_and(x>=-280,x<=80))

# get the area of velocity cells in the SO
a = data_g.variables['area_u'][idxy[0][:],idxx[0][:]] # velocity cell area
a = a.filled(np.nan)

# rate: w m-2 to W
tb_rt = -tb_mean * a
nf_rt = nf_mean * a
g_rt = g_mean * a
ga_rt = ga_mean * a

# mark - Scotia Sea/Drake Passage
lonn4 = -70
lonm4 = -50
latn4 = -62
latm4 = -55

lats4 = [latn4,latm4,latm4,latn4,latn4]
lons4 = [lonn4,lonn4,lonm4,lonm4,lonn4]

ind1 = np.where((lon>=lonn4) & (lon<=lonm4))
ind2 = np.where((lat>=latn4) & (lat<=latm4))

DP_nx = np.size(ind1[0][:])
DP_ny = np.size(ind2[0][:])

DP_g = np.zeros((DP_ny,DP_nx))
DP_ga = np.zeros((DP_ny,DP_nx))
DP_nf = np.zeros((DP_ny,DP_nx))
DP_tb = np.zeros((DP_ny,DP_nx))

DP_g[:,:] = g_mean[ind2[0][:],:][:,ind1[0][:]]
DP_ga[:,:] = ga_mean[ind2[0][:],:][:,ind1[0][:]]
DP_nf[:,:] = -nf_mean[ind2[0][:],:][:,ind1[0][:]]
DP_tb[:,:] = -tb_mean[ind2[0][:],:][:,ind1[0][:]]

print 'Start plotting'

plt.figure(1,figsize=(12,9))

# log - mW/m^2
disp_log_level = np.arange(-4,2+0.1,0.1)
disp_log_ticks = np.arange(-4,2+2,2)

# G2010
plt.subplot(411)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
map.contourf(x,y,g_log,cmap=cmocean.cm.thermal,levels=disp_log_level,extend='both')
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.72,0.76,0.15])
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
plt.title('a) G2010',fontsize=20,loc='left')
print 'Finish 411'

# GA2010
plt.subplot(412)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
map.contourf(x,y,ga_log,cmap=cmocean.cm.thermal,levels=disp_log_level,extend='both')
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.50,0.76,0.15])
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
plt.title('b) GA2010',fontsize=20,loc='left')
print 'Finish 412'

# NF2011
plt.subplot(413)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
map.contourf(x,y,nf_log,cmap=cmocean.cm.thermal,levels=disp_log_level,extend='both')
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.28,0.76,0.15])
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
plt.title('c) NF2011',fontsize=20,loc='left')
print 'Finish 413'

# TBBL
plt.subplot(414)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
mcf = map.contourf(x,y,tb_log,cmap=cmocean.cm.thermal,levels=disp_log_level,extend='both')
cb = plt.colorbar(mcf,ticks=disp_log_ticks,shrink=0.8)
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.06,0.76,0.15])
cb.ax.set_position([0.88,0.18,0.08,0.6])
cb.set_label('Energy conversion rate, log$_{10}$ (mW m$^{-2}$)',y=0.5,fontsize=20)
cb.ax.tick_params(labelsize=16)
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
plt.title('d) TBBL',fontsize=20,loc='left')
print 'Finish 414'

# plt.savefig('/short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_mean.png',dpi=100)
plt.savefig('/short/v45/lxy581/lw/corr_fin_topog/LeeWavesSO10b51yrMean.png',dpi=100)

plt.show()
