import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap
import cmocean

plt.close('all')

# TBBL
data = nc.Dataset('/short/v45/lxy581/pyfiles/decomp/decomp_KE_disp_SO_10_171117_1yr.nc','r')

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

# mean drag [TBBL]
mbblx = data.variables['mbblx'][:,:]
mbbly = data.variables['mbbly'][:,:]

mbbl = np.sqrt(mbblx**2+mbbly**2)
mbbl0 = np.log10(mbbl)

# mean drag [LW]
nf_x = nf_fr4['miwx']
nf_y = nf_fr4['miwy']
nf_x = np.transpose(nf_x)
nf_y = np.transpose(nf_y)
nf_m = np.sqrt(nf_x**2+nf_y**2) # mean drag magnitude
nf_l = np.log10(nf_m) # log 10 scale

g_x = g_fr4['miwx']
g_y = g_fr4['miwy']
g_x = np.transpose(g_x)
g_y = np.transpose(g_y)
g_m = np.sqrt(g_x**2+g_y**2) # mean drag magnitude
g_l = np.log10(g_m) # log 10 scale

ga_x = ga_fr4['miwx']
ga_y = ga_fr4['miwy']
ga_x = np.transpose(ga_x)
ga_y = np.transpose(ga_y)
ga_m = np.sqrt(ga_x**2+ga_y**2) # mean drag magnitude
ga_l = np.log10(ga_m) # log 10 scale

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

DP_g[:,:] = g_m[ind2[0][:],:][:,ind1[0][:]]
DP_ga[:,:] = ga_m[ind2[0][:],:][:,ind1[0][:]]
DP_nf[:,:] = nf_m[ind2[0][:],:][:,ind1[0][:]]
DP_tb[:,:] = mbbl[ind2[0][:],:][:,ind1[0][:]]

print 'Start plotting'

plt.figure(1,figsize=(12,9))
#N m-2
#drag_level = np.arange(0,0.08+0.001,0.001)
#drag_ticks = np.arange(0,0.08+0.02,0.02)
#log10(N m-2)
drag_log_level = np.arange(-4.5,-0.5+0.05,0.05)
drag_log_ticks = np.arange(-4.5,-0.5+1.0,1.0)

# G2010
plt.subplot(411)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
map.contourf(x,y,g_l,cmap=cmocean.cm.thermal,levels=drag_log_level,extend='both')
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.72,0.76,0.15])
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
tx4,ty4 = map(lonn4+2,latn4+0.5)
plt.title('a) G2010',fontsize=20,loc='left')

print 'Finish 411'

# GA2010
plt.subplot(412)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
map.contourf(x,y,ga_l,cmap=cmocean.cm.thermal,levels=drag_log_level,extend='both')
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.50,0.76,0.15])
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
tx4,ty4 = map(lonn4+2,latn4+0.5)
plt.title('b) GA2010',fontsize=20,loc='left')

print 'Finish 412'

# NF2011
plt.subplot(413)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
map.contourf(x,y,nf_l,cmap=cmocean.cm.thermal,levels=drag_log_level,extend='both')
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.28,0.76,0.15])
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
tx4,ty4 = map(lonn4+2,latn4+0.5)
plt.title('c) NF2011',fontsize=20,loc='left')

print 'Finish 413'

# TBBL
plt.subplot(414)
map = Basemap(projection='merc',llcrnrlat=-65,urcrnrlat=-40,llcrnrlon=-280,urcrnrlon=80,resolution='l',fix_aspect=False)
map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color='0.8')
x,y = map(LON,LAT)
mcf = map.contourf(x,y,mbbl0,cmap=cmocean.cm.thermal,levels=drag_log_level,extend='both')
cb = plt.colorbar(mcf,ticks=drag_log_ticks,shrink=0.8)
map.drawparallels(np.arange(-60,-45+5,5),labels=[1,0,0,0],fontsize=16)
map.drawmeridians(np.arange(-240,60+60,60),labels=[0,0,0,1],fontsize=16)
plt.gca().set_position([0.1,0.06,0.76,0.15])
x4,y4 = map(lons4,lats4)
map.plot(x4,y4,linewidth=2, color='k')
tx4,ty4 = map(lonn4+2,latn4+0.5)
cb.ax.set_position([0.88,0.18,0.08,0.6])
cb.set_label('Mean drag, log$_{10}$ (N m$^{-2}$)',y=0.5,fontsize=20)
#cb.set_label('Mean drag, N m$^{-2}$',y=0.5,fontsize=20)
cb.ax.tick_params(labelsize=16)
plt.title('d) TBBL',fontsize=20,loc='left')

print 'Finish 414'
#plt.savefig('/short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_mean_drag_lin.png',dpi=100)
#plt.savefig('/short/v45/lxy581/lw/corr_fin_topog/lw_SO_10_b5_1yr_mean_drag_log.png',dpi=100)
plt.savefig('/short/v45/lxy581/lw/corr_fin_topog/LeeWaveMeanDragSO10b51yrlog.png',dpi=100)

plt.show()

print 'Finished!'
