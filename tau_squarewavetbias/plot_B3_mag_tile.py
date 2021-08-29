# 4/13/2021
# to Correlate loop gain, tau and magnetic pick-up
import mce_data
import numpy as np
import os,sys
sys.path.insert(0, "../tespy")
import lcdata as lc
from get_array_info import get_array_info
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.io import loadmat
import cPickle as pkl
from matplotlib import gridspec
import glob

# get data, plc and loopgain
doplot = True
d=loadmat('/home/data/output/201911magB3_reduc_data/magB3_di_sum.mat')

#outpath = inpath.replace('cryo','output')
onoff = 'on'
outpath = '/home/data/output/20210413/magB3_blk%s'%onoff
if not os.path.isdir(outpath):
        os.makedirs(outpath)

# get fpdata
p,ind = get_array_info('/home/data/fpdatapy/b3_fpdata_2017.mat')
tiles = np.unique(p.tile) 

for tile in tiles:
	# get order
	xc = np.nanmax(p.det_col[np.where(p.tile==tile)[0]])/2.
	yc = np.nanmax(p.det_row[np.where(p.tile==tile)[0]])/2.
	js = np.where(p.tile==tile)[0]
	rs = np.full(len(js), np.inf)
	
	for i,j in enumerate(js):
		if j+1 in ind.gd:
			#x = xc
			#y = 0	
			x = p.det_col[j] - 0.5
			y = p.det_row[j]
		elif j+1 in ind.gl:
			x = p.det_col[j] - 0.5
			y = p.det_row[j]
		else:
			continue
		rs[i] = ((x-xc)**2+(y-yc)**2)**0.5
	podr = js[np.argsort(rs)]
	
	rRbin = np.linspace(0,1,31)
	mag_di_podr = np.full((len(podr),len(rRbin)-1),np.nan)
	mag_dp_podr = np.full((len(podr),len(rRbin)-1),np.nan)
	for i,j in enumerate(podr):
		mag1 =abs(d['di_az_blk%s'%onoff][:,j]*1e6) # uA
                mag2 =mag1*d['R_blk%s'%onoff][:,j]*d['I_blk%s'%onoff][:,j]*1e6 # pW
		rR1  =d['R_blk%s'%onoff][:,j]/d['Rn'][0,j]
		for k in range(len(rRbin)-1):
			mag_di_podr[i,k] = np.nanmedian(np.where((rR1>=rRbin[k])*(rR1<rRbin[k+1]), mag1, np.nan))
			mag_dp_podr[i,k] = np.nanmedian(np.where((rR1>=rRbin[k])*(rR1<rRbin[k+1]), mag2, np.nan))
        
	#ind_podr = np.where(np.sum(~np.isnan(mag_di_podr),1)!=0)[0]
	ind_podr = range(len(podr))
	
	rRi,chni = np.mgrid[0:1:30*1j,0:len(podr[ind_podr]):(len(podr[ind_podr]))*1j]

	fig = plt.figure(figsize=(7, 6))
	cax = plt.pcolormesh(rRi, chni, np.log10(mag_di_podr[ind_podr]).T,vmin=-3,vmax=-1, cmap=plt.cm.rainbow)
	cbar = fig.colorbar(cax, ticks=[-3, -2, -1])
	cbar.ax.set_yticklabels(['0.001', '0.01', '0.1'])
	cbar.ax.set_ylabel(r'magnetic response $\Delta$ I, $\mu$A ', fontsize=14)
	plt.ylabel('detectors, sorted by distance to tile center', fontsize=14)
	plt.xlabel('R/Rn', fontsize=14)
	plt.title('B3, tile %d,'%(tile), fontsize=14)
	
	plt.savefig('%s/magB3di_rR_tile%d.png'%(outpath,tile))

	fig = plt.figure(figsize=(7, 6))
	cax = plt.pcolormesh(rRi, chni, np.log10(mag_dp_podr[ind_podr]).T,vmin=-5,vmax=-2, cmap=plt.cm.rainbow)
	cbar = fig.colorbar(cax, ticks=[-4, -3, -2])
	cbar.ax.set_yticklabels(['0.0001','0.001', '0.01'])
	cbar.ax.set_ylabel(r'magnetic response $\Delta$ P, pW ', fontsize=14)
	plt.ylabel('detectors, sorted by distance to tile center', fontsize=14)
	plt.xlabel('R/Rn', fontsize=14)
	plt.title('B3, tile %d,'%(tile), fontsize=14)
	
	plt.savefig('%s/magB3dp_rR_tile%d.png'%(outpath,tile))






