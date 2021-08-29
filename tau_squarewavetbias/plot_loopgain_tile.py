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
inpath = '/home/data/cryo/20210313/20200120_ba1_full_lc'
lcfn   = 'LC_dark_FPU_300mK_datamode1_runsky2'
mce = 0
# inpath = '/home/data/cryo/20210313/20170106_b3_full_lc'
# mce=3
# lcfn   = 'LC_Ti_MCE%d_sky_329mK_r02'%mce
lcpath = inpath + '/' + lcfn
outpath = inpath.replace('cryo','output') + '/loopgain/' + lcfn
if not os.path.isdir(outpath):
        os.makedirs(outpath)

fnpickle = '%s/Lpg.pkl'%(outpath)
Lpg,lcpath,fitrange,calib,vnames = pkl.load(open(fnpickle,'r'))
plc = lc.get_PLC(lcpath,fitrange,calib,-1)

# get fpdata
p,ind = get_array_info('/home/data/fpdatapy/ba_fpdata_2020.mat')
tiles = np.unique(p.tile[np.where((p.mce==mce))[0]]) 

for tile in tiles:
	# get order
	xc = np.nanmax(p.det_col[np.where((p.tile==tile)*(p.mce==mce))[0]])/2.
	yc = np.nanmax(p.det_row[np.where((p.tile==tile)*(p.mce==mce))[0]])/2.
	js = np.where((p.mce==mce)*(p.tile==tile))[0]
	rs = np.full(len(js), np.inf)
	#plot_order = np.full(len(js), np.nan)
	
	for i,j in enumerate(js):
		if j+1 in ind.gd:
			mce_row = p.mce_row[j]
			mce_col = p.mce_col[j]
			x = xc
			y = 0	
		elif j+1 in ind.gl:
			mce_row = p.mce_row[j]
			mce_col = p.mce_col[j]
			x = p.det_col[j] - 0.5
			y = p.det_row[j]
		else:
			continue
		rs[i] = ((x-xc)**2+(y-yc)**2)**0.5
	podr = js[np.argsort(rs)]
	
	rRbin = np.linspace(0,1,31)
	Lpg_podr = np.full((len(podr),len(rRbin)-1),np.nan)
	for i,j in enumerate(podr):
		mce_row_podr=p.mce_row[j]
		mce_col_podr=p.mce_col[j]
		Lpg1 =abs(Lpg[mce_row_podr,mce_col_podr])
		rR1  =plc.Rdet[mce_row_podr,mce_col_podr]/plc.R[mce_row_podr,mce_col_podr]
		for k in range(len(rRbin)-1):
			Lpg_podr[i,k] = np.nanmedian(np.where((rR1>=rRbin[k])*(rR1<rRbin[k+1]), Lpg1, np.nan))

	ind_podr = np.where(np.sum(~np.isnan(Lpg_podr),1)!=0)[0]
	
	rRi,chni = np.mgrid[0:1:30*1j,0:len(podr[ind_podr]):(len(podr[ind_podr]))*1j]
	fig = plt.figure(figsize=(7, 6))
	cax = plt.pcolormesh(rRi, chni, np.log10(Lpg_podr[ind_podr]).T,vmin=0,vmax=2, cmap=plt.cm.rainbow)
	cbar = fig.colorbar(cax, ticks=[0, 1,  2])
	cbar.ax.set_yticklabels(['1', '10', '100'])
	cbar.ax.set_ylabel('loop gain', fontsize=14)
	plt.ylabel('detectors, sorted by distance to tile center', fontsize=14)
	plt.xlabel('R/Rn', fontsize=14)
	plt.title('BA, mce %d, tile %d,'%(mce,tile), fontsize=14)
	
	plt.savefig('%s/lpg_rR_tile%d.png'%(outpath,tile))






