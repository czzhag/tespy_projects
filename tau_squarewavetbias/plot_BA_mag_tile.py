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
d=loadmat('/home/data/output/202001magBA_reduc_data/magBA_di_sum.mat')

#outpath = inpath.replace('cryo','output')
onoff = 'on'
outpath = '/home/data/output/20210421/magBA_blk%s'%onoff
if not os.path.isdir(outpath):
        os.makedirs(outpath)
fnpickle = '/home/data/output/20210421/loopgain_mag_BA/tile1/Lpg.pkl'
Lpg,lcpath,fitrange,calib,vnames = pkl.load(open(fnpickle,'r'))
plc = lc.get_PLC(lcpath,fitrange,calib,-1)

# get fpdata
p,ind = get_array_info('/home/data/fpdatapy/ba_fpdata_2020.mat')
tiles = np.unique(p.tile) 

def nan_helper(y):
	return np.isnan(y), lambda z: z.nonzero()[0]

for tile in tiles:
	# get order
	xc = np.nanmax(p.det_col[np.where((p.tile==tile)*(p.rx==0))[0]])/2.
	yc = np.nanmax(p.det_row[np.where((p.tile==tile)*(p.rx==0))[0]])/2.
	js = np.where((p.tile==tile)*(p.rx==0))[0]
	rs = np.full(len(js), np.inf)
	
	for i,j in enumerate(js):
		if j+1 in ind.gd:
			x = xc
			y = 0	
			#x = p.det_col[j] - 0.5
			#y = p.det_row[j]
		elif j+1 in ind.gl:
			x = p.det_col[j] - 0.5
			y = p.det_row[j]
		else:
			continue
		rs[i] = ((x-xc)**2+(y-yc)**2)**0.5
	podr = js[np.argsort(rs)]
	
	rRbin = np.linspace(0,1,31)
	rRbinc=(rRbin[:-1]+rRbin[1:])/2.
	mag_di_podr = np.full((len(podr),len(rRbin)-1),np.nan)
	mag_dp_podr = np.full((len(podr),len(rRbin)-1),np.nan)
	Lpg_podr = np.full((len(podr),len(rRbin)-1),np.nan)
	for i,j in enumerate(podr):
		mag1 =abs(d['di_az_blk%s'%onoff][:,j]*1e6) # uA
                mag2 =mag1*d['R_blk%s'%onoff][:,j]*d['I_blk%s'%onoff][:,j]*1e6 # pW
		mce_row_podr=p.mce_row[j]
		mce_col_podr=p.mce_col[j]
		Lpg1 =abs(Lpg[mce_row_podr,mce_col_podr])
		rRlpg=plc.Rdet[mce_row_podr,mce_col_podr]/plc.R[mce_row_podr,mce_col_podr]
		rR1  =d['R_blk%s'%onoff][:,j]/d['Rn'][0,j]
		for k in range(len(rRbin)-1):
			Lpg_podr[i,k] = np.nanmedian(np.where((rRlpg>=rRbin[k])*(rRlpg<rRbin[k+1]), Lpg1, np.nan))
			mag_di_podr[i,k] = np.nanmedian(np.where((rR1>=rRbin[k])*(rR1<rRbin[k+1]), mag1, np.nan))
			mag_dp_podr[i,k] = np.nanmedian(np.where((rR1>=rRbin[k])*(rR1<rRbin[k+1]), mag2, np.nan))
		#mag_di_podr[i] = np.interp(rRbinc,rR1[::-1],mag1[::-1])
		#mag_dp_podr[i] = np.interp(rRbinc,rR1[::-1],mag2[::-1])
		nans, x= nan_helper(mag_di_podr[i])
		if not all(nans):
			mag_di_podr[i][nans]= np.interp(x(nans), x(~nans), mag_di_podr[i][~nans])
		nans, x= nan_helper(mag_dp_podr[i])
		if not all(nans):
			mag_dp_podr[i][nans]= np.interp(x(nans), x(~nans), mag_dp_podr[i][~nans])
        
	ind_podr = np.where(np.sum(~np.isnan(Lpg_podr),1)!=0)[0]
	#ind_podr = range(len(podr))
	
	rRi,chni = np.mgrid[0:1:30*1j,0:len(podr[ind_podr]):(len(podr[ind_podr]))*1j]

	fig = plt.figure(figsize=(7, 6))
	cax = plt.pcolormesh(rRi, chni, np.log10(abs(mag_di_podr[ind_podr])).T,vmin=-4,vmax=-1, cmap=plt.cm.rainbow)
	cbar = fig.colorbar(cax, ticks=[-4,-3, -2, -1])
	cbar.ax.set_yticklabels(['0.0001','0.001', '0.01', '0.1'])
	cbar.ax.set_ylabel(r'magnetic response $\Delta$ I, $\mu$A ', fontsize=14)
	plt.ylabel('detectors, sorted by distance to tile center', fontsize=14)
	plt.xlabel('R/Rn', fontsize=14)
	plt.title('BA, tile %d,'%(tile), fontsize=14)
	
	plt.savefig('%s/magBAdi_rR_tile%d.png'%(outpath,tile))

	fig = plt.figure(figsize=(7, 6))
	cax = plt.pcolormesh(rRi, chni, np.log10(abs(mag_dp_podr[ind_podr])).T,vmin=-5,vmax=-2, cmap=plt.cm.rainbow)
	cbar = fig.colorbar(cax, ticks=[-4, -3, -2])
	cbar.ax.set_yticklabels(['0.0001','0.001', '0.01'])
	cbar.ax.set_ylabel(r'magnetic response $\Delta$ P, pW ', fontsize=14)
	plt.ylabel('detectors, sorted by distance to tile center', fontsize=14)
	plt.xlabel('R/Rn', fontsize=14)
	plt.title('BA, tile %d,'%(tile), fontsize=14)
	
	plt.savefig('%s/magBAdp_rR_tile%d.png'%(outpath,tile))




