import mce_data
import numpy as np
import os,sys
sys.path.insert(0, "../tespy")
import lcdata as lc
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import cPickle as pkl
from matplotlib import gridspec
import glob 

def main():

	doplot = True
	inpath = '/home/data/cryo/20210309mce1'
	lcfn   = '1615425034_plc_bias'
	#inpath = '/home/data/cryo/20210309'
	#lcfn   = '1615274088_plc_bias'
	lcpath = inpath + '/' + lcfn
	outpath = inpath.replace('cryo','output') + '/loopgain/' + lcfn
	if not os.path.isdir(outpath):
		os.makedirs(outpath)
	
	# setup
	import calib_keck as calib_keck
	calib = calib_keck.calib_keck()
	fitrange = {}
	fitrange['rnti_low'] = 8000
	fitrange['rnti_hgh'] = 9500
	xtail = -500

	print "Process PLC and calculate loop gain."
	plc = lc.get_PLC(lcpath,fitrange,calib,1)
	Lpg = lc.get_loopgain(plc.Pdet,plc.Rdet,win=200)

	Pflat = np.full_like(plc.Pdet,np.nan)
	Rflat = np.full_like(plc.Rdet,np.nan)
	Bflat = np.full_like(plc.Rdet,np.nan)

	fnpickle = '/home/data/output/20210309mce1/fasttaurx1_1/taus.pkl'
	biases,taus = pkl.load(open(fnpickle,'r'))
	#tau0 = np.full_like(taus, np.nan)
	tau0 = 100 #ms

	colors = plt.cm.jet(np.linspace(0,1,int(max(biases)+1)))
	#for row in [12]:
	for row in range(33):
	#	for col in [2]:
		for col in range(16):

			if not doplot:
				continue
			if all(np.isnan(plc.Pdet[row,col,:xtail])):
				continue	
			
			fig = plt.figure(figsize=(6, 8))
			gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
			ax0 = plt.subplot(gs[0])
			ax0.plot(plc.Rdet[row,col,:xtail]/plc.R[row,col],plc.Pdet[row,col,:xtail]*1e12,color='b')
			ax0.set_xlim(0,1)
			ax0.set_ylabel('Psat, pW')
			ax0.set_xlabel(' ')
			ax0.set_title('P/loop gain vs R/Rn, row %d, col%d'%(row,col))
		
			c,b=np.histogram(plc.Pdet[row,col,:xtail]*1e12,bins=50)
			b1 = np.sum(np.where(c==max(c),1,0)*b[:-1])
			b2 = np.sum(np.roll(np.where(c==max(c),1,0),1)*b[:-1])
			psat = np.nanmedian(np.where(((plc.Pdet[row,col,:xtail]*1e12<b2)*(plc.Pdet[row,col,:xtail]*1e12>b1)),plc.Pdet[row,col,:xtail]*1e12,np.nan))
			b1 = psat-0.05
			b2 = psat+0.05
			Rflat[row,col,:xtail] = np.where(((plc.Pdet[row,col,:xtail]*1e12<b2)*(plc.Pdet[row,col,:xtail]*1e12>b1)),plc.Rdet[row,col,:xtail],np.nan)
			Pflat[row,col,:xtail] = np.where(((plc.Pdet[row,col,:xtail]*1e12<b2)*(plc.Pdet[row,col,:xtail]*1e12>b1)),plc.Pdet[row,col,:xtail],np.nan)
			Bflat[row,col,:xtail] = np.where(((plc.Pdet[row,col,:xtail]*1e12<b2)*(plc.Pdet[row,col,:xtail]*1e12>b1)),plc.B[col,:xtail],np.nan)

			ax0.axvline(x=np.nanmin(Rflat[row,col])/plc.R[row,col],color='k',linestyle='--',label='flat PR')
			ax0.axvline(x=np.nanmax(Rflat[row,col])/plc.R[row,col],color='k',linestyle='--')
			ax0.set_ylim(0,25)
			plt.legend()

			ax1 = plt.subplot(gs[1])
			ax1.plot(plc.Rdet[row,col,:xtail]/plc.R[row,col],abs(Lpg[row,col,:xtail]),color='b', label='loop gain from P-R')
			Rtaus=np.interp(biases[::-1],plc.B[col,::-1],plc.Rdet[row,col,::-1])
			ax1.scatter(Rtaus/plc.R[row,col],tau0/taus[0,::-1,row,col], s=80,c='r', marker='^',label='tau0(%d ms)/tau'%tau0)	
			ax1.axvline(x=np.nanmin(Rflat[row,col])/plc.R[row,col],color='k',linestyle='--')
			ax1.axvline(x=np.nanmax(Rflat[row,col])/plc.R[row,col],color='k',linestyle='--')
			ax1.set_ylim(1,1e4)
			ax1.set_xlim(0,1)
			ax1.set_yscale('log')
			ax1.set_ylabel('Loop gain')
			ax1.set_xlabel('R/Rn')
			plt.grid()	
			plt.legend(loc=1)
			plt.tight_layout()
			
			print '%s/lpgrow%d_col%d.png'%(outpath,row,col)
			plt.savefig('%s/lpgrow%d_col%d.png'%(outpath,row,col))	

if __name__=='__main__':
	main()
	
