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

def main():

	doplot = True
	inpath = '/home/data/cryo/20210313/20170106_b3_full_lc'
	mce=2
	lcfn   = 'LC_Ti_MCE%d_sky_328mK_r02'%mce
	lcpath = inpath + '/' + lcfn
	outpath = inpath.replace('cryo','output') + '/loopgain_tau/' + lcfn
	if not os.path.isdir(outpath):
		os.makedirs(outpath)
	
	# setup
	import calib_B3 as calib_b3
	calib = calib_b3.calib_b3()
	fitrange = {}
	fitrange['rnti_low'] = 800
	fitrange['rnti_hgh'] = 900
	xtail = -100

	print "Process PLC and calculate loop gain."
	plc = lc.get_PLC(lcpath,fitrange,calib,-1)
	
	fnpickle = '%s/Lpg_.pkl'%(outpath) # add _ at the end for plotting
	if os.path.exists(fnpickle):
		print 'Results exist. Load %s.'%fnpickle
		Lpg,lcpath,fitrange,calib,vnames = pkl.load(open(fnpickle,'r'))
	else:
		Lpg = lc.get_loopgain(plc.Pdet,plc.Rdet,win=40)

		fnpickletau = '/home/data/output/20210324_b3_mce%d/fasttau_b3_mce%d/taus.pkl'%(mce,mce)
		biases,taus = pkl.load(open(fnpickletau,'r'))
		tau0 = 100 #ms

		#colors = plt.cm.jet(np.linspace(0,1,int(max(biases)+1)))
		#for col in range(30):
		for col in range(6):
			for row in range(1,22):
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
		
				ax0.set_ylim(0,10)
				plt.legend()
	
				ax1 = plt.subplot(gs[1])
				Rdet1 = np.where((plc.Pdet[row,col,:]>=0.5e-12), plc.Rdet[row,col,:], np.inf)
				Lpg[row,col,:]= np.where((plc.Pdet[row,col,:]>=(plc.Pdet[row,col,np.argmin(Rdet1)]*0.8)), Lpg[row,col,:], np.nan)
				Lpg1 = abs(Lpg[row,col,:xtail])
				rR1  = plc.Rdet[row,col,:xtail]/plc.R[row,col]
				ax1.plot(rR1, Lpg1,color='b', label='loop gain from P-R')
				rRbin = np.linspace(0,1,31)
				for k in range(len(rRbin)-1):
					Lpgbinmedian = np.nanmedian(np.where((rR1>=rRbin[k])*(rR1<rRbin[k+1]), Lpg1, np.nan))
					ax1.plot([rRbin[k], rRbin[k+1]],[Lpgbinmedian, Lpgbinmedian],color='orange')
				Rtaus=np.interp(biases[::-1],plc.B[col,::-1],plc.Rdet[row,col,::-1])
				ax1.scatter(Rtaus/plc.R[row,col],tau0/taus[0,::-1,row,col], s=80,c='r', marker='^',label='tau0(100 ms)/tau')	
				ax1.set_ylim(1,1e4)
				ax1.set_xlim(0,1)
				ax1.set_yscale('log')
				ax1.set_ylabel('Loop gain')
				ax1.set_xlabel('R/Rn')
				plt.grid()	
				plt.legend()
				plt.tight_layout()
		
				print '%s/lpgrow%d_col%d.png'%(outpath,row,col)		
				plt.savefig('%s/lpgrow%d_col%d.png'%(outpath,row,col))	
		
		#pkl.dump((Lpg,lcpath,fitrange,calib,'Lpg,lcpath,fitrange,calib'),open(fnpickle,'w'))

	p,ind = get_array_info('/home/data/fpdatapy/b3_fpdata_2017.mat')
	
	# do sum plots:
	# do per col
	for mce_col in range(30):
		xc = np.nanmax(p.det_col)/2.
		yc = np.nanmax(p.det_row)/2.
		R  = (xc**2+yc**2)**0.5

		fig = plt.figure(figsize=(12, 6))
		gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
		for mce_row in range(int(np.nanmax(p.mce_row))):
			if all(np.isnan(plc.Pdet[mce_row,mce_col,:xtail])):
				continue	
			j = np.where((p.mce==mce)*(p.mce_col==mce_col)*(p.mce_row==mce_row))[0][0]
			if (not j in ind.gd) and (not j in ind.gl):
				continue 
			x = p.det_col[j] - 0.5
			y = p.det_row[j] - 1.0
			r = ((x-xc)**2+(y-yc)**2)**0.5
			c = plt.cm.jet(r/R)[0]

			ax0 = plt.subplot(gs[0])
			ax0.plot(plc.Rdet[mce_row,mce_col,:xtail]/plc.R[mce_row,mce_col],
				plc.Pdet[mce_row,mce_col,:xtail]*1e12,
				color=c)
			
			ax1 = plt.subplot(gs[1])
			ax1.plot(plc.Rdet[mce_row,mce_col,:]/plc.R[mce_row,mce_col],
				abs(Lpg[mce_row,mce_col]),
				color=c)
		ax0 = plt.subplot(gs[0])
		ax0.set_xlim(0,1)
		ax0.set_ylim(0,10)
		ax0.set_ylabel('Psat, pW')
		ax0.set_xlabel('R/Rn')
		ax0.set_title('P vs R/Rn')
		ax1 = plt.subplot(gs[1])
		ax1.set_ylim(1,1e3)
		ax1.set_xlim(0,1)
		ax1.set_yscale('log')
		ax1.set_ylabel('Loop gain')
		ax1.set_xlabel('R/Rn')
		ax1.set_title('loop gain vs R/Rn, B3 mce%d col%d'%(mce,mce_col))
		plt.grid()
		plt.legend()
		plt.tight_layout()

		print '%s/lpgsum_col%d.png'%(outpath,mce_col)
        	plt.savefig('%s/lpgsum_col%d.png'%(outpath,mce_col))	

		fig = plt.figure(figsize=(12, 6))
		gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
		for mce_row in range(int(np.nanmax(p.mce_row))):
			if all(np.isnan(plc.Pdet[mce_row,mce_col,:xtail])):
				continue	
			j = np.where((p.mce==mce)*(p.mce_col==mce_col)*(p.mce_row==mce_row))[0][0]
			if (not j in ind.gd) and (not j in ind.gl):
				continue 
			x = p.det_col[j] - 0.5
			y = p.det_row[j] - 1.0
			r = ((x-xc)**2+(y-yc)**2)**0.5
			c = plt.cm.jet(r/R)[0]

			ax0 = plt.subplot(gs[0])
			ax0.plot(plc.Rdet[mce_row,mce_col,:xtail]*1e3,
				plc.Pdet[mce_row,mce_col,:xtail]*1e12,
				color=c)
			
			ax1 = plt.subplot(gs[1])
			ax1.plot(plc.Rdet[mce_row,mce_col,:]*1e3,
				abs(Lpg[mce_row,mce_col]),
				color=c)
		ax0.set_xlim(0,100)
		ax0.set_ylim(0,10)
		ax0.set_ylabel('Psat, pW')
		ax0.set_xlabel('R, mOhm')
		ax0.set_title('P vs R')
		ax1.set_ylim(1,1e3)
		ax1.set_xlim(0,100)
		ax1.set_yscale('log')
		ax1.set_ylabel('Loop gain')
		ax1.set_xlabel('R, mOhm')
		ax1.set_title('loop gain vs R, B3 mce%d col%d'%(mce,mce_col))
		plt.grid()
		plt.legend()
		plt.tight_layout()

		print '%s/lpgsum_col%d_Ohm.png'%(outpath,mce_col)
        	plt.savefig('%s/lpgsum_col%d_Ohm.png'%(outpath,mce_col))	
	

if __name__=='__main__':
	main()
	
