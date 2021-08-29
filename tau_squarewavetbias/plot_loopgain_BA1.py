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
	inpath = '/home/data/cryo/20210313/20200120_ba1_full_lc'
	lcfn   = 'LC_dark_FPU_300mK_datamode1_runsky2'
	lcpath = inpath + '/' + lcfn
	outpath = inpath.replace('cryo','output') + '/loopgain/' + lcfn
	if not os.path.isdir(outpath):
		os.makedirs(outpath)
	
	# setup
	import calib_BA1 as calib_ba1
	calib = calib_ba1.calib_ba1()
	fitrange = {}
	fitrange['rnti_low'] = 850
	fitrange['rnti_hgh'] = 950
	xtail = -10

	print "Process PLC and calculate loop gain."
	plc = lc.get_PLC(lcpath,fitrange,calib,-1)

	Pflat = np.full_like(plc.Pdet,np.nan)
	Rflat = np.full_like(plc.Rdet,np.nan)
	Bflat = np.full_like(plc.Rdet,np.nan)


	
	fnpickle = '%s/Lpg.pkl'%(outpath)
	if os.path.exists(fnpickle):
		print 'Results exist. Load %s.'%fnpickle
		Lpg,lcpath,fitrange,calib,vnames = pkl.load(open(fnpickle,'r'))
	else:
		Lpg = lc.get_loopgain(plc.Pdet,plc.Rdet,win=40)

		#fnpickletau = '/home/data/output/20210219/fasttau3/taus.pkl'
		#biases,taus = pkl.load(open(fnpickle,'r'))
		tau0 = 100 #ms

		#colors = plt.cm.jet(np.linspace(0,1,int(max(biases)+1)))
		for row in range(1,33):
			for col in range(24):
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
				ax0.set_ylim(0,7)
				plt.legend()
	
				ax1 = plt.subplot(gs[1])
				Rdet1 = np.where((plc.Pdet[row,col,:]>=0.01e-12), plc.Rdet[row,col,:], np.inf)
				Lpg[row,col,:]= np.where((plc.Pdet[row,col,:]>=(plc.Pdet[row,col,np.argmin(Rdet1)]*0.8)), Lpg[row,col,:], np.nan)
				Lpg1 = abs(Lpg[row,col,:xtail])
				rR1  = plc.Rdet[row,col,:xtail]/plc.R[row,col]
				ax1.plot(rR1, Lpg1,color='b', label='loop gain from P-R')
				rRbin = np.linspace(0,1,31)
				for k in range(len(rRbin)-1):
					Lpgbinmedian = np.nanmedian(np.where((rR1>=rRbin[k])*(rR1<rRbin[k+1]), Lpg1, np.nan))
					ax1.plot([rRbin[k], rRbin[k+1]],[Lpgbinmedian, Lpgbinmedian],color='orange')
				ax1.axvline(x=np.nanmin(Rflat[row,col])/plc.R[row,col],color='k',linestyle='--',label='flat PR')
				ax1.axvline(x=np.nanmax(Rflat[row,col])/plc.R[row,col],color='k',linestyle='--')
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
		
		pkl.dump((Lpg,lcpath,fitrange,calib,'Lpg,lcpath,fitrange,calib'),open(fnpickle,'w'))

	p,ind = get_array_info('/home/data/fpdatapy/ba_fpdata_2020.mat')
	
	# do sum plots:
	# do per col
	for mce_col in range(24):
		xc = np.nanmax(p.det_col[np.where((p.tile==int(mce_col/2)+1)*(p.rx==0))[0]])/2.
		yc = np.nanmax(p.det_row[np.where((p.tile==int(mce_col/2)+1)*(p.rx==0))[0]])/2.
		R  = (xc**2+yc**2)**0.5


		fig = plt.figure(figsize=(12, 6))
		gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
		for mce_row in range(int(np.nanmax(p.mce_row))):
			if all(np.isnan(plc.Pdet[mce_row,mce_col,:xtail])):
				continue	
			j = np.where((p.mce_col==mce_col)*(p.mce_row==mce_row)*(p.rx==0))[0][0]
			if j+1 in ind.gd:
				x = xc
				y = 0
			elif j+1 in ind.gl:
				x = p.det_col[j][0] - 0.5
				y = p.det_row[j][0] 
			else:
				continue
			r = ((x-xc)**2+(y-yc)**2)**0.5
			c = plt.cm.jet(r/R)

			ax0 = plt.subplot(gs[0])
			ax0.plot(plc.Rdet[mce_row,mce_col,:xtail]/plc.R[mce_row,mce_col],
				plc.Pdet[mce_row,mce_col,:xtail]*1e12,
				color=c)
			
			ax1 = plt.subplot(gs[1])
			ax1.plot(plc.Rdet[mce_row,mce_col,:]/plc.R[mce_row,mce_col],
				abs(Lpg[mce_row,mce_col]),
				color=c)
		ax0.set_xlim(0,1)
		ax0.set_ylim(0,7)
		ax0.set_ylabel('Psat, pW')
		ax0.set_xlabel('R/Rn')
		ax0.set_title('P vs R/Rn')
		ax1.set_ylim(1,1e3)
		ax1.set_xlim(0,1)
		ax1.set_yscale('log')
		ax1.set_ylabel('Loop gain')
		ax1.set_xlabel('R/Rn')
		ax1.set_title('loop gain vs R/Rn, BA1 col%d'%(mce_col))
		plt.grid()
		plt.legend()
		plt.tight_layout()

		try:
			print '%s/lpgsum_col%d.png'%(outpath,mce_col)
        		plt.savefig('%s/lpgsum_col%d.png'%(outpath,mce_col))	
		except:
			continue	

		fig = plt.figure(figsize=(12, 6))
		gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
		for mce_row in range(int(np.nanmax(p.mce_row))):
			if all(np.isnan(plc.Pdet[mce_row,mce_col,:xtail])):
				continue	
			j = np.where((p.mce_col==mce_col)*(p.mce_row==mce_row)*(p.rx==0))[0][0]
			if j+1 in ind.gd:
				x = xc
				y = 0
			elif j+1 in ind.gl:
				x = p.det_col[j][0] - 0.5
				y = p.det_row[j][0] 
			else:
				continue
			r = ((x-xc)**2+(y-yc)**2)**0.5
			c = plt.cm.jet(r/R)

			ax0 = plt.subplot(gs[0])
			ax0.plot(plc.Rdet[mce_row,mce_col,:xtail]*1e3,
				plc.Pdet[mce_row,mce_col,:xtail]*1e12,
				color=c)
			
			ax1 = plt.subplot(gs[1])
			ax1.plot(plc.Rdet[mce_row,mce_col,:]*1e3,
				abs(Lpg[mce_row,mce_col]),
				color=c)
		ax0.set_xlim(0,200)
		ax0.set_ylim(0,7)
		ax0.set_ylabel('Psat, pW')
		ax0.set_xlabel('R, mOhm')
		ax0.set_title('P vs R')
		ax1.set_ylim(1,1e3)
		ax1.set_xlim(0,200)
		ax1.set_yscale('log')
		ax1.set_ylabel('Loop gain')
		ax1.set_xlabel('R, mOhm')
		ax1.set_title('loop gain vs R, BA1 col%d'%(mce_col))
		plt.grid()
		plt.legend()
		plt.tight_layout()

		try:
			print '%s/lpgsum_col%d_Ohm.png'%(outpath,mce_col)
        		plt.savefig('%s/lpgsum_col%d_Ohm.png'%(outpath,mce_col))	
		except:
			continue	


if __name__=='__main__':
	main()
	
