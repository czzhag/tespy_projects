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
	#tiles = [1,4,5,9,11,12,13,16]
	tiles = [2, 10, 15]
	p,ind = get_array_info('/home/data/fpdatapy/b3_fpdata_2017.mat')


	current_biases = [[420,420,420,420,420,420,290,290,290,290,290,290,350,350,350,350,350,350,420,420,420,420,420,420,330,330,330,330,330,330],
			  [400,400,400,  0,400,400,300,300,300,  0,300,300,330,  0,330,330,330,330,330,330,330,330,330,330,350,315,350,315,350,350],
			  [350,350,350,350,350,350,380,380,380,380,380,380,315,315,315,315,315,315,325,325,325,325,325,325,340,340,340,340,340,340],
			  [340,340,340,340,340,  0,  0,350,350,350,350,350,290,290,290,290,290,290,360,360,360,360,360,360,400,400,400,400,  0,  0],
			] 
	bl = [800,757,715,673,631,589,547,505,
		480,463,440,421,400,378,360,336,
		320,310,300,290,280,270,260,250,
		240,230,220,210,200,190,180,170,
		160,150,140,120,100, 60, 42,  0]

	for tile in tiles:
		mce=np.unique(p.mce[np.where(p.tile==tile)][0])[0]
		if mce in [0,1]:
			lcfn   = 'LC_Ti_MCE%d_sky_330mK_r02'%mce
		elif mce==2:
			lcfn   = 'LC_Ti_MCE%d_sky_328mK_r02'%mce
		elif mce==3:
			lcfn   = 'LC_Ti_MCE%d_sky_329mK_r02'%mce
		else:
			print "Error: Wrong MCE %20d"%mce
			return
			
		lcpath = inpath + '/' + lcfn
		outpath = '/home/data/output/20210413/loopgain_tau_mag/tile%d'%tile
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
	
		fnpickle = '%s/Lpg.pkl'%(outpath) # add _ at the end for plotting
		if os.path.exists(fnpickle):
			print 'Results exist. Load %s.'%fnpickle
			Lpg,lcpath,fitrange,calib,vnames = pkl.load(open(fnpickle,'r'))
		else:
			Lpg = lc.get_loopgain(plc.Pdet,plc.Rdet,win=40)
			pkl.dump((Lpg,lcpath,fitrange,calib,'Lpg,lcpath,fitrange,calib'),open(fnpickle,'w'))

		fnpickletau = '/home/data/output/20210324_b3_mce%d/fasttau_b3_mce%d/taus.pkl'%(mce,mce)
		biases,taus = pkl.load(open(fnpickletau,'r'))
		tau0 = 100 #ms

		dmat=loadmat('/home/data/output/201911magB3_reduc_data/magB3_di_sum.mat')
		
		mce_cols=np.unique(p.mce_col[np.where(p.tile==tile)])

		for col in mce_cols:
			
			bias0 = current_biases[mce][col]

			for row in range(1,22):
				if not doplot:
					continue
				#if all(np.isnan(plc.Pdet[row,col,:xtail])):
				#	continue	

				j = np.where((p.mce==mce)*(p.mce_col==col)*(p.mce_row==row))[0][0]
				dtype = p.type[j][0][0].encode()
				if not (dtype in ['L','D']):
					continue
				dcol = int(p.det_col[j][0])
				drow = int(p.det_row[j][0])
				pol  = p.pol[j][0][0].encode()
	
				# di/tau/loop gain	
				fig = plt.figure(figsize=(6, 8))
				gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
				ax0 = plt.subplot(gs[0])
				ax0.plot(plc.Rdet[row,col,:xtail]/plc.R[row,col],plc.Pdet[row,col,:xtail]*1e12,color='b')
				rR0 = np.interp(bias0,bl[::-1],(dmat['R_blkon'][::-1,j]/dmat['Rn'][0,j]))
				ax0.axvline(x=rR0, color='gray', linestyle='-.')
				ax0.set_xlim(0,1)
				ax0.set_ylabel('Psat, pW')
				ax0.set_xlabel(' ')
				ax0.set_title('P/loop gain vs R/Rn, row %d, col%d'%(row,col))
				ax0.set_ylim(0,10)

				ax1 = plt.subplot(gs[1])
				Rdet1 = np.where((plc.Pdet[row,col,:]>=0.5e-12), plc.Rdet[row,col,:], np.inf)
				Lpg[row,col,:]= np.where((plc.Pdet[row,col,:]>=(plc.Pdet[row,col,np.argmin(Rdet1)]*0.8)), Lpg[row,col,:], np.nan)
				Lpg1 = abs(Lpg[row,col,:xtail])
				rR1  = plc.Rdet[row,col,:xtail]/plc.R[row,col]
				ax1.plot(rR1, Lpg1,color='b', label='loop gain from P-R')

				Rtaus=np.interp(biases[::-1],plc.B[col,::-1],plc.Rdet[row,col,::-1])
				ax1.scatter(Rtaus/plc.R[row,col],tau0/taus[0,::-1,row,col], s=80,c='r', marker='x',label='tau0(100 ms)/tau')	
				ax1.set_ylim(1,1e4)
				ax1.set_xlim(0,1)
				ax1.set_yscale('log')
				ax1.set_ylabel('Loop gain')
				ax1.set_xlabel('R/Rn')
				plt.grid()	
				plt.legend(loc=2)
	
				ax2=ax1.twinx()
				#ax2.scatter(dmat['R_blkon'][:,j]/dmat['Rn'][0,j], dmat['di_az_blkon'][:,j]*1e6, s=80, c='k', marker='^', label='dI - azscan withplate on')
				ax2.plot(dmat['R_blkon'][::-1,j]/dmat['Rn'][0,j], dmat['di_az_blkon'][::-1,j]*1e6, 'xk--', label='dI - azscan withplate on') 
				ax2.axvline(x=rR0, color='gray', linestyle='-.', label='CMB bias')
				ax2.set_ylim(0, 0.12)
				ax2.set_ylabel('magnetic response dI, uA')

				plt.legend()
				plt.tight_layout()

				print '%s/lpgtaudi_row%d_col%d.png'%(outpath,row,col)		
				plt.savefig('%s/lpgtaudi_row%d_col%d.png'%(outpath,row,col))	
				plt.savefig('%s/lpgtaudi_tile%d_drow%d_dcol%d_%s.png'%(outpath,tile,drow,dcol,pol))	

				# dp/tau/loop gain	
				fig = plt.figure(figsize=(6, 8))
				gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
				ax0 = plt.subplot(gs[0])
				ax0.plot(plc.Rdet[row,col,:xtail]/plc.R[row,col],plc.Pdet[row,col,:xtail]*1e12,color='b')
				ax0.set_xlim(0,1)
				ax0.set_ylabel('Psat, pW')
				ax0.set_xlabel(' ')
				ax0.set_title('P/loop gain vs R/Rn, row %d, col%d'%(row,col))
				ax0.set_ylim(0,10)
				plt.legend()

				ax1 = plt.subplot(gs[1])
				Rdet1 = np.where((plc.Pdet[row,col,:]>=0.5e-12), plc.Rdet[row,col,:], np.inf)
				Lpg[row,col,:]= np.where((plc.Pdet[row,col,:]>=(plc.Pdet[row,col,np.argmin(Rdet1)]*0.8)), Lpg[row,col,:], np.nan)
				Lpg1 = abs(Lpg[row,col,:xtail])
				rR1  = plc.Rdet[row,col,:xtail]/plc.R[row,col]
				ax1.plot(rR1, Lpg1,color='b', label='loop gain from P-R')

				Rtaus=np.interp(biases[::-1],plc.B[col,::-1],plc.Rdet[row,col,::-1])
				ax1.scatter(Rtaus/plc.R[row,col],tau0/taus[0,::-1,row,col], s=80,c='r', marker='x',label='tau0(100 ms)/tau')	
				ax1.set_ylim(1,1e4)
				ax1.set_xlim(0,1)
				ax1.set_yscale('log')
				ax1.set_ylabel('Loop gain')
				ax1.set_xlabel('R/Rn')
				plt.grid()	
				plt.legend(loc=2)
	
				ax2=ax1.twinx()
				#ax2.scatter(dmat['R_blkon'][:,j]/dmat['Rn'][0,j], dmat['di_az_blkon'][:,j]*dmat['R_blkon'][:,j]*dmat['I_blkon'][:,j]*1e12, s=80, c='k', marker='^', label='dP - azscan withplate on')
				ax2.plot(dmat['R_blkon'][::-1,j]/dmat['Rn'][0,j], dmat['di_az_blkon'][::-1,j]*dmat['R_blkon'][:,j]*dmat['I_blkon'][:,j]*1e12, 'xk--', label='dP - azscan withplate on') 
				ax2.axvline(x=rR0, color='gray', linestyle='-.', label='CMB bias')
				ax2.set_ylim(0, 0.1)
				ax2.set_ylabel('magnetic response dP, pW')

				plt.legend()
				plt.tight_layout()

				print '%s/lpgtaudp_row%d_col%d.png'%(outpath,row,col)		
				plt.savefig('%s/lpgtaudp_row%d_col%d.png'%(outpath,row,col))	
				plt.savefig('%s/lpgtaudp_tile%d_drow%d_dcol%d_%s.png'%(outpath,tile,drow,dcol,pol))	


if __name__=='__main__':
	main()
	
