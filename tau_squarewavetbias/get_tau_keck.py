import mce_data
import numpy as np
import os,sys
sys.path.insert(0, "../tespy")
import tsdata as ts
import histogram as hist
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import cPickle as pkl
from matplotlib import gridspec


def edge1(t,tau,V,Vc):
	return Vc-V*np.exp(-t/tau)

def edge2(t,tau,V,Vc):
	return V*np.exp(-t/tau) + Vc

def main():

	doplot = False 
	#Lexp = 100000
	#period = 10000
	Nstep = 10
	inpath = '/home/data/cryo/20210309mce1/fasttaurx1_1'
	outpath= inpath.replace('cryo','output')
	if not os.path.isdir(outpath):
		os.makedirs(outpath)

	if not os.path.exists('%s/taus.pkl'%outpath):
		
		#biases = [2500]
		biases  = [5000,4000,3500,3000,2750,2500,2250,2000,1750,1500,1250,1000,500,2]
		biases = np.array(biases)
		nbias = len(biases)
		taus = np.full((2,nbias,33,16),np.nan) # up/down-bias-row-col
	
	
		for row in range(33):
		#for row in [12]:
			for col in range(16):
			#for col in [2]:
					
				print 'row %d, col%d'%(row,col)

				if row in [2,12,22,32]:
					doplot = True
				else:
					doplot = False

				for ib,bias in enumerate(biases):
					fn = '%s/fastdata_row%d_offset0_bias%d'%(inpath,row,bias)
					tsc= ts.get_raw_timestream(fn, calfn='calib_SK')
					y0  = tsc.fb[0,col]

					if all(y0==0) or all(np.isnan(y0)):
						continue

					L = len(y0)
					period = int(L/Nstep)
					Lexp = 	period*Nstep

					fsamp = tsc.info['freq']
					t  = np.array(range(period))*1/fsamp*1000 # t in ms
					t0 = np.array(range(L))*1/fsamp*1000	
					
					y  = np.reshape(y0[:Lexp],(Nstep,period))
					y  = np.transpose(y)
					
					ym1= np.nanmean(y[:,[1,3,5,7,9]],1)
					ym2= np.nanmean(y[:,[0,2,4,6,8]],1)

					ds = 10# ignore the spike
					de =-50# ignore tail

					p1 = np.polyfit(t[int(period*2/3):de], ym1[int(period*2/3):de], 1)
					ym1 = ym1 - p1[0]*t
					p2 = np.polyfit(t[int(period*2/3):de], ym2[int(period*2/3):de], 1)
					ym2 = ym2 - p2[0]*t
					y0 = y0 - (p1[0]+p2[0])*t0/2			
			
					try:	
						popt1,pcov1 = curve_fit(edge1,t[:-ds+de],ym1[ds:de],bounds=([0,0,-8000],[60,3000,8000]))
						if popt1[0]>59 or popt1[1]<1.:
							popt1 = [np.nan, np.nan, np.nan]
						popt2,pcov2 = curve_fit(edge2,t[:-ds+de],ym2[ds:de],bounds=([0,0,-8000],[60,3000,8000]))
						if popt2[0]>59 or popt2[1]<1.:
							popt2 = [np.nan, np.nan, np.nan]
						if popt1[2]<popt2[2]:
							popt1 = [np.nan, np.nan, np.nan]
							popt2 = [np.nan, np.nan, np.nan]
					except:
						popt1 = [np.nan, np.nan, np.nan]
						popt2 = [np.nan, np.nan, np.nan]
					tau1 = popt1[0]
					tau2 = popt2[0]
	
					print '%d : tau1 = %.2f, tau2 = %.2f'%(bias,tau1,tau2)
					taus[0,ib,row,col] = tau1
					taus[1,ib,row,col] = tau2
			
					# plot
					if not doplot:
						continue
	
					fig = plt.figure(figsize=(6, 8)) 
					gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
					ax0 = plt.subplot(gs[0])
					ax0.plot(t0,y0,color='b')
					ax0.set_title('row %d, col %d, bias %d'%(row, col, bias))
					ax0.set_ylabel('sq1 fb')
					ax1 = plt.subplot(gs[1])
					ax1.plot(t[ds:],ym1[ds:],color='b',label='averaged over all cycles')
					ax1.plot(t[ds:],ym2[ds:],color='b')
					if not any(np.isnan(popt1)):
						ax1.plot(t[ds:],edge1(t[:-ds],*popt1),color='r',label='fit')
					if not any(np.isnan(popt2)):
			                	ax1.plot(t[ds:],edge2(t[:-ds],*popt2),color='r')
					ax1.text(0.65,0.52,'tau1 = %.2f ms'%tau1,transform=ax1.transAxes)
					ax1.text(0.65,0.48,'tau2 = %.2f ms'%tau2,transform=ax1.transAxes)
					ax1.set_ylabel('sq1 fb')
					ax1.set_xlabel('time, ms')
					ax1.legend()
			
					plt.tight_layout()
	
					#plt.show()
					plt.savefig('%s/fig_row%d_col%d_bias%d.png'%(outpath,row,col,bias))
	
		fnpickle = '%s/taus.pkl'%(outpath)	
		pkl.dump((biases,taus),open(fnpickle,'w'))			

	else:
		fnpickle = '%s/taus.pkl'%(outpath)	
		print 'Results exist. Load %s.'%fnpickle
		biases,taus = pkl.load(open(fnpickle,'r'))

	# plot summary
	colors = plt.cm.jet(np.linspace(0,1,int(max(biases)+1)))
	for row in range(33):
		for col in range(16):
			fig = plt.figure(figsize=(6, 8))
			gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
			ax1 = plt.subplot(gs[1])
			ax1.set_title('tau vs bias, row %d, col %d'%(row, col))
			taus0 = taus[0,::-1,row,col]
			taus1 = taus[1,::-1,row,col]
			ax1.plot(biases[::-1],taus0,c='k',linestyle='--')
			ax1.scatter(biases[::-1],taus0,s=80,c=colors[biases[::-1]],marker='^',label='rising')
			ax1.plot(biases[::-1],taus1,c='k',linestyle='--')
			ax1.scatter(biases[::-1],taus1,s=80,c=colors[biases[::-1]],marker='v',label='falling')
			ax1.set_xlabel('bias, ADU')
			ax1.set_ylabel('tau, ms')
			ax1.set_xlim(0,6000)
			ax1.legend()
			
			plt.tight_layout()
			
			#plt.show()
			plt.savefig('%s/fig_row%d_col%d.png'%(outpath,row,col))	

	# plot histogram at 0.5 Rn
	fplc = open('/home/data/output/20210219/fasttau2/plc0.txt','r')
	plc = fplc.read()
	fplc.close()
	dplc = pkl.loads(plc)
	rR2000 = np.full((33,16),np.nan) 
	rR3000 = np.full((33,16),np.nan) 
	taus2000 = np.full((2,33,16),np.nan) 
	taus3000 = np.full((2,33,16),np.nan) 
	
	for row in range(33):
		for col in range(16):
			rR2000[row,col] = np.interp(2000,dplc['B'][8,::-1],dplc['Rdet'][row,col,::-1])/dplc['R'][row,col]
			rR3000[row,col] = np.interp(3000,dplc['B'][8,::-1],dplc['Rdet'][row,col,::-1])/dplc['R'][row,col]
			if rR2000[row,col]>0.01 and rR2000[row,col]<0.99:
				taus2000[0,row,col] = np.interp(2000,biases[::-1],taus[0,::-1,row,col])
				taus2000[1,row,col] = np.interp(2000,biases[::-1],taus[1,::-1,row,col])
			if rR3000[row,col]>0.01 and rR3000[row,col]<0.99:
				taus3000[0,row,col] = np.interp(3000,biases[::-1],taus[0,::-1,row,col])
				taus3000[1,row,col] = np.interp(3000,biases[::-1],taus[1,::-1,row,col])
	taus2000 = np.nanmean(taus2000,0)
	taus3000 = np.nanmean(taus3000,0)

	hist.plot_1Dhist((taus2000[:, 8:10]).reshape(2*33), outpath, '/hist_taus_bias2000',
			nbins = 20,
                        maintitle='tau(in transition), bias=2000',
                        xlabel='tau',
                        xunit='ms',
                        binrange=[0., 40.])

	hist.plot_1Dhist((taus3000[:, 8:10]).reshape(2*33), outpath, '/hist_taus_bias3000',
			nbins = 20,
                        maintitle='tau(in transition), bias=3000',
                        xlabel='tau',
                        xunit='ms',
                        binrange=[0., 40.])

	hist.plot_1Dhist((taus3000[:, 8:10]).reshape(2*33), outpath, '/hist_taus_bias3000_zoom',
			nbins = 20,
                        maintitle='tau(in transition), bias=3000, zoomed',
                        xlabel='tau',
                        xunit='ms',
                        binrange=[0., 10.])

	fig = plt.figure(figsize=(6, 6), dpi=80)
	ax1 = plt.subplot(1,1,1)
	ax1.set_title('tau vs R/Rn')
	ax1.scatter(rR2000[:,8:10].reshape(2*33),taus2000[:, 8:10].reshape(2*33),s=80,c='b',marker='^',label='bias=2000')
	ax1.scatter(rR3000[:,8:10].reshape(2*33),taus3000[:, 8:10].reshape(2*33),s=80,c='r',marker='^',label='bias=3000')
	ax1.set_xlabel('R/Rn',fontsize=15)
	ax1.set_ylabel('tau, ms',fontsize=15)
	ax1.set_xlim(0,0.6)
	ax1.grid()
	ax1.legend(fontsize=14)
	plt.tick_params(labelsize=14)
	
	plt.savefig('%s/scatter_taus.png'%(outpath))
	


if __name__=='__main__':
	main()
