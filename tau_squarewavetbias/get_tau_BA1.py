import mce_data
import numpy as np
import os,sys
sys.path.insert(0, "../tespy")
import tsdata as ts
import histogram as hist
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter as sgfltr
import cPickle as pkl
from matplotlib import gridspec


def edge1(t,tau,V,Vc): # going up
	return Vc-V*np.exp(-t/tau)

def edge2(t,tau,V,Vc): # going down
	return V*np.exp(-t/tau) + Vc

def main():

	doplot = True 
	dostats= False
	Nstep = 10
	inpath = '/home/data/cryo/20210309/fasttau10'
	outpath= inpath.replace('cryo','output')
	if not os.path.isdir(outpath):
		os.makedirs(outpath)

	if not os.path.exists('%s/taus.pkl'%outpath):
		
		#biases  = [6000,5600,5200,4800,4400,4000,3800,3600,3400,3200,3000,2800,2600,2400,2200,2000,1800,1600,1500,1400,1300,1200,1100,1000,0]
		biases  = [750,500,450,350,280,240,200,180,160,120,2]
		biases = np.array(biases)
		nbias = len(biases)
		taus = np.full((2,nbias,33,24),np.nan) # up/down-bias-row-col
		for row in range(33):
		#for row in [25]:
			for col in range(16):
	
				print 'row %d, col%d'%(row,col)
				rc     = int(col/8.) + 1

				for ib,bias in enumerate(biases):
					fn = '%s/fastdata_row%d_rc%d_offset_bias%d'%(inpath,row,rc,bias)
					tsc= ts.get_raw_timestream(fn, calfn='calib_BA1')
					y0 = tsc.fb[0,col-rc*8]
					y0s = sgfltr(y0,101,5)
					
					if all(y0==0) or all(np.isnan(y0)):
						continue

					L = len(y0)
					period = int(0.4*tsc.info['freq'])
					y = np.full((Nstep-1,period), np.nan)

					for i in range(Nstep-2):
						#s = np.argmin(y0s[0:period]) + i*period
						s = 3950 + i*period
						y[i,:period-100] = y0[s:s+period-100]
					y  = np.transpose(y)
					
					fsamp = tsc.info['freq']
					t  = np.array(range(period))*1/fsamp*1000 # t in ms
					t0 = np.array(range(L))*1/fsamp*1000	

					ym1= np.nanmean(y[:,[0,2,4,6]],1)
					ym2= np.nanmean(y[:,[1,3,5,7]],1)
			
			
					ds = 10 # ignore the spike
					try:	
					#if 1:
						popt1,pcov1 = curve_fit(edge1,t[:-ds-100],ym1[ds:-100],bounds=([0,0,-5000],[400,5000,5000]))
						print(popt1)
						if popt1[0]>200 or popt1[1]<5.:
							popt1 = [np.nan, np.nan, np.nan]
						popt2,pcov2 = curve_fit(edge2,t[:-ds-100],ym2[ds:-100],bounds=([0,0,-5000],[400,5000,5000]))
						print(popt2)
						if popt2[0]>200 or popt2[1]<5.:
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
					ax0.plot(t0[:],(y0),color='b')
					ax0.set_title('row %d, col %d, bias %d'%(row, col, bias))
					ax0.set_ylabel('sq1 fb')
					ax1 = plt.subplot(gs[1])
					ax1.plot(t[ds:],ym1[ds:],color='b',label='averaged over all cycles')
					ax1.plot(t[ds:],ym2[ds:],color='b')
					ax1.plot(t[ds:],edge1(t[:-ds],*popt1),color='r',label='fit')
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
		#pkl.dump((biases,taus),open(fnpickle,'w'))			

	else:
		fnpickle = '%s/taus.pkl'%(outpath)	
		print 'Results exist. Load %s.'%fnpickle
		biases,taus = pkl.load(open(fnpickle,'r'))

	# plot summary
	colors = plt.cm.jet(np.linspace(0,1,int(max(biases)+1)))
	for row in range(33):
		for col in range(24):
			fig = plt.figure(figsize=(6, 8))
			gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
			ax1 = plt.subplot(gs[1])
			ax1.set_title('tau vs bias, row %d, col %d'%(row, col))
			taus0_ = taus[0,::-1,row,col]
			taus0 = taus0_[np.where(~np.isnan(taus0_))]
			bias0 = biases[::-1][np.where(~np.isnan(taus0_))]
			taus1_ = taus[1,::-1,row,col]
			taus1 = taus1_[np.where(~np.isnan(taus1_))]
			bias1 = biases[::-1][np.where(~np.isnan(taus1_))]
			ax1.plot(bias0,taus0,c='k',linestyle='--')
			ax1.scatter(biases[::-1],taus0,s=80,c=colors[biases[::-1]],marker='^',label='rising')
			ax1.plot(bias1,taus1,c='k',linestyle='--')
			ax1.scatter(biases[::-1],taus1,s=80,c=colors[biases[::-1]],marker='v',label='falling')
			ax1.set_xlabel('bias, ADU')
			ax1.set_ylabel('tau, ms')
			ax1.set_xlim(0,1000)
			ax1.legend()
			
			plt.tight_layout()
			
			#plt.show()
			plt.savefig('%s/fig_row%d_col%d.png'%(outpath,row,col))	

	if dostats:
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
			for col in [8,9]:
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
