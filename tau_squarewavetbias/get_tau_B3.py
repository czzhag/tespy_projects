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
import scipy.io as sio

def edge1(t,tau,V,Vc):
	return Vc-V*np.exp(-t/tau)

def edge2(t,tau,V,Vc):
	return V*np.exp(-t/tau) + Vc

def main():

	doplot = True
	current_biases = [[420,420,420,420,420,420,290,290,290,290,290,290,350,350,350,350,350,350,420,420,420,420,420,420,330,330,330,330,330,330],
			  [400,400,400,  0,400,400,300,300,300,  0,300,300,330,  0,330,330,330,330,330,330,330,330,330,330,350,315,350,315,350,350],
			  [350,350,350,350,350,350,380,380,380,380,380,380,315,315,315,315,315,315,325,325,325,325,325,325,340,340,340,340,340,340],
			  [340,340,340,340,340,  0,  0,350,350,350,350,350,290,290,290,290,290,290,360,360,360,360,360,360,400,400,400,400,  0,  0],
			] 
	for mce in [0,1,2,3]:
		inpath = '/home/data/cryo/20210324_b3_mce%d/fasttau_b3_mce%d/'%(mce,mce)
		outpath= inpath.replace('cryo','output')
		if not os.path.isdir(outpath):
			os.makedirs(outpath)
		current_bias = current_biases[mce]
	
		if not os.path.exists('%s/taus.pkl'%(outpath)):
			
			biases  = [800,700,600,500,450,400,360,320,300,280,260,240,220,200,180,160,140,120,100,2]
			#biases  = [800,500,400,320,280,240,200,140,100,2]
			#biases  = [800,500,400]
			biases = np.array(biases)
			nbias = len(biases)
			taus = np.full((2,nbias,22,30),np.nan) # up/down-bias-row-col
		
		
			for row in range(22):
			#for row in [10]:
				for col in range(30):
				#for col in [16]:
		
					print 'row %d, col%d'%(row,col)
					offset = int(col/10.)*10
	
					for ib,bias in enumerate(biases):
						fn = '%s/fastdata_row%d_rcs_offset%d_bias%d'%(inpath,row,offset,bias)
						if not os.path.exists(fn):
							print 'Missing '+fn
							continue
						try:
							tsc= ts.get_raw_timestream(fn, calfn='calib_BA1')
						except:
							continue
						y  = tsc.fb[0,col]
				
						if len(y)<101000:
							print len(y)
							print 'Skip: Short Time Stream.'
							continue
						y0  = y[20:101020]
						fsamp = tsc.info['freq']
						t  = np.array(range(10100))*1/fsamp*1000 # t in ms
						t0 = np.array(range(101000))*1/fsamp*1000	
	
						fitind1 = []
						fitind2 = []
						for ii in range(5):
							fitind1 = fitind1 + range(ii*20200+100, ii*20200+10000)
							fitind2 = fitind2 + range(ii*20200+10200, ii*20200+20100)
							
						
						p3 = np.poly1d(np.polyfit(t0[fitind2],y0[fitind2],3))
						y0 = y0 - p3(t0)
						y  = np.reshape(y0,(10,10100))
						y  = np.transpose(y)
						
						ym2= np.nanmean(y[:,[1,3,5,7,9]],1)
						ym1= np.nanmean(y[:,[2,4,6,8]],1)
				
				
						ds = 10 # ignore the spike
						try:	
							popt1,pcov1 = curve_fit(edge1,t[:-ds],ym1[ds:],bounds=([0,0,-2000],[100,2000,2000]))
							if popt1[0]>=100 or popt1[1]<2.:
								popt1 = [np.nan, np.nan, np.nan]
							popt2,pcov2 = curve_fit(edge2,t[:-ds],ym2[ds:],bounds=([0,0,-2000],[100,2000,2000]))
							if popt2[0]>=100 or popt2[1]<2.:
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
						if not row==2:
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
			pkl.dump((biases,taus),open(fnpickle,'w'))			
	
		else:
			fnpickle = '%s/taus.pkl'%(outpath)	
			print 'Results exist. Load %s.'%fnpickle
			biases,taus = pkl.load(open(fnpickle,'r'))
	
		# plot summary
		colors = plt.cm.jet(np.linspace(0,1,int(max(biases)+1)))
		current_taus = np.full((22,30),np.nan)
		for row in range(22):
		#for row in [10]:
			for col in range(30):
			#for col in [16]:
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
				current_tau0 = np.interp(current_bias[col],biases[::-1],taus0)
				current_tau1 = np.interp(current_bias[col],biases[::-1],taus1)
				current_taus[row,col] = np.nanmean([current_tau0,current_tau1])
				ax1.axvline(x=current_bias[col],color='r',linestyle='--')
				ax1.set_xlabel('bias, ADU')
				ax1.set_ylabel('tau, ms')
				ax1.set_xlim(0,900)
				ax1.legend()
				plt.tight_layout()
				plt.savefig('%s/fig_row%d_col%d.png'%(outpath,row,col))	
	
		fnpickle = '%s/taus_atcmbbias.pkl'%(outpath)	
		pkl.dump((current_bias,current_taus),open(fnpickle,'w'))
		fnmat    = '%s/taus_mce%d.mat'%(outpath,mce)
		sio.savemat(fnmat,{"tau":current_taus, "bias":current_bias})
				

	'''
	# plot histogram at 0.5 Rn
	fplc = open('/home/data/output/20210219/fasttau2/plc0.txt','r')
	plc = fplc.read()
	fplc.close()
	dplc = pkl.loads(plc)
	rR2000 = np.full((22,16),np.nan) 
	rR3000 = np.full((22,16),np.nan) 
	taus2000 = np.full((2,22,16),np.nan) 
	taus3000 = np.full((2,22,16),np.nan) 
	
	for row in range(22):
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

	hist.plot_1Dhist((taus2000[:, 8:10]).reshape(2*22), outpath, '/hist_taus_bias2000',
			nbins = 20,
                        maintitle='tau(in transition), bias=2000',
                        xlabel='tau',
                        xunit='ms',
                        binrange=[0., 40.])

	hist.plot_1Dhist((taus3000[:, 8:10]).reshape(2*22), outpath, '/hist_taus_bias3000',
			nbins = 20,
                        maintitle='tau(in transition), bias=3000',
                        xlabel='tau',
                        xunit='ms',
                        binrange=[0., 40.])

	hist.plot_1Dhist((taus3000[:, 8:10]).reshape(2*22), outpath, '/hist_taus_bias3000_zoom',
			nbins = 20,
                        maintitle='tau(in transition), bias=3000, zoomed',
                        xlabel='tau',
                        xunit='ms',
                        binrange=[0., 10.])

	fig = plt.figure(figsize=(6, 6), dpi=80)
	ax1 = plt.subplot(1,1,1)
	ax1.set_title('tau vs R/Rn')
	ax1.scatter(rR2000[:,8:10].reshape(2*22),taus2000[:, 8:10].reshape(2*22),s=80,c='b',marker='^',label='bias=2000')
	ax1.scatter(rR3000[:,8:10].reshape(2*22),taus3000[:, 8:10].reshape(2*22),s=80,c='r',marker='^',label='bias=3000')
	ax1.set_xlabel('R/Rn',fontsize=15)
	ax1.set_ylabel('tau, ms',fontsize=15)
	ax1.set_xlim(0,0.6)
	ax1.grid()
	ax1.legend(fontsize=14)
	plt.tick_params(labelsize=14)
	
	plt.savefig('%s/scatter_taus.png'%(outpath))
	
	'''

if __name__=='__main__':
	main()
