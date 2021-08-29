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
	inpath = '/home/data/cryo/20210302/20210302b3'
	outpath= inpath.replace('cryo','output') + '/test2/'
	if not os.path.isdir(outpath):
		os.makedirs(outpath)
	flist = glob.glob('%s/*plc*'%inpath)

	# setup
	import calib_SK as calib_sk
	calib = calib_sk.calib_sk()
	fitrange = {}
	fitrange['rnti_low'] = 1000
	fitrange['rnti_hgh'] = 1200

	for lcfn in flist:
		plc = lc.get_PLC(lcfn,fitrange,calib,-1)
		bias= min(plc.B[8])
		#if not bias==0:
		#	continue
	
		fnpickle = '%s/plc%d.txt'%(outpath,int(bias))
		f = open(fnpickle,'w')
		f.write(pkl.dumps(plc.__dict__))
		f.close()
		print 'Save plc to %s'%fnpickle
	
		if not doplot:
			continue

		outpathsub= outpath + '/plc_bias%d'%bias
		if not os.path.isdir(outpathsub):
			os.makedirs(outpathsub)

		for row in range(22):
			for col in range(30):
				fig = plt.figure(figsize=(6, 8))
				ax1 = plt.subplot(2,1,1)
				#ax1.plot(plc.B[col],plc.Idet[row,col]*1e6,color='b')
				ax1.plot(plc.B[col],plc.Y[row,col],color='b')
				ax1.grid()
				#ax1.text(0.7,0.9 ,'Rn = %.3f mOhm'%plc.R[row,col],transform=ax.transAxes)
				#ax1.text(0.7,0.85,'Rop= %.3f mOhm'%plc.R0[row,col],transform=ax.transAxes)
				ax1.set_ylabel(r'$I_{tes}, \mu$A')
				ax1.set_xlabel('bias, ADU')
				ax1.set_title('PLC row %d, col %d, bias %d'%(row, col, bias))
				#ax1.set_xlim(0,10000)
				#ax1.set_ylim(0,40)

				ax2 = plt.subplot(2,1,2)
				ax2.plot(plc.Rdet[row,col]*1e3,plc.Pdet[row,col]*1e12,color='b')
				ax2.grid()
				ax2.text(0.7,0.9 ,'Rn = %.3f mOhm'%plc.R[row,col],transform=ax2.transAxes)
				ax2.text(0.7,0.85,'Rop= %.3f mOhm'%plc.R0[row,col],transform=ax2.transAxes)
				ax2.set_ylabel('P, pW')
				ax2.set_xlabel(r'R, m$\Omega$')
				#ax2.set_xlim(0,300)
				#ax2.set_ylim(0,25)

				plt.savefig('%s/row%d_col%d_bias%d.png'%(outpathsub,row,col,bias))

				# make sum
				if bias==0:
					biases = [6000,5600,5200,4800,4400,4000,3800,3600,3400,3200,3000,2800,2600,2400,2200,2000,1800,1600,1500,1400,1300,1200,1100,1000,0]
					biases = np.array(biases)
					colors = plt.cm.jet(np.linspace(0,1,int(max(biases)+1)))
					for b in biases:
						if b>=2800:
							ax1.axvline(x=b,color=colors[b],linestyle='--',label='bias %d'%b)
						else:
							ax1.axvline(x=b,color=colors[b],linestyle='--')
						r = np.interp(b,plc.B[col,::-1],plc.Rdet[row,col,::-1]*1e3)
						if b>=2800:
							ax2.axvline(x=r,color=colors[b],linestyle='--')
						else:
							ax2.axvline(x=r,color=colors[b],linestyle='--',label='bias %d'%b)
						ax1.legend(loc=1)
						ax2.legend(loc=1)
					plt.savefig('%s/row%d_col%d_sum.png'%(outpathsub,row,col))

if __name__=='__main__':
	main()

