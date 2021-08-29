import os,sys
sys.path.insert(0,'../tespy')
import numpy as np
import cPickle as pickle
import histogram as hist
import shutil
import easyplots as esp
import ba30_ModuleMapping_N3 as bm

# inputs
run = '20210423_SK_BA30N6N7'
inpath = '/home/data/output/20210423/%s/'%run
fn = '%s_dpdt_rnti.pkl'%run
inpath_al = '/home/data/output/20210423/%s_Al/'%run
fn_al = '%s_Al_dpdt_rnti.pkl'%run
ims=[0,4]
modules=['N6','N7']
cmap = 'jet'

infile = inpath + fn
d = pickle.load(open(infile,'r'))
infile_al = inpath_al + fn_al
d_al = pickle.load(open(infile_al,'r'))

# output
outpath = inpath
if not os.path.isdir(outpath):
        os.makedirs(outpath)
shutil.copy2(os.path.realpath(__file__), outpath + (os.path.realpath(__file__).split("/")[-1]).replace(".py",".txt"))


rntil = np.full((24, 33), float('nan'))
dpdtl = np.full((24, 33), float('nan'))
rnall = np.full((24, 33), float('nan'))
rntid = np.full((24, 33), float('nan'))
dpdtd = np.full((24, 33), float('nan'))
rnald = np.full((24, 33), float('nan'))

for col in range(24):
	for row in range(33):
		im,dc,dr,pol=bm.mce2det(col,row)
		if pol=='D':
			try:
				rntid[col, row] = d[1][col][row]
				dpdtd[col, row] = d[0][col][row]
				rnald[col, row] = d_al[1][col][row]-rntid[col, row]
			except:
				continue
		elif pol in ['A','B']:
			try:
				rntil[col, row] = d[1][col][row]
				dpdtl[col, row] = d[0][col][row]
				rnall[col, row] = d_al[1][col][row]-rntil[col, row]
			except:
				continue
		else:
			continue

for jj,im in enumerate(ims):
	module=modules[jj]
	print (rntid[im*2:im*2+2, :]).reshape(2*33) 
	hist.plot_1Dhist2((rntil[im*2:im*2+2, :]).reshape(2*33), (rntid[im*2:im*2+2, :]).reshape(2*33),
			outpath, '%s_rnti_hist'%module,
			maintitle=r'$R_n(Ti)$, %s'%module,
			xlabel=r'$R_n(Ti)$', 
			xunit=r'm$\Omega$',
			binrange=[50., 180.])
	
	hist.plot_1Dhist2((rnall[im*2:im*2+2, :]).reshape(2*33), (rnald[im*2:im*2+2, :]).reshape(2*33),
			outpath, '%s_rnal_hist'%module,
			maintitle=r'$R_n(Al)$, %s'%module,
			xlabel=r'$R_n(Al)$', 
			xunit=r'm$\Omega$',
			binrange=[0., 1000.])
	
	hist.plot_1Dhist2((dpdtl[im*2:im*2+2, :]).reshape(2*33), (dpdtd[im*2:im*2+2, :]).reshape(2*33),
			outpath, '%s_dpdt_hist'%module,
	                maintitle=r'dP/dT, %s'%module,
	                xlabel=r'dP/dT', 
	                xunit=r'pW/K',
			nbins=20,
	                binrange=[.0, 0.075])
