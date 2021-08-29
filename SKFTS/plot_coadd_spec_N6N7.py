import os,sys
sys.path.insert(0,'../tespy')
import numpy as np
import cPickle as pickle
import shutil
import matplotlib.pyplot as plt
import pylab as pl

# inputs
modules = ['N6','N7']
ims=[0,4]
inpath='/home/data/output/20210511'
fn='SK_N7B6_fts.pkl'
lowp=0
highp=50

infile='%s/%s'%(inpath,fn)
d=pickle.load(open(infile,'r'))
interfg=d[6]
spcos=d[7]
spsin=d[8]
v_mirror=d[9]
nSample=d[10]

c = 299792458.
delta_step = 2*v_mirror/nSample
fNyq = c/2/delta_step/1e9 #in GHz

# output
outfigdir=inpath + '/coadd_spec'
if not os.path.isdir(outfigdir):
	os.makedirs(outfigdir)
shutil.copy2(os.path.realpath(__file__), outfigdir + (os.path.realpath(__file__).split("/")[-1]).replace(".py",".txt"))

detlist = spcos.keys()
nlen_hlf_intf = len(spcos[detlist[0]])
freq_ = np.linspace(1,fNyq,nlen_hlf_intf+1)
freq = freq_[:-1]

# For 30/40, we practically can only have upto 2 modules in SK
# But considering read out ability, we can have up to 8 modules. 16(cols)/2(col-per-module) = 8 (modules)
css = np.full((8,len(detlist),nlen_hlf_intf),np.nan) 
sss = np.full((8,len(detlist),nlen_hlf_intf),np.nan) 

for idet,det in enumerate(detlist):
	col = int(det.split(('c'))[1])	
	css[int(col/2),idet] = spcos[det]
	sss[int(col/2),idet] = spsin[det]

for im in ims:
	if im == 0:
		module = modules[0]
	elif im == 4:
		module = modules[1]
	else:
		module = ''
		print 'WARN: Wrong module!!'
	css_ = css[im]
	sss_ = sss[im]
	cs = np.nanmean(css_,(0))
	ss = np.nanmean(sss_,(0))

	pl.figure(figsize=(18,15), dpi=80)
	for idet in range(len(detlist)):
		if np.all(np.isnan(css_[idet])):
			continue
		plt.plot(freq,css_[idet],'grey')
	plt.plot(freq,ss,'g',label='sin')
	plt.plot(freq,cs,'r',label='cos')
	plt.xlim(lowp,highp)
	plt.ylim(-0.05,1)
	plt.xlabel('freq GHz',fontsize=25)
	plt.ylabel('efficiency',fontsize=25)
	plt.title('%s, coadded spec'%module,fontsize=25)
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)

	fn = '%s/%s_coadded_spec_freqrange%d.png'%(outfigdir,module,highp)
	pl.savefig(fn)
