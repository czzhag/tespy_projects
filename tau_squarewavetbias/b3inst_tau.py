# plot current tau in B3

import numpy as np
import os,sys
import matplotlib.pyplot as plt
sys.path.insert(0, "../tespy")
import histogram as hist
import cPickle as pkl

mce = 0
inpath = '/home/data/output/20210324_b3_mce%d/fasttau_b3_mce%d/'%(mce,mce)
fnpickle = '%s/taus_atcmbbias.pkl'%(inpath)
biases,taus = pkl.load(open(fnpickle,'r'))
taus = (taus.T).reshape(-1)

outpath = '/home/cheng/Pole2019/analysispy/tau_squarewavetbias/b3inst/'
if not os.path.isdir(outpath):
	os.makedirs(outpath)

from get_array_info import get_array_info
p,ind = get_array_info('/home/data/fpdatapy/b3_fpdata_2017.mat')

mask = np.full((22*30),np.nan)
for ii in range(22*30):
	if ii+1 in ind.rgl:
		mask[ii] = 1

hist.plot_1Dhist(taus*mask, outpath, '/hist_currenttaus',
	nbins = 30,
	maintitle=r'$\tau$, cmb bias, mce %d'%mce,
	xlabel=r'$\tau$',
	xunit=r'${\rm ms}$',
	binrange=[0., 5.])



