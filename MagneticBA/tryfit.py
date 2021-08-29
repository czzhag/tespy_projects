import sys,os
sys.path.insert(0, "../tespy")
import numpy as np
import scipy.io as sio
import lcdata
from lcdata import get_LC_calibed as lc
from lcdata import get_PR_Ti as pr
import calib_BA1
calib = calib_BA1.calib_ba1()
import easyplots as esp
from tsdata import remove_baseline as RemoveBaseline
modules=['N2','M4','N3','M5','M6','N1','M8','Mx2','M3','M7','N4','M9']
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit
import cPickle as pkl

info=sio.loadmat('/home/czhang/analysispy/ba_fpdata.mat')
d = pkl.load(open('tryfit.pkl','r'))
azs = d[0]
azd = d[1]

def cosmodel(x,a0,a1,a2,a3,a4,a5):
    return a0+a1*x+a2*np.cos(2*np.pi/2087*x+a3)+a4*np.cos(2*np.pi/1040*x+a5)

popt, pcov = curve_fit(cosmodel, np.linspace(0,len(azd[300:]),len(azd[300:])), azd[300:])
plt.close()
plt.plot(azd[300:])
plt.plot(np.linspace(0,len(azd[300:]),len(azd[300:])), cosmodel(np.linspace(0,len(azd[300:]),len(azd[300:])), *popt))
plt.show()
