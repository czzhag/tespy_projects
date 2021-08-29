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

print('loading...')
info=sio.loadmat('/home/czhang/analysispy/ba_fpdata.mat')
d=sio.loadmat('/home/czhang/pipeline/202001magBA/reduc_data/magBA_py.mat')
ampfit = pkl.load(open('MagAmps.pkl','r')) 
indl=info['ind']['l'][0][0][0]
indd=info['ind']['d'][0][0][0]
indo=info['ind']['o'][0][0][0]
out_path = './outfig/magfit/'
if not os.path.isdir(out_path):
    os.makedirs(out_path)

a0s=ampfit[0]
a1s=ampfit[1]
a2s=ampfit[2]
a3s=ampfit[3]
a4s=ampfit[4]
a5s=ampfit[5]

def cosmodel(x,a0,a1,a2,a3,a4,a5):
    return a0+a1*x+a2*np.cos(2*np.pi/2080*x+a3)+a4*np.cos(2*np.pi/1040*x+a5)

for ib in range(len(d['bl'][0])):
    print('do bias %d'%ib)
    for idet in range(792):
        azd=d['azdata'][ib,300:,idet]
        ind=np.where(~np.isnan(azd))[0]
        azd=azd[ind]
        x  =np.linspace(0,len(azd),len(azd))
        popt=[a0s[ib,idet], a1s[ib,idet], a2s[ib,idet], a3s[ib,idet], a4s[ib,idet], a5s[ib,idet]]
        
        plt.close()
        plt.plot((azd-popt[0]-x*popt[1])*1e12, label='timestream')
        plt.plot(x, (cosmodel(x, *popt)-popt[0]-x*popt[1])*1e12, label='fit')
        plt.xlabel('sample', fontsize=14)
        plt.ylabel('pA', fontsize=14)
        if idet+1 in indl:
            plt.title('az scanset, gcp %d, bias %d, %s, detcol %d, detrow %d, pol %s'%(idet, 
                d['bl'][0][ib], 
                modules[int(idet/66)], 
                int(info['p']['det_col'][0][0][idet][0]),
                int(info['p']['det_row'][0][0][idet][0]),
                str(info['p']['pol'][0][0][idet][0][0])))
        elif idet+1 in indd:
            plt.title('az scanset, gcp %d, bias %d, %s, D'%(idet, d['bl'][0][ib], modules[int(idet/66.)]))
        elif idet+1 in indo:
            plt.title('az scanset, gcp %d, bias %d, %s, O'%(idet, d['bl'][0][ib], modules[int(idet/66.)]))
        else:
            continue
        plt.grid()
        plt.savefig(out_path+'gcp%d_bias%d.png'%(idet, d['bl'][0][ib]))
