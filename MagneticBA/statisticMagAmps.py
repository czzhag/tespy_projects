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
import histogram as hist

NETselectbias=[80,80,120,160,60,80,160,80,120,100,80,120,220,200,160,120,100,140,160,160,120,240,300,240]
MAGselectbias=[200,180,200,200,240,220,260,160,180,180,300,300,440,400,380,300,240,240,320,320,300,300,640,640]
print('loading...')
info=sio.loadmat('/home/czhang/analysispy/ba_fpdata.mat')
d=sio.loadmat('/home/czhang/pipeline/202001magBA/reduc_data/magBA_py.mat', variable_names='bl')
bl=d['bl'][0]
ampfit = pkl.load(open('MagAmps.pkl','r'))
indl=info['ind']['l'][0][0][0]
indd=info['ind']['d'][0][0][0]
indo=info['ind']['o'][0][0][0]
out_path = './outfig/magstat/'
if not os.path.isdir(out_path):
    os.makedirs(out_path)

#D-39*792
a0s=ampfit[0]
a1s=ampfit[1]
a2s=ampfit[2]
a3s=ampfit[3]
a4s=ampfit[4]
a5s=ampfit[5]

#elnod ref
delnod=sio.loadmat('/home/czhang/pipeline/202001magBA/reduc_data/enBA.mat')
elnodgain=delnod['en']['g'][0][0] #len(bl)bias*2xxxdets*4(offset, gain, gof, ?)
elnodam  =delnod['endatx'] #len(bl)*ndets*1200
elnodfb  =delnod['endaty']
elnodel  =delnod['endatel']
bll = np.full((len(bl),33), np.nan)

Acos=np.full((24,len(bl),33), np.nan)
Acos2=np.full((24,len(bl),33), np.nan)
Aelnod=np.full((24,len(bl),33), np.nan)

for col in range(24):
    Acos[col]=2*a2s[:,col*33:(col+1)*33] #len(bl) by 33
    Acos[col,np.where(bl==40)[0], :]=[np.nan]*33
    Acos2[col]=2*a4s[:,col*33:(col+1)*33] #len(bl) by 33
    Acos2[col,np.where(bl==40)[0], :]=[np.nan]*33
    Aelnod[col]=np.full((len(bl),33),np.nan)
    for row in range(33):
        # mask out not responsive channels. condition is max elnod gain above 1000 
        if max(elnodgain[:,col*33+row, 1])<1000:
            Acos[col,:,row] = [np.nan]*len(bl)
            Acos2[col,:,row] = [np.nan]*len(bl)
        else:
            # get elnod amplitude for working channels
            for ib in range(len(bl)):
                #Aelnod[col,ib,row]=elnodgain[ib,col*33+row,1]*(np.nanmax(elnodam[ib,col*33+row,:])-np.nanmin(elnodam[ib,col*33+row,:]))*calib['FB_CAL'][0]
                Aelnod[col,ib,row]=(np.nanmax(elnodfb[ib,col*33+row,:])-np.nanmin(elnodfb[ib,col*33+row,:]))*calib['FB_CAL'][0]
                if np.isnan(Aelnod[col,ib,row]) or np.isinf(Aelnod[col,ib,row]):
                    Acos[col,ib,row] = np.nan
                    Acos2[col,ib,row]= np.nan
    plt.close()
    idet = np.array(range(33))
    xi, yi = np.mgrid[bl.min():bl.max():len(bl)*1j, idet.min():idet.max():33*1j]
    plt.pcolormesh(xi.T, yi.T, np.log10(abs(Aelnod[col,::-1,:])*1e12).T, vmin=1., vmax=5., cmap="coolwarm")
    cbar = plt.colorbar()
    cbar.set_label(r'$log10(elnod p2p, pA)$')
    #for ic in range(2):
    #    plt.axvline(x=33*ic, color='k', linestyle='--')
    #    plt.axhline(y=33*ic, color='k', linestyle='--')
    #plt.text(60, 60, modules[im], fontsize=14)
    plt.axvline(x=NETselectbias[col], color='k', linestyle='--',label='bias selected by NET script')
    plt.axvline(x=MAGselectbias[col], color='k', linestyle='-',label='bias adjusted by magnetic analysis')
    plt.xlabel("bias, ADU", fontsize=14)
    plt.legend()
    plt.ylabel("row", fontsize=14)
    plt.title('elnod amplitude, col %d, %s'%(col, modules[int(col/2)]))
    plt.savefig("%scol%d_elnod.png"%(out_path, col))

    plt.close()
    plt.pcolormesh(xi.T, yi.T, np.log10(abs(Acos[col,::-1,:])*1e12).T, vmin=1., vmax=5., cmap="coolwarm")
    cbar = plt.colorbar()
    cbar.set_label(r'$log10(a2, pA)$')
    #for ic in range(2):
    #    plt.axvline(x=33*ic, color='k', linestyle='--')
    #    plt.axhline(y=33*ic, color='k', linestyle='--')
    #plt.text(60, 60, modules[im], fontsize=14)
    plt.axvline(x=NETselectbias[col], color='k', linestyle='--', label='bias selected by NET script')
    plt.axvline(x=MAGselectbias[col], color='k', linestyle='-',label='bias adjusted by magnetic analysis')
    plt.legend()
    plt.xlabel("bias, ADU", fontsize=14)
    plt.ylabel("row", fontsize=14)
    plt.title('azscan amplitude (mag), mode1, col %d, %s'%(col, modules[int(col/2)]))
    plt.savefig("%scol%d_magaz_mode1.png"%(out_path, col))

    plt.close()
    plt.pcolormesh(xi.T, yi.T, np.log10(abs(Acos[col,::-1,:])/abs(Aelnod[col,::-1,:])).T, vmin=-2., vmax=0., cmap="coolwarm")
    cbar = plt.colorbar()
    cbar.set_label(r'$log10(a2/elnod)$')
    #for ic in range(2):
    #    plt.axvline(x=33*ic, color='k', linestyle='--')
    #    plt.axhline(y=33*ic, color='k', linestyle='--')
    #plt.text(60, 60, modules[im], fontsize=14)
    plt.axvline(x=NETselectbias[col], color='k', linestyle='--', label='bias selected by NET script')
    plt.axvline(x=MAGselectbias[col], color='k', linestyle='-',label='bias adjusted by magnetic analysis')
    plt.legend()
    plt.xlabel("bias, ADU", fontsize=14)
    plt.ylabel("row", fontsize=14)
    plt.title('azscan amplitude (mag) mode1/elnod, col %d, %s'%(col, modules[int(col/2)]))
    plt.savefig("%scol%d_magaz_mode1elnodratio.png"%(out_path, col))

# do some histogram
AcosNETbias = np.full(792,np.nan)
AelnNETbias = np.full(792,np.nan)
rNETbias = np.full(792,np.nan)
AcosMAGbias = np.full(792,np.nan)
AelnMAGbias = np.full(792,np.nan)
rMAGbias = np.full(792,np.nan)

for col in range(24):
    for row in range(33):
        ind=np.where(~np.isnan(Acos[col,::-1,row]))[0]
        if len(ind):
            AcosNETbias[33*col+row] = np.interp(NETselectbias[col], bl[::-1][ind], abs(Acos[col,::-1,row][ind]))
            AelnNETbias[33*col+row] = np.interp(NETselectbias[col], bl[::-1][ind], abs(Aelnod[col,::-1,row][ind]))
            rNETbias[33*col+row] = np.interp(NETselectbias[col], bl[::-1][ind], abs(Acos[col,::-1,row][ind])/abs(Aelnod[col,::-1,row][ind]))
            AcosMAGbias[33*col+row] = np.interp(MAGselectbias[col], bl[::-1][ind], abs(Acos[col,::-1,row][ind]))
            AelnMAGbias[33*col+row] = np.interp(MAGselectbias[col], bl[::-1][ind], abs(Aelnod[col,::-1,row][ind]))
            rMAGbias[33*col+row] = np.interp(MAGselectbias[col], bl[::-1][ind], abs(Acos[col,::-1,row][ind])/abs(Aelnod[col,::-1,row][ind]))

hist.plot_1Dhist(np.log10(AcosNETbias*1e12), out_path, 'AcosNETbias_hist',
    maintitle='azscan amplitude, bias selected by NET script',
    xlabel=r'$log_{10}(a2)$',
    xunit=r'pA',
    nbins=30,
    binrange=[1., 6.])

hist.plot_1Dhist(np.log10(AelnNETbias*1e12), out_path, 'AelnodNETbias_hist',
    maintitle='elnod amplitude, bias selected by NET script',
    xlabel=r'$log_{10}(elnod, pA)$',
    xunit=r'pA',
    nbins=30,
    binrange=[1., 6.])

hist.plot_1Dhist(np.log10(rNETbias), out_path, 'ratioNETbias_hist',
    maintitle='az/elnod, bias selected by NET script',
    xlabel=r'$log_{10}(az/elnod)$',
    xunit=r'pA',
    nbins=30,
    binrange=[-3., 1.])

hist.plot_1Dhist(np.log10(AcosMAGbias*1e12), out_path, 'AcosMAGbias_hist',
    maintitle='azscan amplitude, bias adjusted by magnetic analysis',
    xlabel=r'$log_{10}(a2, pA)$',
    xunit=r'pA',
    nbins=30,
    binrange=[1., 6.])

hist.plot_1Dhist(np.log10(AelnMAGbias*1e12), out_path, 'AelnodMAGbias_hist',
    maintitle='elnod amplitude, bias adjusted by magnetic analysis',
    xlabel=r'$log_{10}(elnod, pA)$',
    xunit=r'pA',
    nbins=30,
    binrange=[1., 6.])

hist.plot_1Dhist(np.log10(rMAGbias), out_path, 'ratioMAGbias_hist',
    maintitle='az/elnod, bias adjusted by magnetic analysis',
    xlabel=r'$log_{10}(az/elnod)$',
    xunit=r'pA',
    nbins=30,
    binrange=[-3., 1.])
