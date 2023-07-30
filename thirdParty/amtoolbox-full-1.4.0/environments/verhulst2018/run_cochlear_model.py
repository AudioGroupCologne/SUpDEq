# Part of the verhulst2018 model licensed under UGent license, see /licenses/ugent.txt
# Based on https://github.com/HearingTechnology/Verhulstetal2018Model (30 April 2020). 
# Alejandro Osses (2020) <ale.a.osses@gmail.com> Adaptations for the AMT
# Piotr Majdak (2021): Integrated in the AMT 1.0

import numpy as np
import scipy as sp
from scipy import signal
import scipy.io as sio
import cochlear_model
import os

import multiprocessing as mp
import ctypes as c
import time
import warnings

# Original comment from UGent code: this relates to python 3.6 on ubuntu
#there is one future warning related to "scipy.signal.decimate" in this file
#there is one runtime warning related to firwin "scipy.signal.decimate" in ic_cn2017.py (not important)
#so we suppress these warnings here
warnings.filterwarnings("ignore")

par=sio.loadmat('out/input.mat')

version_year = par['version_year']
version_year = version_year[0][0]

probes=np.array(par['probes']) 
probe_points=probes

storeflag_in=np.array(par['storeflag'],dtype=str)
storeflag=storeflag_in[0]

fs=par['fs']
fs=fs[0][0]  # Sampling frequency [Hz]

stim=par['stim']
channels=par['channels']
channels=int(channels[0][0])

subjectNo=int(par['subject'])
lgt=len(stim[0])

IrrPct=par['IrrPct']
IrrPct=IrrPct[0][0]
non_lin_type=np.array(par['non_linear_type'])

### Fixed values (that won't be changed from MATLAB)
sectionsNo=1000
compression_slope=0.4 # dB/dB

sheraPo=par['sheraPo']
if version_year == 2012:
    sheraPo=sheraPo[0][0] # default: sheraPo = 0.06

elif version_year == 2018 or version_year == 2015:
    if(max(np.shape(sheraPo))==1):
        sheraPo=sheraPo[0][0]
    else:
        sheraPo=sheraPo[:,0]
###
irr_on=np.array(par['irregularities'])
d=len(stim[0].transpose())

#if version_year == 2012:
#    print("Running cochlear simulation 2012: Verhulst, Dau, Shera (2012)")
#elif version_year == 2015:
#    print("Running human auditory model 2015: Verhulst, Bharadwaj, Mehraei, Shera, Shinn-Cunningham (2015)")
#elif version_year == 2018:
#    print("Running human auditory model 2018: Verhulst, Altoe, Vasilkov (2018),")
#    print("    in its version 1.2 (Osses and Verhulst, 2019)")
#else:
#    print("Running cochlear simulation")

insig=stim
#                                          model: [1] ,     [2]      [3]
cochlear_list=[ [cochlear_model.cochlea_model(),insig[i],irr_on[0][i],i] for i in range(channels)]

def solve_one_cochlea(model): #definition here, to have all the parameter implicit
    coch  = model[0]
    insig = model[1]
    irr_on= model[2] # Zweig_irregularities, on/off
    ii    = model[3]

    ### Solve the TL cochlear model:
    #    Important - defaults will be used for 'non_linearity_type','KneeVar','low_freq_irregularities',
    #                'IrrPct' because they are not indicated as inputs for coch.init_model (next line)
    coch.init_model(insig,fs, sectionsNo, probe_points, sheraPo, compression_slope, irr_on, subject=subjectNo, IrrPct=IrrPct, non_linearity_type=non_lin_type,version_year=version_year) 
    coch.solve()
    f=open("out/cf"+str(ii+1)+".np",'wb')
    np.array(coch.cf,dtype='=d').tofile(f)
    f.close()

    f=open("out/v"+str(ii+1)+".np",'wb')
    np.array(coch.Vsolution,dtype='=d').tofile(f)
    f.close()

    if 'y' in storeflag:
        f=open("out/y"+str(ii+1)+".np",'wb')
        np.array(coch.Ysolution,dtype='=d').tofile(f)
        f.close()

    if 'e' in storeflag:
        f=open("out/e"+str(ii+1)+".np",'wb')
        np.array(coch.oto_emission,dtype='=d').tofile(f)
        f.close()

    return

if __name__ == "__main__":
    p=mp.Pool(int(mp.cpu_count()),maxtasksperchild=1)
    p.map(solve_one_cochlea,cochlear_list)
    p.close()
    p.join()
    # print("Cochlear simulation: done")
