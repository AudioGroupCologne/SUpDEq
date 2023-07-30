#The AMT Team

# This file is licensed unter the GNU General Public License (GPL) either 
# version 3 of the license, or any later version as published by the Free Software 
# Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
# at <https://www.gnu.org/licenses/gpl-3.0.html>. 
# You can redistribute this file and/or modify it under the terms of the GPLv3. 
# This file is distributed without any warranty; without even the implied warranty 
# of merchantability or fitness for a particular purpose.

#!/usr/bin/env python3

import numpy as np
import scipy as sp
from scipy import signal
import scipy.io as sio
import cochlear_model
import os
#def erb2hz(x):
#    return 229*(10.**(x/21.3)-1)
#def next_pow2(x):
#    return 2**(np.ceil(np.log2(x)))
import multiprocessing as mp
import ctypes as c
import time

#import sys
#print (sys.version)
Oversampling=1
sectionsNo=1000
p0=float(2e-5)
par=sio.loadmat('out/input.mat')

probes=np.array(par['probes']) 
probe_points=probes
Fs=par['Fs']
Fs=Fs[0][0]
stim=par['stim']
spl=par['spl'] 
spl=np.array(spl[0])
channels=par['channels']
channels=int(channels[0][0])
subjectNo=int(par['subject'])
lgt=len(stim[0])
norm_factor=p0*10.**(spl/20.)
sheraPo=0.06
irr_on=np.array(par['irregularities'])
d=len(stim[0].transpose())
print("running cochlear simulation")

for i in range(channels):
    stim[i]=stim[i]*norm_factor[i]
    
sig=stim

cochlear_list=[ [cochlear_model.cochlea_model(),sig[i],irr_on[0][i],i] for i in range(channels)]

def solve_one_cochlea(model): #definition here, to have all the parameter implicit
    i=model[3]
    coch=model[0]
    coch.init_model(model[1],Oversampling*Fs,sectionsNo,probe_points,Zweig_irregularities=model[2],sheraPo=sheraPo,subject=subjectNo) #model needs to be init here because if not pool.map crash
    coch.solve()
    f=open("out/v"+str(i+1)+".np",'wb')
    np.array(coch.Vsolution,dtype='=d').tofile(f)
    f.close()
    f=open("out/y"+str(i+1)+".np",'wb')
    np.array(coch.Ysolution,dtype='=d').tofile(f)
    f.close()
    f=open("out/E"+str(i+1)+".np",'wb')
    np.array(coch.oto_emission,dtype='=d').tofile(f)
    f.close()
    f=open("out/F"+str(i+1)+".np",'wb')
    np.array(coch.cf,dtype='=d').tofile(f)
    f.close()
    return

if __name__ == "__main__":
    p=mp.Pool(int(mp.cpu_count()),maxtasksperchild=1)
    p.map(solve_one_cochlea,cochlear_list)
    p.close()
    p.join()

    
