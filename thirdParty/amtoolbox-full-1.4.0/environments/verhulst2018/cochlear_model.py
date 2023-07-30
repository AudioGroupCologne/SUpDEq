# Part of the verhulst2018 model licensed under UGent license, see /licenses/ugent.txt
# Based on https://github.com/HearingTechnology/Verhulstetal2018Model (30 April 2020). 
# Alejandro Osses (2020) <ale.a.osses@gmail.com> Adaptations for the AMT
# Piotr Majdak (2021): Integrated in the AMT 1.0

    # -*- coding: utf-8 -*-
import numpy as np
import time
from scipy.integrate import ode
from scipy import signal
import ctypes
import os
import sys

DOUBLE = ctypes.c_double
INT=ctypes.c_int
PINT=ctypes.POINTER(ctypes.c_int)
PLONG=ctypes.POINTER(ctypes.c_long)
PDOUBLE = ctypes.POINTER(ctypes.c_double)

class tridiag_matrix(ctypes.Structure):
    _fields_=[("aa",ctypes.POINTER(ctypes.c_double)),
              ("bb",ctypes.POINTER(ctypes.c_double)),
              ("cc",ctypes.POINTER(ctypes.c_double))]
    
# Load C library
os.path.dirname(os.path.abspath(__file__))
tridiagLib = 'tridiag.dll' if 'win32' in sys.platform else 'tridiag.so'
try:
	libtrisolv=np.ctypeslib.load_library(tridiagLib,os.path.dirname(os.path.abspath(__file__)))
except IOError as err:
	raise RuntimeError("VERHULST2018 library (tridiag.dll or tridiag.so) not loaded. Run amt_mex.") from err

# Load tridiagonal solver function and defines input
libtrisolv.solve_tridiagonal.restype=None
libtrisolv.solve_tridiagonal.argtypes=[ctypes.POINTER(tridiag_matrix), #aa
                                       PDOUBLE, #vv
                                       PDOUBLE, #solution
                                       INT] #nrows

libtrisolv.delay_line.restype=None 
libtrisolv.delay_line.argtypes=[PDOUBLE,#in_matrix
                                PINT,    # delay1
                                PINT,    # delay2
                                PINT,    # delay3
                                PINT,    # delay4
                                PDOUBLE, # dev
                                PDOUBLE, # YZweig
                                INT,     # delay_buffer_length
                                INT]     # n

#######################################################
### Definition of the function
def TLsolver(t,y,model): #y''=dv/dt y'=v
    n=model.n+1
    frac=(t-model.lastT)/model.dt
    a=model.interplPoint1
    b=model.interplPoint2
    c=model.interplPoint3
    d=model.interplPoint4
    cminusb=c-b
    F0=b + frac * (cminusb - 0.1666667* (1.-frac) * ((d - a - 3.0 * cminusb) * frac + (d + 2.0*a - 3.0*b))) #fast cubic interpolation
    model.Vtmp=y[0:n]
    model.Ytmp=y[n:2*n]
    if model.version_year == 2012 or model.version_year == 2015:
        bDo = model.use_Zweig # either 0 (do not do) or 1 (do)
    elif model.version_year == 2018:
        bDo = model.non_linearity # either 0 (do not do) or 1 or 2 (do)
    else:
        bDo = model.non_linearity # same criterion as in version_year 2018

    # Updating poles if needed:
    SheraP_current = model.SheraP
    model.polecalculation(bDo)
    SheraP = model.SheraP # these are the new pole values
    # Updates non-linear parameters here only if the pole displacement is larger than 1 percent:
    if( np.max(abs(SheraP[1:n]-SheraP_current[1:n])/abs(SheraP_current[1:n]))>0.01):
        # model.SheraP=SheraP # new pole is updated into the array
        model.SheraParameters()
        model.ZweigImpedance()
        model.current_t=t
    else:
        model.SheraP=SheraP_current
    # End updating poles
         
    model.Dev[0:n]=model.Dev[0:n]+frac
    libtrisolv.delay_line( model.Ybuffer_pointer, model.Zrp_pointer, 
                           model.Zrp1_pointer, model.Zrp2_pointer, model.Zrp3_pointer, 
                           model.Dev_pointer, model.YZweig_pointer, 
                           ctypes.c_int(model.YbufferLgt), ctypes.c_int(model.n+1))

    model.Dev[0:n]=model.Dev[0:n]-frac
    model.calculate_g()
    model.calculate_right(F0)

    libtrisolv.solve_tridiagonal(ctypes.byref(model.tridata),model.r_pointer,model.Qpointer,ctypes.c_int(n)) # compute q
    zero_val=(model.RK4_0*model.Qsol[0]+model.RK4G_0*(model.g[0]+model.p0x*F0)) # ME equation

    Vderivative=(model.Qsol-model.g)  
    Vderivative[0]=zero_val
    solution=np.concatenate([Vderivative,model.Vtmp]) #output the velocity input as displacement derivative output. 
    return solution     
       
###########################################################
class cochlea_model ():
   
    #initialise constants
    def __init__(self):
        self.ttridiag=0
        self.calling_function=0
        self.cochleaLength=.035
        self.bmMass=0.5
        self.bmImpedanceFactor=1. 
        self.scalaWidth=0.001
        self.scalaHeight=0.001
        self.helicotremaWidth=0.001
        self.rho=float(1e3)
        # self.Normal_Q=20. # Not used
        self.Greenwood_A=20682.
        self.Greenwood_alpha=61.765
        self.Greenwood_B=140.6
        self.stapesArea=float(3e-6)
        self.EardrumArea=float(60e-6)
        self.MiddleEarResonanceFrequency=float(2e3)
        self.MiddleEarQualityFactor=0.4
        # self.SpecificAcousticImpedanceOfAir=415. # Not used
        # self.middleEarTransformer=30. # not used
        # self.damping_coupler=float(140e5) # Not used
        # self.mass_coupler=float(43.2e2) # Not used
        # self.stiffness_coupler=1./float(2.28e-11) # Not used
        # self.p0=float(2e-5) # Not used
        # self.ZweigQ=1/0.0606 # Not used
        # self.ZweigFactor=1.7435 # Not used
        # self.ZweigQAtBoundaries=20. # Not used
        # self.ZweigBeta=10000. # Not used
        # self.ZweigGamma=6200. # Not used
        self.ZweigN=1.5
        self.RMSref=0.6124
        self.Rme=float(0.3045192500000000e12) #TODO setRme function: np.sqrt(self.ZweigOmega_co**2*self.ZweigMpo*self.ZweigMso)
        self._is_init=0 #variable to check if the model is intialize before calling the solver
        self.interplPoint1=0
        self.interplPoint2=0
        self.interplPoint3=0
        self.interplPoint4=0
    
    #######################################################
    # Function to intitialise all the parameters        
    def init_model(self, stim, fs, sections, probe_freq, sheraPo=0.061, compression_slope=0.4, Zweig_irregularities=1, non_linearity_type="vel", KneeVar=1., low_freq_irregularities=1, subject=1, IrrPct=0.05, version_year=2012):

        self.version_year = version_year
        if self.version_year == 2012 or self.version_year == 2015:
            self.SheraMuMax=4.3
            self.SheraPo=sheraPo #can be vector or single value

        elif version_year == 2018:
            self.SheraMuMax = 3.
            self.SheraPo = np.zeros(sections+1)
            self.SheraPo = self.SheraPo+sheraPo # sheraPo can be either a vector or a single value

        self.low_freq_irregularities=low_freq_irregularities
        self.KneeVar=(KneeVar)
        self.IrrPct=IrrPct
        if(Zweig_irregularities==0):
            self.use_Zweig=0
        else:
            self.use_Zweig=1 # i.e., use "random BM impedance irregularities" (Verhulst et al. 2012)
        if(non_linearity_type=="disp"):
            self.non_linearity=1
        elif(non_linearity_type=="vel"):
            self.non_linearity=2
        else:
            self.non_linearity=0 #linear model
        self.n=sections
        self.fs=fs
        self.dt=1./self.fs
        self.probe_freq=probe_freq

        self.initCochlea()
        self.initMiddleEar()
        self.SetDampingAndStiffness()
        self.initZweig()
        self.initGaussianElimination()
        self.compression_slope_param(compression_slope)
        self.is_init=1
        self.lastT=0
        self.seed=subject # Change here the seed
        np.random.RandomState(self.seed)
        np.random.seed(self.seed)
        self.Rth=2*(np.random.random(self.n+1)-0.5) # random numbers between -1 and 1
        self.Rth_norm=10**(self.Rth/20./self.KneeVar)
        lf_limit=self.ctr
        if(self.use_Zweig==0):
            lf_limit=0
            print('   No irregularities for at least one of the signals')

        if self.version_year == 2012 or self.version_year == 2015:
            bDo = self.use_Zweig
        elif self.version_year == 2018:
            bDo = 1 # i.e., the next lines are always compiled (cochlear non-linearity cannot be bypassed)

        if bDo:
            factorA=self.factorA # from compression_slope_param
            n=self.n+1
            Rth=self.Rth
            Rth_norm=self.Rth_norm
            # Normalised RTH, it saves some computation
            self.RthY1=self.Yknee1*Rth_norm 
            self.RthY2=self.Yknee2*Rth_norm
            self.RthV1=self.Vknee1*Rth_norm
            self.RthV2=self.Vknee2*Rth_norm
            
            Rndm=self.IrrPct*Rth/2.
            self.PoleS=(1+Rndm)*self.SheraPo

            self.RthY1[lf_limit:n]=self.Yknee1
            self.RthY2[lf_limit:n]=self.Yknee2
            self.RthV1[lf_limit:n]=self.Vknee1
            self.RthV2[lf_limit:n]=self.Vknee2
            if version_year == 2012:
                self.PoleS[lf_limit:n]=self.SheraPo
            elif version_year == 2018 or version_year == 2015:
                self.PoleS[lf_limit:n] = self.SheraPo[lf_limit:n]
    
            ratioV2_V1 = (self.RthV2/self.RthV1)
            Theta=0.5*np.arctan(((self.PoleE-self.PoleS)*factorA)/(ratioV2_V1-1.)) # Eq. 5, verhulst2018

            F_p = (self.PoleS * factorA) / ratioV2_V1 # F_p in Eq. 8, verhulst2018
            self.Sb=F_p*np.sin(Theta) # b_p in Eq. 8, verhulst2018
            self.Sa=F_p*np.cos(Theta) # a_p in Eq. 8, verhulst2018
            self.const_nl1=np.cos(Theta)/np.cos(2*Theta)
            self.Theta = Theta
        
        # PURIAM1 FILTER removed, implemented in MATLAB in 
        #     middleearfilter.m (flags 'verhulst2012', or 'verhulst2018')
        self.stim=stim
      
    # Adapted from Fortran script: intializeCochlea.f90
    def initCochlea(self):
        self.bm_length = self.cochleaLength - self.helicotremaWidth
        self.bm_width = self.scalaWidth
        self.bm_mass = self.bmMass*self.bmImpedanceFactor
        self.ZweigMso = 2.* self.rho / ( self.bm_width * self.scalaHeight )
        self.ZweigL = 1. / (2.3030 * self.Greenwood_alpha)

        if self.version_year == 2012 or self.version_year == 2015:
            self.ZweigOmega_co = 2.0 * np.pi * self.Greenwood_A # w_c0 in Table I of Verhulst2012
        elif self.version_year == 2018:
            # Small deviation in 2018-code (this change is not reported):
            self.ZweigOmega_co = 2.0 * np.pi * self.Greenwood_A-self.Greenwood_B

        self.ZweigMpo = (self.ZweigMso * (self.ZweigL**2)) / ( (4 * self.ZweigN)**2) # M_p0, acoustical mass at base
        self.Ko = self.ZweigMpo * (self.ZweigOmega_co**2) # K_0, stiffness constant
        self.x=np.array(np.linspace(0,self.bm_length,self.n+1),order='C') 
        self.dx=self.bm_length/(1.*self.n)
        
        ### Other variables to be initialised here for practical purposes:        
        self.g=np.zeros_like(self.x)
        self.Vtmp=np.zeros_like(self.x)
        self.Ytmp=np.zeros_like(self.x)
        if self.version_year == 2018:
            self.Atmp= np.zeros_like(self.x)
        self.right=np.zeros_like(self.x)
        self.r_pointer=self.right.ctypes.data_as(PDOUBLE)
        self.zerosdummy=np.zeros_like(self.x)
        self.gamma=np.zeros_like(self.x)
        self.Qsol=np.zeros_like(self.x)
        self.Qpointer=self.Qsol.ctypes.data_as(PDOUBLE)

    def initMiddleEar(self):
        self.q0_factor=self.ZweigMpo*self.bm_width
        self.p0x=self.ZweigMso*self.dx/(1.*self.ZweigMpo*self.bm_width)
        self.d_m_factor= -self.p0x*self.stapesArea*self.Rme
        self.RK4_0=-(self.bm_width*self.ZweigMpo)/(self.stapesArea)
        self.RK4G_0=(self.ZweigMpo*self.bm_width)/(self.ZweigMso*self.stapesArea*self.dx)

    def SetDampingAndStiffness(self):
        self.f_resonance=self.Greenwood_A*10**(-self.Greenwood_alpha*self.x)-self.Greenwood_B # Greenwood1990, Eq. (1) but from base-to-apex
        self.ctr=np.argmin(np.abs(self.f_resonance-100.)) # finds idx of bin with freq=100 Hz
        if(self.low_freq_irregularities):
            self.ctr=self.n+1
        self.onek=np.argmin(np.abs(self.f_resonance-1000.))  # finds idx of bin with freq=1000 Hz
        self.omega=2.*np.pi*self.f_resonance
        if self.version_year == 2018:
            self.omega[0]=self.ZweigOmega_co
        self.omega2=self.omega**2
        self.Sherad_factor=np.array(self.omega)
        self.SheraP=np.zeros_like(self.x)
        self.SheraD=np.zeros_like(self.x)
        self.SheraRho=np.zeros_like(self.x)
        self.SheraMu=np.zeros_like(self.x)
        self.SheraP=self.SheraPo+self.SheraP
        self.c=120.8998691636393
        
        ### PROBE POINTS:
        if(self.probe_freq=='all'):
            self.probe_points=np.zeros(len(self.f_resonance)-1)
            for i in range(len(self.f_resonance)-1):
                self.probe_points[i]=i+1
            self.probe_points=(self.probe_points)
            if self.version_year == 2018:
                self.cf=(self.f_resonance[1:len(self.f_resonance)])  # slight difference in Verhulst2018 model implementation
            else:
                self.cf=(self.f_resonance[0:len(self.f_resonance)])
        elif(self.probe_freq=='half'):
            self.probe_points=np.array(range(1,self.n+1,2))
            self.cf=(self.f_resonance[range(1,self.n+1,2)])
        elif(self.probe_freq=='abr'):
            self.probe_points=np.array(range(110,911,2))
            self.cf=(self.f_resonance[range(110,911,2)])
        else:
            self.probe_points=np.zeros_like(self.probe_freq)
            for i in range(len(self.probe_freq)):
                self.probe_points[i]=np.argmin(abs(self.f_resonance-self.probe_freq[i]))
            self.cf=self.f_resonance[self.probe_points]
    def initZweig(self):
        n=self.n+1
        self.exact_delay=self.SheraMuMax/(self.f_resonance*self.dt)
        self.delay=np.floor(self.exact_delay)+1
        self.YbufferLgt=int(np.amax(self.delay))
        self.Ybuffer=np.zeros([n,self.YbufferLgt]) #Ybuffer implemented here as a dense matrix (python for cycles are slow...)
        self.Ybuffer=np.array(self.Ybuffer,order='C',ndmin=2,dtype=float)
        self.Ybuffer_pointer=self.Ybuffer.ctypes.data_as(PDOUBLE)
        self.ZweigSample1=np.zeros_like(self.exact_delay)
        self.Zwp=int(0)
        self.ZweigSample1[0]=1.
        self.ZweigSample2=self.ZweigSample1+1
        #init buffers etc...
        self.Dev=np.zeros_like(self.x)
        self.Dev_pointer=self.Dev.ctypes.data_as(PDOUBLE)
        self.YZweig=np.zeros_like(self.x)
        self.YZweig_pointer=self.YZweig.ctypes.data_as(PDOUBLE)
        self.Zrp=np.array(np.zeros(n),dtype=np.int32,order='C')
        self.Zrp_pointer=self.Zrp.ctypes.data_as(PINT)
        self.Zrp1=np.array(np.zeros(n),dtype=np.int32,order='C')
        self.Zrp1_pointer=self.Zrp1.ctypes.data_as(PINT)
        self.Zrp2=np.array(np.zeros(n),dtype=np.int32,order='C')
        self.Zrp2_pointer=self.Zrp2.ctypes.data_as(PINT)
        self.Zrp3=np.array(np.zeros(n),dtype=np.int32,order='C')
        self.Zrp3_pointer=self.Zrp3.ctypes.data_as(PINT)

    def initGaussianElimination(self): #set tridiagonal matrix values for trasmission line
        n=self.n+1
        self.ZweigMs=(self.ZweigMso*self.ZweigOmega_co)/self.omega
        self.ZweigMp=self.Ko/(self.ZweigOmega_co*self.omega) 
        self.ZASQ=np.zeros_like(self.x)
        self.ZASC=np.zeros_like(self.x)
        self.ZAL=np.zeros_like(self.x)
        self.ZAH=np.zeros_like(self.x)
        #init values of transimission line
        self.ZASQ[0]=1.
        self.ZASC[0]=1+self.ZweigMso*self.dx
        if self.version_year == 2012 or self.version_year == 2015:
            self.ZAH[0]=-1.
            self.ZAL[1:n]=-self.ZweigMs[1:n]
            self.ZAH[1:n-1]=-self.ZweigMs[0:n-2]

        elif self.version_year == 2018:
            self.ZAH[0] = -1.*self.ZweigOmega_co/self.omega[1]
            self.ZAL[1:n] = -self.ZweigMs[1:n]*self.omega[1:n]/self.omega[0:n-1]
            self.ZAH[1:n-1] =-self.ZweigMs[0:n-2]*self.omega[1:n-1]/self.omega[2:n]

        self.ZASQ[1:n]=self.omega[1:n]*self.ZweigMs[1:n]*self.ZweigMs[0:n-1]*(self.dx**2)/(self.ZweigOmega_co*self.ZweigMpo)
        self.ZASC[1:n]=self.ZASQ[1:n]+self.ZweigMs[1:n]+self.ZweigMs[0:n-1]
        self.tridata=tridiag_matrix()
        self.tridata.aa=self.ZAL.ctypes.data_as(PDOUBLE)
        self.tridata.bb=self.ZASC.ctypes.data_as(PDOUBLE)
        self.tridata.cc=self.ZAH.ctypes.data_as(PDOUBLE)
        
    def calculate_g(self): # Same as in the fortran code
        n=self.n+1
        self.g[0]=self.d_m_factor*self.Vtmp[0]
        dtot=self.Sherad_factor*self.SheraD
        stot=(self.omega2)*(self.Ytmp+(self.SheraRho*self.YZweig))
        self.g[1:n]=(dtot[1:n]*self.Vtmp[1:n])+stot[1:n]
        self.passive=(dtot[0:n]*self.Vtmp+self.SheraRho*self.YZweig*self.omega2)
        
    def calculate_right(self,F0): # Same as in the fortran code
        n=self.n+1
        self.right[0] = self.g[0] +self.p0x*F0
        self.right[1:n]=self.ZASQ[1:n]*self.g[1:n]
        
    def SheraParameters(self): # Same as in the fortran code
        a=(self.SheraP+np.sqrt((self.SheraP**2.)+self.c*(1.0-self.SheraP**2)))/self.c
        self.SheraD=2.0*(self.SheraP-a)
        self.SheraMu=1./(a)
        self.SheraRho=2.* a * np.sqrt(1.-(self.SheraD/2.)**2.) * np.exp(-self.SheraP/a)
        
    
    def ZweigImpedance(self):
        n=self.n+1
        MudelayExact=self.SheraMu/(self.omega*self.dt)
        Mudelay=np.floor(MudelayExact)+1.
        self.Dev[:]=Mudelay-MudelayExact
        self.Zrp1[0:n]=((self.Zwp+self.YbufferLgt)-Mudelay[0:n])%self.YbufferLgt
        const=self.YbufferLgt-1
        self.Zrp[0:n]=(self.Zrp1[0:n]+const)%self.YbufferLgt
        self.Zrp2[0:n]=(self.Zrp1[0:n]+1)%self.YbufferLgt
        self.Zrp3[0:n]=(self.Zrp2[0:n]+1)%self.YbufferLgt

    #######################################################
    def compression_slope_param(self,slope):

        factorA = 100.
        self.factorA = factorA
        self.Yknee1=float(6.9183e-10)
        if self.version_year == 2012 or self.version_year == 2015:
            self.Vknee1=float(4.3652e-6)
            self.PoleE=np.zeros_like(self.x)
            if(slope==0.2):
                self.Yknee2=float(3.228e-9) 
                self.Vknee2=float(2.037e-5) 
                Ax=(0.7-0.06)/(80.59-30.)
                Bx=self.SheraPo-Ax*30
                self.PoleE=self.PoleE+Ax*80.59+Bx
            elif(slope==0.3):
                self.Yknee2=float(7.015e-9) 
                self.Vknee2=float(4.426e-5) 
                Ax=(0.7-0.06)/(87.77-30.)
                Bx=self.SheraPo-Ax*30.
                self.PoleE=self.PoleE+Ax*87.77+Bx
            elif(slope==0.5):
                self.Yknee2=float(1.766e-8) 
                self.Vknee2=float(1.114e-4) 
                Ax=(0.7-0.06)/(97.82-30.)
                Bx=self.SheraPo-Ax*30.
                self.PoleE=self.PoleE+Ax*97.82+Bx
            else:
                self.Yknee2=float(1.5488e-8) 
                self.Vknee2=float(9.7836e-5) 
                Ax=(0.7-0.06)/(97.4-30.)
                Bx=self.SheraPo-Ax*30.
                self.PoleE=self.PoleE+Ax*97.4+Bx
        elif self.version_year == 2018: 
            self.Yknee2 = float(1.5488e-8)

            # 'New compression' by Alessandro Altoe (comments by Alejandro Osses):
            self.PoleE = np.zeros_like(self.x)+0.31 # Saturating pole
            gain_m10 = 1./3. # gain at approximately -10 dB, it is exactly at -9.54 dB 
                             # (the compression will actually start at an instant level of 30.46 dB, not at 30 dB as expected)
            v1=gain_m10*0.6807e-08/np.sqrt(2) # velocity at "-10 dB" (-9.54 dB) for starting Pole
            v2=gain_m10*26.490e-11/np.sqrt(2) # velocity at "-10 dB" (-9.54 dB) for saturating pole

            K1dB=40.; # Knee point of the first linear regime in dB (you can select it from here now). 
                      # AO: my hypothesis is that if you are at -10 you need 40 dB to reach 30 dB
            K1L=10**(K1dB/20) # Knee point in linear scale
            Vknee1=K1L*v1
            self.Vknee1=Vknee1

            vst1dB=20*np.log10(Vknee1) # Velocity with the two poles when the compression starts
            vst2dB=20*np.log10(v2)+K1dB
            K2dB=(vst1dB-vst2dB)/(1-slope)+K1dB # Intersection in dB re Knee1. '1-slope' because the level in the linear scale is our incognita
            Vknee2 = v2*10**(K2dB/20)
            self.Vknee2=v2*10**(K2dB/20)

    #######################################################    
    def polecalculation(self,bDo=-1):

        factorA=self.factorA
        # lf_limit=self.ctr
        n=self.n+1
        if bDo == -1: # i.e., if default input is given
            if self.version_year == 2012 or self.version_year == 2015:
                bDo = self.use_Zweig
            elif self.version_year == 2018:
                bDo = 1
        
        if bDo:
            if(self.non_linearity==1): #To check
                # Non-linearity DISPLACEMENT ('disp') costs about three times more than in fortran (Not implemented now)
                Yknee1CST=self.RthY1*self.omega[self.onek]
                Yknee2CST=self.RthY2*self.omega[self.onek]
                Yknee1F=Yknee1CST/self.omega
                Yknee2F=Yknee2CST/self.omega
                Yvect=np.abs(self.Ytmp/Yknee1F)
                ratioY2_Y1=Yknee2F/Yknee1F
                Theta=0.5*np.arctan(((self.PoleE-self.PoleS)/(ratioY2_Y1-1.)))

                cos_Theta0=2*np.cos(Theta)**2-1
                F_p=self.PoleS*factorA*ratioY2_Y1
                Sb=F_p*np.sin(Theta)
                Sa=F_p*np.cos(Theta)
                Sxp=(Yvect-1.)*cos_Theta/cos_Theta0
                Syp=Sb*np.sqrt(1+(Sxp/Sa)**2)
                Sy  = Sxp*np.sin(Theta)+Syp*np.cos(Theta)
                self.SheraP=self.PoleS+Sy/factorA
		    
            elif(self.non_linearity==2):  #non-linearity VEL
                Vvect=np.abs(self.Vtmp)/self.RthV1
                Sxp=(Vvect-1.)*self.const_nl1
                Syp=self.Sb*np.sqrt(1+(Sxp/self.Sa)**2)
                Sy=Sxp*np.sin(self.Theta)+Syp*np.cos(self.Theta)
                self.SheraP=self.PoleS+Sy/factorA
            else:
                print('linear')
                self.SheraP = self.PoleS

        self.SheraP=np.fmin(self.SheraP,self.PoleE)

    #######################################################
    def solve(self):
        n=self.n+1
        tstart=time.time()
        if not(self.is_init):
           print("Error: the model needs to be initialised")
        length=np.size(self.stim) - 2
        time_length=length*self.dt

        # Each probe point signal in a row:
        self.Vsolution=np.zeros([length+2,len(self.probe_points)]) 
        self.Ysolution=np.zeros([length+2,len(self.probe_points)])
        self.oto_emission=np.zeros(length+2)

        # Check next two variables (Asolution,organ_of_corti_acceleration): Maybe unused in '2018' and '2012', respectively... 
        if self.version_year == 2018: 
            self.Asolution= np.zeros([length + 2, len(self.probe_points)])
        elif self.version_year == 2012 or self.version_year == 2015:
            self.organ_of_corti_acceleration=np.zeros([length+2,len(self.probe_points)]) 

        self.time_axis=np.linspace(0,time_length,length)

        # Constructor to ODE:
        r=ode(TLsolver).set_integrator('dopri5',rtol=1e-2,atol=1e-13)
        r.set_f_params(self)
        r.set_initial_value(np.concatenate([np.zeros_like(self.x),np.zeros_like(self.x)]))
        r.t=0
        j=0
        self.last_t=0.0
        self.current_t=r.t
        self.polecalculation()
        self.SheraParameters()
        self.ZweigImpedance()
        self.V1=np.zeros_like(self.x) # in older codes 'V1 was named Vtmp'

        while(j<length):
            # Assign the stimulus points and interpolation parameters:
            if(j>0):
                self.interplPoint1=self.stim[j-1]
            self.interplPoint2=self.stim[j]
            self.interplPoint3=self.stim[j+1] 
            self.interplPoint4=self.stim[j+2]
            r.integrate(r.t+self.dt)
            self.lastT=r.t
            self.V1=r.y[0:n]
            self.V1[0]=self.Vtmp[0]
            self.Y1=r.y[n:2*n] #Non linearities HERE
            if self.version_year == 2018:
                self.Atmp=self.Qsol-self.g

            self.Zwp=(self.Zwp+1)%self.YbufferLgt #update Zweig Buffer
            self.Ybuffer[:,self.Zwp]=self.Y1
            self.ZweigImpedance()
            self.current_t=r.t

            if(self.probe_freq=='all'):
                self.Vsolution[j,:]=self.V1[1:n] #storing all probe points (Debug)
                self.Ysolution[j,:]=self.Y1[1:n]
                if self.version_year == 2012 or self.version_year == 2015:
                    self.organ_of_corti_acceleration[j,:]=self.passive[1:n]
            elif(self.probe_freq=='half'):
                self.Vsolution[j,:]=self.V1[range(1,n,2)]
                self.Ysolution[j,:]=self.Y1[range(1,n,2)]
                if self.version_year == 2012 or self.version_year == 2015:
                    self.organ_of_corti_acceleration[j,:]=self.passive[range(1,n,2)]
            elif(self.probe_freq=='abr'):
                self.Vsolution[j,:]=self.V1[range(110,911,2)]
                self.Ysolution[j,:]=self.Y1[range(110,911,2)]
                if self.version_year == 2012 or self.version_year == 2015:
                    self.organ_of_corti_acceleration[j,:]=self.passive[range(110,911,2)]
            else:
                self.Vsolution[j,:] = self.V1[self.probe_points]  # storing the decided probe points
                self.Ysolution[j,:] = self.Y1[self.probe_points]
                if self.version_year == 2012 or self.version_year == 2015:
                    self.organ_of_corti_acceleration[j,:]=self.passive[self.probe_points]
                # ### Old code:
                # k = 0
                # for i in self.probe_points:
                #     self.Vsolution[j,k]=self.V1[i] 
                #     self.Ysolution[j,k]=self.Y1[i]
                #     self.organ_of_corti_acceleration[j,k]=self.passive[i]
                #     k=k+1
            self.oto_emission[j]=self.Qsol[0]
            j=j+1
            # End of 'while'

        #### filter out the otoacoustic emission (Reverse filter) ####
        fs=self.fs
        if self.version_year == 2012 or self.version_year == 2015:
            b,a=signal.butter(1,[600./(fs/2.),3000./(fs/2.)],'bandpass')
        elif self.version_year == 2018:
            b,a=signal.butter(1,[600./(fs/2.),4000./(fs/2.)],'bandpass')
        self.oto_emission=signal.lfilter(b*self.q0_factor,a,self.oto_emission)
        elapsed = time.time()-tstart

#END
