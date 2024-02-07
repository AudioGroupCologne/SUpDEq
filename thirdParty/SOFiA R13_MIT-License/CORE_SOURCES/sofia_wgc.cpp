/**********************************************************************
%  
% /// ASAR/MARA Research Group
%  
% Technology Arts Sciences TH Köln
% Technical University of Berlin
% Deutsche Telekom Laboratories
% University of Rostock
% WDR Westdeutscher Rundfunk
% IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis
% 
% W/G/C Wave Generator Core R13-0306
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com                        
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
%  
% [Pnm, kr] = sofia_wgc(N, r, ac, FS, NFFT, AZ, EL, 
%                             t, c, wavetype, ds, lSegLim, uSegLim, SeqN)
% ------------------------------------------------------------------------
% Pnm      Spatial Fourier Coefficients 
%          Columns : nm coeff
%          Rows    : FFT bins 
% kr       kr-Vector 
%          Can also be a matrix [krm; krs] for rigid sphere configurations:
%          [1,:] => krm referring to the Microphone Radius
%          [2,:] => krs referring to the Sphere Radius (Scatterer)
%  
% ------------------------------------------------------------------------   
% N        Maximum transform order
% r        Microphone Radius 
%          Can also be a vector for rigid/dual sphere configurations:
%          [1,1] => rm  Microphone Radius 
%          [2,1] => rs  Sphere Radius or Microphone2 Radius
%          ! If only one radius (rm) is given using a Rigid/Dual Sphere  
%            Configuration: rs = rm and only one kr-vector is returned!
% ac       Array Configuration 
%          0  Open Sphere with pressure Transducers (NO plc!)
%          1  Open Sphere with cardioid Transducers
%          2  Rigid Sphere with pressure Transducers
%          3  Rigid Sphere with cardioid Transducers (Thx to Nils Peters!)
%          4  Dual Open Sphere with pressure Transducers (Thx to Nils Peters!)
% FS       Sampling Frequency
% NFFT     FFT Order (Number of bins) should be 2^x, x=1,2,3,... 
% AZ       Azimuth   angle in [DEG] 0-2pi            
% EL       Elevation angle in [DEG] 0-pi                   
% t        Time Delay in s. The delay has: (t*FS) Samples
% c        Speed of sound in [m/s] (Default: 343m/s)
% wavetype Type of the Wave. 0: Plane Wave (default) 1: Spherical Wave 
% ds       Distance of the source in [m] (For wavetype = 1 only)
%          Warning: If NFFT is smaller than the time the wavefront 
%          needs to travel from the source to the array, the impulse 
%          response will by cyclically shifted (cyclic convolution). 
% ---
% lSegLim  (Lower Segment Limit) Used by the S/W/G wrapper
% uSegLim  (Upper Segment Limit) Used by the S/W/G wrapper
% SegN     (Sement Order)        Used by the S/W/G wrapper
%
@ end of header
%
% CONTACT AND LICENSE INFORMATION:
% 
% /// ASAR/MARA Research Group 
%  
%     [1] Technology Arts Sciences TH Köln
%     [2] Technical University of Berlin 
%     [3] Deutsche Telekom Laboratories 
%     [4] University of Rostock
%     [5] WDR Westdeutscher Rundfunk 
%     [6] IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis toolbox
% 
% Copyright 2011-2017 Benjamin Bernschütz et al.(§)  
% 
% Contact ------------------------------------
% Technology Arts Sciences TH Köln 
% Institute of Communications Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
% 
% phone       +49 221 8275 -2496 
% cell phone  +49 171 4176069 
% mail        rockzentrale 'at' me.com 
% --------------------------------------------
% 
% This file is part of the SOFiA sound field analysis toolbox
%
% Licence Type: MIT License
%
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
%
% (§) Christoph Pörschmann [1]     christoph.poerschmann 'at' th-koeln.de
%     Sascha Spors         [2,3,4] sascha.spors 'at' uni-rostock.de  
%     Stefan Weinzierl     [2]     stefan.weinzierl 'at' tu-berlin.de
% 
**********************************************************************/

#include <mex.h>
#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <sofia_radial.h>   

#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif
       
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
   using namespace std;
   using namespace boost::math;
           
   unsigned int           N, ac, FS, F_NFFT, NFFT, wavetype; 
   unsigned int           NMLocatorSize, lowerSegLim, upperSegLim, SegN; 
   int                    n, m, nor;
   long long unsigned int ctr, ctrb, ctrc, f;
   double                 az, el, rm, rs, t, c, ds, maxw, maxk, iteratew, iteratek,iteratem, iterates, iterated;
   double                 *r, *w, *k, *krm, *krs, *kds, *ReturnReal, *ReturnImag, *Return_kr;
   complex<double>        timeShift, SHresult;  
   complex<double>        **OutputArray, **rfArray;

          
   if(nrhs<7) {mexErrMsgIdAndTxt("SOFiA:WGC:notEnoughInputs",
                 "Minimum 7 Inputs required: N, r, ac, FS, NFFT, AZ, EL");}
   
  
  // Get ARGS
   
    N        = (int)mxGetScalar(prhs[0]);
    r        =      mxGetPr(prhs[1]);    
    ac       = (int)mxGetScalar(prhs[2]);
    FS       = (int)mxGetScalar(prhs[3]);
    F_NFFT   = (int)mxGetScalar(prhs[4]);
    az       = (double)mxGetScalar(prhs[5]);
    el       = (double)mxGetScalar(prhs[6]);
       
   
    NFFT = F_NFFT/2+1;
    NMLocatorSize=(N+1)*(N+1); //All n,m up to N included   
    
    // Get ARGS / set defaults
    if(nrhs<14) SegN = N;
    else        SegN = (int)mxGetScalar(prhs[13]);
    
    if(nrhs<13) upperSegLim = NFFT-1;
    else        upperSegLim = (int)mxGetScalar(prhs[12])-1;
    
    if(nrhs<12) lowerSegLim = 0;
    else        lowerSegLim = (int)mxGetScalar(prhs[11])-1;
    

    if (upperSegLim<lowerSegLim) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Invalid segment limits");}
    
    if (upperSegLim>(NFFT-1) || upperSegLim<0) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Invalid segment limits");}
    
    if (lowerSegLim<0 || lowerSegLim>(NFFT-1)) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Invalid segment limits");}
    
    if (SegN>N) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Invalid segment order");}
    
    #ifndef DEBUG
      if(lowerSegLim == 0 && upperSegLim == (NFFT-1)) mexPrintf("SOFiA W/G/C - Wave Generator Core R13-0306\n");
      else if (lowerSegLim == 0 && upperSegLim != (NFFT-1)) mexPrintf("SOFiA W/G/C - Wave Generator Core R11-1220\n      ");
      else mexPrintf("|");
      if (lowerSegLim != 0 && upperSegLim == (NFFT-1)) mexPrintf("\n");
    #endif 
        
    
    if(nrhs<11) ds=1.0;
    else        ds = (double)mxGetScalar(prhs[10]);
        
    if(nrhs<10) wavetype=0;
    else        wavetype = (int)mxGetScalar(prhs[9]);
    
    if (wavetype != 0 && wavetype !=1) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Invalid wavetype. 0:Plane Wave, 1:Spherical Wave");}
        
    if(nrhs<9)  c = 343.0;
    else        c = (double)mxGetScalar(prhs[8]);
    
    if (ac!= 0 && ac!=1 && ac!=2 && ac!=3 && ac!=4) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Invalid array configuration. ac must be in [0,1,2,3,4].");}
    
    if(nrhs<8)  t = 0.0;
    else        t = (double)mxGetScalar(prhs[7]);
    
    if (t*FS > F_NFFT/2) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Delay too large / NFFT too low. Choose t < NFFT/(2*FS), including a NFFT/2 guard interval to avoid cyclic convolution.");}
    
    w      = new double [NFFT];
    k      = new double [NFFT];
    krm    = new double [NFFT]; 
    krs    = new double [NFFT]; 
    kds    = new double [NFFT];
    
    rm     = r[0];
    if      (mxGetM(prhs[1])==2) {rs = r[1]; nor=2;}
    else if (mxGetN(prhs[1])==2) {rs = r[1]; nor=2;}
    else                         {rs = r[0]; nor=1;}
    
    if (nor == 2 && (ac == 0 || ac == 1)) nor=1;
     
    if (wavetype==1 && ds<=rm) {mexErrMsgIdAndTxt("SOFiA:WGC:invalidInputArg",
                 "Invalid source distance. Source must be outside the microphone radius.");}
        
    maxw     = M_PI*FS;
    maxk     = M_PI*FS/c;
    
    iteratew = maxw/(NFFT-1);
    iteratek = maxk/(NFFT-1);
    iteratem = maxk*rm/(NFFT-1);
    iterates = maxk*rs/(NFFT-1);
    iterated = maxk*ds/(NFFT-1);
    
    w[0]   = 0;
    k[0]   = 0;
    krm[0] = 0;
    krs[0] = 0;
    kds[0] = 0;
    
    for(ctr=1; ctr<NFFT; ctr++)
    {
        w[ctr]      = w[ctr-1]+iteratew;
        k[ctr]      = k[ctr-1]+iteratek;
        krm[ctr]    = krm[ctr-1]+iteratem;
        krs[ctr]    = krs[ctr-1]+iterates; 
        kds[ctr]    = kds[ctr-1]+iterated;
    }
    
    krm[0] = krm[1];
    krs[0] = krs[1];
    kds[0] = kds[1];    
    //w[0]   = w[1]
    //k[0]   = k[1];
         
    try
     {
        //Allocate Output Array Size [(N+1)^2]*[FFTBins]
        OutputArray = new complex<double>*[NMLocatorSize];

        for(ctr = 0; ctr < NMLocatorSize; ctr++)		
           {OutputArray[ctr] = new complex<double>[NFFT];}


        //Initializate Output Array
        for(ctr = 0; ctr < NMLocatorSize; ctr++)
           {
            for(ctrb = 0; ctrb < NFFT; ctrb++)
               {
                    OutputArray[ctr][ctrb] = complex<double>(0,0);                
               }
           }        
        
        
        //Allocate bn Array Size [N+1]*[NFFT]
        rfArray = new complex<double>*[N+1];

        for(ctr = 0; ctr < (N+1); ctr++)		
           {rfArray[ctr] = new complex<double>[NFFT];}
        
        
        //Initializate bn Array
        for(ctr = 0; ctr < (N+1); ctr++)
           {
            for(ctrb = 0; ctrb < NFFT; ctrb++)
               {
                    rfArray[ctr][ctrb] = complex<double>(0,0);                
               }
           }   
     }
     catch(...) 
     { mexErrMsgIdAndTxt("SOFiA:WGC:OutArrayAllocation",
                 "Not able to allocate memory for the output Matrix. Maybe to large?");
     }
     
     //Precalculate radial filters 

      if (wavetype == 0)
      {   
        for(f = lowerSegLim; f <= upperSegLim; f++) //PLANE WAVE
        {
           timeShift = exp(complex<double>(0,-1)*complex<double>(w[f]*t,0));
                   
           for(n = 0; n <= SegN; n++)
           {
              rfArray[n][f] = bn(n, krm[f], krs[f], ac)*timeShift;
           }
        }  
      }
     
      else //(wavetype == 1)
      {
        for(f = lowerSegLim; f <= upperSegLim; f++) //SPHERICAL WAVE
        {
           timeShift = exp(complex<double>(0,-1)*complex<double>(w[f]*t,0)); 
            
           for(n = 0; n <= SegN; n++)
           {
               rfArray[n][f] = bn_npf(n, krm[f], krs[f], ac)*complex<double>(4*M_PI,0)*complex<double>(0,-1)*complex<double>(k[f],0)*sphankel(n,kds[f])*timeShift;     
           }
        }          
      }     

            
            
   
//GENERATOR CORE --------------------------------------------------------
         
    ctr=0;
    for(n=0; n<=SegN; n++)
     {      
        for(m=-n; m<=n; m++)
         {               
                 SHresult=conj(spherical_harmonic(n,m,el,az));
                 for(f=lowerSegLim; f<=upperSegLim; f++)
                 {    
                      OutputArray[ctr][f]=SHresult*rfArray[n][f];                      
                 }                 
                 ctr++;                
         }
     } 
    
     
     //Declare return Matrix ARG0
      plhs[0] = mxCreateDoubleMatrix(NMLocatorSize,NFFT,mxCOMPLEX);
      ReturnReal = mxGetPr(plhs[0]);
      ReturnImag = mxGetPi(plhs[0]);
           
     //Return ARG0 Matrix Fill 
      ctrc=0;
      for(ctr=0;ctr<NFFT;ctr++)
         {
         for(ctrb=0;ctrb<NMLocatorSize;ctrb++)
              {
                 ReturnReal[ctrc]=real(OutputArray[ctrb][ctr]);
                 ReturnImag[ctrc]=imag(OutputArray[ctrb][ctr]);
                 ctrc++;
               }
         }
      
      
      //Declare return Matrix ARG1
      w[0]   = 0;
      k[0]   = 0;
      krm[0] = 0;
      krs[0] = 0;
      kds[0] = 0;
      
      plhs[1] = mxCreateDoubleMatrix(nor,NFFT,mxREAL);
      Return_kr = mxGetPr(plhs[1]);
      
      ctrb=0;
      for(ctr=0;ctr<NFFT;ctr++)
         {          
           Return_kr[ctrb]=krm[ctr];
           ctrb++;           
           if (nor == 2)
           {
                Return_kr[ctrb]=krs[ctr];
                ctrb++;
           }           
         }
      
      delete [] w;
      delete [] k;
      delete [] krm;
      delete [] krs;
      delete [] kds;
      for(ctr = 0; ctr < NMLocatorSize; ctr++) delete OutputArray[ctr];
      for(n = 0; n < N+1; n++) delete rfArray[n]; 
}

