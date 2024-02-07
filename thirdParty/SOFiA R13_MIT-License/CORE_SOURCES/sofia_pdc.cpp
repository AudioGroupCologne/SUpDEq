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
% P/D/C Fast pwd|beamfoming core R13-0306
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com                        
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
%  
% Y = ASAR_PDC(N, OmegaL, Pnm, dn, [cn]) 
% ------------------------------------------------------------------------     
% Y      MxN Matrix of the decomposed wavefield 
%        Col - Look Direction as specified in OmegaL
%        Row - kr bins
% ------------------------------------------------------------------------              
% N      Decomposition Order
% 
% OmegaL Look Directions (Vector) 
%        Col - L1, L2, ..., Ln 
%        Row - AZn ELn
% 
% Pnm    Spatial Fourier Coefficients from SOFiA S/T/C
% 
% dn     Modal Array Filters from SOFiA M/F
% 
% cn     (Optional) Weighting Function
%        Can be used for N=0...N weigths:
%        Col - n...N
%        Row - 1
%        Or n(f)...N(f) weigths:
%        Col - n...N
%        Row - kr bins   
%        If cn is not specified a PWD will be done
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

#include "mex.h"
#include <iostream>
#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif

        
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    using namespace std;
    using namespace boost::math;
    int n, m, N, omegactr;       
    int numberOfAngles, NMLocatorSize, NMDeliveredSize, Ndn, Ncn, pwdflag, cnnofreqflag;
    long long unsigned int f, ctr, ctrb, ctrc, FFTBlocklengthPnm, FFTBlocklengthdn, FFTBlocklengthcn;
    double gaincorrection; // Correction is done in OutMTX 
    double *OmegaL, *Elevation, *Azimut, *PnmGETr, *PnmGETi, *dnGETr, *dnGETi, *cnGETr, *cnGETi;
    double *ReturnReal,*ReturnImag; 
    complex<double> Ynm;   
    complex<double> **Pnm, **dn, **cn, **OutputArray;
    
    #ifndef DEBUG
        mexPrintf("SOFiA P/D/C - pwd|beamfoming core R13-0306\n");
    #endif
            
    if(nrhs<4)
    {
        mexErrMsgIdAndTxt("SOFiA:PDC:NotEnoughArgruments",
                 "Arguments missing. At least 4 args required (N, OmegaL, Pnm, dn, [cn]).");
    }    

    //Get ARGS: (N, OmegaL[], Pnm, dn, cn)
    N       = (int)mxGetScalar(prhs[0]); //max Order 
    if(N<0) {N=0;} 

    
    OmegaL  = mxGetPr(prhs[1]);          //Requested_BeamDirections
   
    if (mxGetM(prhs[1])<1 || mxGetN(prhs[1])<2 || mxGetN(prhs[1])>2)
    {
        mexErrMsgIdAndTxt("SOFiA:PDC:InputArgError",
                 "Angle Matrix OmegaL is not valid. Must consist of AZ/EL pairs in one column [AZ1 EL1; AZ2 EL2; ... ; AZn ELn].\nRemember: All angles are in RAD.");
    }    
    
    
            
    numberOfAngles=mxGetM(prhs[1]);

    Azimut    = new double [numberOfAngles]; 
    Elevation = new double [numberOfAngles]; 
    
    for(ctr=0; ctr<(mxGetM(prhs[1])); ctr++)
    {
        Azimut[ctr]    = (double)OmegaL[ctr];
        Elevation[ctr] = (double)OmegaL[ctr+numberOfAngles];        
    }
    
    
    
    if(mxGetN(prhs[2])!=mxGetN(prhs[3]))
    {
        mexErrMsgIdAndTxt("SOFiA:PDC:InputArgError",
                 "FFT Blocksizes of Pnm and dn are not consistent.");
    }    
    
    if(nrhs>4)                          
    {
         if(mxGetN(prhs[2])!=mxGetN(prhs[4]) && mxGetN(prhs[4])!=1)
            {
                mexErrMsgIdAndTxt("SOFiA:PDC:InputArgError",
                         "FFT Blocksize of cn is not consistent to Pnm and dn.");
            }
    }
             
    PnmGETr  = mxGetPr(prhs[2]);          //Spatial Fourier Coefficients
    if(mxIsComplex(prhs[2])==1)
    {    
        PnmGETi  = mxGetPi(prhs[2]);
    }
    NMDeliveredSize = mxGetM(prhs[2]);
    FFTBlocklengthPnm  = mxGetN(prhs[2]);
    
    
    
    dnGETr  = mxGetPr(prhs[3]);           //Modal Filters
    
    if(mxIsComplex(prhs[3])==1)
    {    
        dnGETi  = mxGetPi(prhs[3]);
    }
     
        
    Ndn = mxGetM(prhs[3]);
    FFTBlocklengthdn = mxGetN(prhs[3]);
    
    pwdflag=1;
    if(nrhs>4)                            //If no Coeff given, do PWD
    {
        pwdflag=0;
        cnGETr  = mxGetPr(prhs[4]);           //Beam Shape Coefficients
        if(mxIsComplex(prhs[4])==1)
        { 
            cnGETi  = mxGetPi(prhs[4]);
        }

        Ncn = mxGetM(prhs[4]);
        FFTBlocklengthcn = mxGetN(prhs[4]);
        
        cnnofreqflag=0;
        if(mxGetN(prhs[4])==1){cnnofreqflag=1;}
    }
    
 
    NMLocatorSize=(N+1)*(N+1); //All n,m up to N included, N USER REQUESTED 
    
    if(NMLocatorSize>NMDeliveredSize)
      {
        mexPrintf("\nWARNING: The requested order N=%d cannot be achieved.\n         The Pnm coefficients deliver a maximum of N=%d.\n         Will decompose on maximum available order.\n\n",N, (int)sqrt((float)NMDeliveredSize)-1);
        N=(int)sqrt((float)NMDeliveredSize)-1;
       }
    
    if(N>(Ndn-1))
      {
        mexPrintf("\nWARNING: The requested order N=%d cannot be achieved.\n         The dn filters deliver a maximum of N=%d.\n         Will decompose on maximum available order.\n\n",N, (int)Ndn-1);
        N=(int)Ndn-1;
       }
    if(pwdflag==0)
    {
        if(N>(Ncn-1))
          {
            mexPrintf("\nWARNING: The requested order N=%d cannot be achieved.\n         The cn beam filters are defined for a maximum of N=%d.\n         Will decompose on maximum available order.\n\n",N, (int)Ncn-1);
            N=(int)Ncn-1;
           }
     }
    
    gaincorrection=(4*M_PI)/pow(((double)N+1),2);
    
    try
    {    
        // ------------------------------------------------- Pnm
        //Allocate Pnm Array Size [(N+1)^2]*[FFTBins]           
        Pnm = new complex<double>*[NMDeliveredSize];
        for(ctr = 0; ctr < NMDeliveredSize; ctr++)		
           {Pnm[ctr] = new complex<double>[FFTBlocklengthPnm];}
        
        if(mxIsComplex(prhs[2])==1)
        { 
        //Write Pnm Data to Array: Pnm[nm][f]
        for(ctr = 0; ctr < FFTBlocklengthPnm; ctr++)
           {
            for(ctrb = 0; ctrb < NMDeliveredSize; ctrb++)
               {
                    Pnm[ctrb][ctr] = complex<double>(PnmGETr[ctrb+ctr*NMDeliveredSize],PnmGETi[ctrb+ctr*NMDeliveredSize]);            
               }
           }  
        }//END IF COMPLEX
        else
             { 
        //Write Pnm Data to Array: Pnm[nm][f]
        for(ctr = 0; ctr < FFTBlocklengthPnm; ctr++)
           {
            for(ctrb = 0; ctrb < NMDeliveredSize; ctrb++)
               {
                    Pnm[ctrb][ctr] = complex<double>(PnmGETr[ctrb+ctr*NMDeliveredSize],0);            
               }
           }  
        }//END ELSE COMPLEX
        
        
        // ------------------------------------------------- dn
        //Allocate dn Array Size 
        dn = new complex<double>*[Ndn];
        for(ctr = 0; ctr < Ndn; ctr++)		
           {dn[ctr] = new complex<double>[FFTBlocklengthdn];}
        
        if(mxIsComplex(prhs[3])==1)
        {
            //Write dn Data to Array: dn[n][f]
            for(ctr = 0; ctr < FFTBlocklengthdn; ctr++)
               {
                for(ctrb = 0; ctrb < Ndn; ctrb++)
                   {
                        dn[ctrb][ctr] = complex<double>(dnGETr[ctrb+ctr*Ndn],dnGETi[ctrb+ctr*Ndn]);            
                   }
               }
        }
        else
            {
            //Write dn Data to Array: dn[n][f]
            for(ctr = 0; ctr < FFTBlocklengthdn; ctr++)
               {
                for(ctrb = 0; ctrb < Ndn; ctrb++)
                   {
                        dn[ctrb][ctr] = complex<double>(dnGETr[ctrb+ctr*Ndn],0);            
                   }
               }
        }
        
        // ------------------------------------------------- cn 
        
        if(pwdflag==0) 
        {
              //Allocate dn Array Size 
              cn = new complex<double>*[Ncn];
              for(ctr = 0; ctr < Ncn; ctr++)		
                 {cn[ctr] = new complex<double>[FFTBlocklengthPnm];} 
                
              if(mxIsComplex(prhs[4])==1)
              { 
                  if(cnnofreqflag==0)
                  {
                  //Write dn Data to Array: dn[n][f]
                  for(ctr = 0; ctr < FFTBlocklengthPnm; ctr++) //Yes, Pnm is ok!
                     {
                      for(ctrb = 0; ctrb < Ncn; ctrb++)
                         {
                            cn[ctrb][ctr] = complex<double>(cnGETr[ctrb+ctr*Ncn],cnGETi[ctrb+ctr*Ncn]);            
                         }
                     }
                  }
                  else //Expand single filters to all FFT Bins
                  {
                     for(ctr = 0; ctr < FFTBlocklengthPnm; ctr++) //Yes, Pnm is ok! 
                     {
                      for(ctrb = 0; ctrb < Ncn; ctrb++)
                         {
                            cn[ctrb][ctr] = complex<double>(cnGETr[ctrb],cnGETi[ctrb]);            
                         }
                     }                    
                  }    
              }//END IF COMPLEX
              else
                  { 
                  if(cnnofreqflag==0)
                  {
                  //Write dn Data to Array: dn[n][f]
                  for(ctr = 0; ctr < FFTBlocklengthPnm; ctr++) //Yes, Pnm is ok!
                     {
                      for(ctrb = 0; ctrb < Ncn; ctrb++)
                         {
                            cn[ctrb][ctr] = complex<double>(cnGETr[ctrb+ctr*Ncn],0);            
                         }
                     }
                  }
                  else //Expand single filters to all FFT Bins
                  {
                     for(ctr = 0; ctr < FFTBlocklengthPnm; ctr++) //Yes, Pnm is ok! 
                     {
                      for(ctrb = 0; ctrb < Ncn; ctrb++)
                         {
                            cn[ctrb][ctr] = complex<double>(cnGETr[ctrb],0);            
                         }
                     }                    
                  }    
              }//END ELSE COMPLEX
        }
        
        
        
        
        
        
        // ------------------------------------------------- OutputArray
        //Allocate Output Array Size [numberOfAngles]*[FFTBins]
        OutputArray = new complex<double>*[numberOfAngles];

        for(ctr = 0; ctr < numberOfAngles; ctr++)		
           {OutputArray[ctr] = new complex<double>[FFTBlocklengthPnm];}

        //Initializate Output Array
        for(ctr = 0; ctr < numberOfAngles; ctr++)
           {
            for(ctrb = 0; ctrb < FFTBlocklengthPnm; ctrb++)
               {
                    OutputArray[ctr][ctrb] = complex<double>(0,0);                
               }
           }
        
        }//TRY allocation of memory      
     catch(...) 
     { mexErrMsgIdAndTxt("SOFiA:PDC:OutArrayAllocation",
                 "Not able to allocate enough memory.");
     }
           
  
    
    // ----------------------------------------------------------- PWD-CORE

    ctr=0;
    Ynm=complex<double>(1,0);
    if(pwdflag==1) 
    {
        for(n=0; n<=N; n++)
        {
           for(m=-n; m<=n; m++)
           { 
              for(omegactr=0; omegactr<numberOfAngles; omegactr++) 
              {
                   Ynm=spherical_harmonic(n,m,Elevation[omegactr],Azimut[omegactr]);   

                  for(f=0; f<FFTBlocklengthPnm; f++)
                  {
                      OutputArray[omegactr][f]+=Ynm*Pnm[ctr][f]*dn[n][f];
                  }              
              }
              ctr++;
           }
        }
    }//END PWDFLAG
     // --------------------------------------------------BEAMFORMING-CORE
    
    else
    {

        for(n=0; n<=N; n++)
        {
           for(m=-n; m<=n; m++)
           { 
              for(omegactr=0; omegactr<numberOfAngles; omegactr++) 
              {
                   Ynm=spherical_harmonic(n,m,Elevation[omegactr],Azimut[omegactr]);   

                  for(f=0; f<FFTBlocklengthPnm; f++)
                  {
                      OutputArray[omegactr][f]+=Ynm*Pnm[ctr][f]*dn[n][f]*cn[n][f];
                  }              
              }
              ctr++;
           }
        }
    } //ELSE PWD
          
           
     //Declare return ARG Matrix 
     plhs[0] = mxCreateDoubleMatrix(numberOfAngles,FFTBlocklengthPnm,mxCOMPLEX);
     ReturnReal = mxGetPr(plhs[0]);
     ReturnImag = mxGetPi(plhs[0]);
     
     //Return ARG Matrix Fill 
     ctrc=0;     
     for(ctr=0;ctr<FFTBlocklengthPnm;ctr++)
        {
        for(ctrb=0;ctrb<numberOfAngles;ctrb++)
             {
                ReturnReal[ctrc]=real(OutputArray[ctrb][ctr])*gaincorrection;
                ReturnImag[ctrc]=imag(OutputArray[ctrb][ctr])*gaincorrection;
                ctrc++;
              }
        }         
  
      delete[] Elevation;
      delete[] Azimut;
      for(ctr = 0; ctr < NMDeliveredSize; ctr++) delete Pnm[ctr];
      for(ctr = 0; ctr < Ndn; ctr++)             delete dn[ctr];
      if (pwdflag==0) 
        { for(ctr = 0; ctr < Ncn; ctr++)         delete cn[ctr];}
      for(ctr = 0; ctr < numberOfAngles; ctr++)  delete OutputArray[ctr];      
     
}

