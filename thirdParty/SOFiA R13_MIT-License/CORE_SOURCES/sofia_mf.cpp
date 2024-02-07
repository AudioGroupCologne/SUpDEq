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
% M/F Modal radial filters R13-0306
%     Soft amplification limiting
%     On-axis powerloss compensation with
%     N0plc to N0 interpolation.
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%                        and Nils Peters, nils 'AT' icsi.berkeley.edu 
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
% 
% [dn, beam] = SOFIA_MF(N, kr, ac, [a_max], [plc], [fadeover])
% ------------------------------------------------------------------------   
% dn          Vector of modal 0-N frequency domain filters
% beam        Expected free field On-Axis kr-response 
% ------------------------------------------------------------------------
% N           Maximum Order
% kr          Vector or Matrix of kr values
%             First Row   (M=1) N: kr values Microphone Radius
%             Second Row  (M=2) N: kr values Sphere/Microphone2 Radius 
%             [kr_mic;kr_sphere] for Rigid/Dual Sphere Configurations
%             ! If only one kr-vector is given using a Rigid/Dual Sphere  
%             Configuration: kr_sphere = kr_mic 
% ac          Array Configuration: 
%             0  Open Sphere with pressure Transducers (NO plc!)
%             1  Open Sphere with cardioid Transducers
%             2  Rigid Sphere with pressure Transducers
%             3  Rigid Sphere with cardioid Transducers (Thx to Nils Peters!)
%             4  Dual Open Sphere with pressure Transducers (Thx to Nils Peters!)
% a_max       Maximum modal amplification limit in [dB]
% plc         OnAxis powerloss-compensation: 
%             0  Off
%             1  Full kr-spectrum plc
%             2  Low kr only -> set fadeover             
% fadeover    Number of kr values to fade over +/- around min-distance 
%             gap of powerloss compensated filter and normal N0 filters.
%             0 = auto fadeover
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
%     Nils Peters                  nils 'at' icsi.berkeley.edu
% 
**********************************************************************/

#include <mex.h>
#include <complex>
#include <math.h>
#include <sofia_radial.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
   using namespace std;
   using namespace boost::math;

   int N, plc, ac, noplcflag, limiteronflag, normalizebeam, zeroflag;
   long long unsigned int ctr, ctrb, ctrc, numberOfkrValues, locatemindistance, fadeover;   
   long long signed int ctrd;
   double *kr, *krm, *krs;        
   double *ReturnReal,*ReturnImag,*ReturnRealB,*ReturnImagB;
   double a_max, a_maxdB, amplicalc, mindistance, mix, filtergap;
   complex<double> bnval; 
   complex<double> *xi, *beamresponse;
   complex<double> **OutputArray;   

   #ifndef DEBUG
    mexPrintf("SOFiA M/F - Modal radial filters R13-0306\n");
   #endif
   
   if(nrhs<3) {mexErrMsgIdAndTxt("SOFiA:MF:notEnoughInputs",
                 "Minimum 3 Inputs required: N, kr, ac, [a_max], [plc], [fadeover]");}
    
    //Get first 3 ARGS
    N           =           (int)mxGetScalar(prhs[0]);        
    kr          =                    mxGetPr(prhs[1]);
    ac          =           (int)mxGetScalar(prhs[2]);
    
    if(nrhs<4) 
    {
        a_max=1;
        a_maxdB=0;
        limiteronflag=0;
    }
    else
    {
        a_maxdB =(double)mxGetScalar(prhs[3]);        
        limiteronflag=1; 
        a_max=pow(10,(a_maxdB/(double)20)); //dB-Value to double skalar
    }
   
    if(nrhs<5) {plc=0;}
    else
    {
       plc = (int)mxGetScalar(prhs[4]);
    }
   
    if(nrhs<6) {fadeover=0;}
    else
    {
       fadeover = (long long int)mxGetScalar(prhs[5]);
    }
  
      
    numberOfkrValues=mxGetN(prhs[1]);
      
    
    if(mxGetM(prhs[1])>2)
    {
        {mexErrMsgIdAndTxt("SOFiA:MF:InputArgError",
                 "Vector size not valid: NxM with N=1 or 2 and M=[...] expected.");}        
    }
    
    if(N<0 || N>40)
    {
        {mexErrMsgIdAndTxt("SOFiA:MF:InputArgError",
                 "N: Order not valid.");}        
    }
    if(a_maxdB<-20.0f)
    {
        {mexErrMsgIdAndTxt("SOFiA:MF:InputArgError",
                 "a_max: Amplification limit too low.");}        
    }
    
    if(ac!=0 && ac!=1 && ac!=2 && ac!=3 && ac!=4)
    {
        {mexErrMsgIdAndTxt("SOFiA:MF:InputArgError",
                 "ac: Array configuration not valid. ac=[0,1,2,3,4]");}        
    }
    
    if(plc!=0 && plc!=1 && plc!=2)
    {
        {mexErrMsgIdAndTxt("SOFiA:MF:InputArgError",
                 "plc: Choice not valid. plc=[0,1,2]");}        
    }
    
    try
     {
        //Allocate Output Array Size [N]*[size(kr)]
        OutputArray = new complex<double>*[N+1];

        for(ctr = 0; ctr <= N; ctr++)		
           {OutputArray[ctr] = new complex<double>[numberOfkrValues];}


        //Initializate Output Array
        for(ctr = 0; ctr <= N; ctr++)
           {
            for(ctrb = 0; ctrb < numberOfkrValues; ctrb++)
               {
                    OutputArray[ctr][ctrb] = complex<double>(0,0);                
               }
           }
        
     }
     catch(...) 
     { mexErrMsgIdAndTxt("SOFiA:MF:OutArrayAllocation",
                 "Not able to allocate memory for the output Matrix. Maybe to large?");
     }
    
    

    //Make Dynamic Arrays 
    krm    = new double [mxGetN(prhs[1])]; 
    krs    = new double [mxGetN(prhs[1])]; 

    
    for(ctr=0; ctr<(mxGetN(prhs[1])); ctr++)// ctr=ctr+mxGetM(prhs[1]))
    {
        krm[ctr]    = (double)kr[ctr*mxGetM(prhs[1])];
        krs[ctr]    = (double)kr[ctr*mxGetM(prhs[1])+mxGetM(prhs[1])-1];        
    }
        
  
     //Check for Zero-Elements in kr-Vector
     for(ctr = 0; ctr < numberOfkrValues; ctr++)
     {
        zeroflag=0;
        if (krm[ctr]<=0)
        {
            //mexPrintf("\nWarning: kr (mic) contains zero element or negative value.\n");
            zeroflag=1;
            
            for(ctrb = ctr; ctrb < numberOfkrValues; ctrb++) //Try to repair with next valid element
            {
                if(krm[ctrb]>0)
                {
                    krm[ctr]=krm[ctrb];
                    //mexPrintf("         Replaced by next valid kr element [%f].\n", krm[ctrb]);
                    zeroflag=0;
                    break;
                }               
            }    
            if(zeroflag==1)
              {                
                 mexErrMsgIdAndTxt("SOFiA:MF:DataNotValid",
                "SOFiA MF Error: kr vector is not valid.");
              }
           
        }
        
        zeroflag=0;
        if (krs[ctr]<=0)
        {   
            //mexPrintf("\nWarning: kr (sphere) contains zero element or negative value.\n");
            zeroflag=1;
            
            for(ctrb = ctr; ctrb < numberOfkrValues; ctrb++) //Try to repair with next valid element
            {
                if(krs[ctrb]>0)
                {
                    krs[ctr]=krs[ctrb];
                    //mexPrintf("         Replaced by next valid kr element [%f].\n", krs[ctrb]);
                    zeroflag=0;
                    break;
                }               
            }    
            if(zeroflag==1)
              {                
                 mexErrMsgIdAndTxt("SOFiA:MF:DataNotValid",
                "SOFiA MF Error: kr vector is not valid.");
              }
           
        }
     }
       
    
    //DO bn-Filter Calculation    
     for(ctr = 0; ctr <= N; ctr++)
           {

            for(ctrb = 0; ctrb < numberOfkrValues; ctrb++)
               {      
                    bnval=bn(ctr, krm[ctrb], krs[ctrb], ac);
                    if(limiteronflag==1)
                    {
                        amplicalc=(2*a_max/M_PI)*abs(bnval)*atan(M_PI/(2*a_max*abs(bnval))); 
                    }
                    else {amplicalc=1;}
                    OutputArray[ctr][ctrb]=complex<double>(amplicalc,0)/bnval;                   
               }
           }   
    
     
     if(numberOfkrValues<32 && plc!=0)
      {mexPrintf("\nWARNING: Not enough kr values for PLC fading. PLC disabled.\n");
          plc=0;}
    
    
     //POWERLOSS COMPENSATION FILTER
     noplcflag=0;
     if(plc!=0)
     {        
         xi= new complex<double>[numberOfkrValues]; //Array for xi-Term

         for(ctr=0; ctr<numberOfkrValues; ctr++)
         {
            xi[ctr]=0;               
             for(ctrb=0; ctrb<=N; ctrb++)
             {
                 xi[ctr]+=complex<double>(2*ctrb+1,0)*(complex<double>(1,0)-OutputArray[ctrb][ctr]*bn(ctrb,krm[ctr],krs[ctr],ac));           
             }
            xi[ctr]*=complex<double>(1,0)/bn(0,krm[ctr],krs[ctr],ac); 
            xi[ctr]+=OutputArray[0][ctr];            
         }
     }//plc!=0 flag
     
     
         
     if(plc==1) // -------------------- low kr only
     {       
         
         //Find minimum distance
         mindistance=(double)abs(OutputArray[0][0]-xi[0]);
         locatemindistance=0;
         
         for(ctr=0; ctr<numberOfkrValues; ctr++)
         {
             if((double)abs(OutputArray[0][ctr]-xi[ctr]) < mindistance)
             {
                 mindistance=(double)abs(OutputArray[0][ctr]-xi[ctr]);
                 locatemindistance=ctr;                 
             }                  
    
         }
         
         filtergap=20*log10(1/abs(OutputArray[0][locatemindistance]/xi[locatemindistance]));
         mexPrintf("\nFilter fade gap: %4.2f dB",filtergap);
                          
         if (filtergap>(double)20.0f || filtergap<(double)-20.0f)
         {mexPrintf("\nWARNING: Filtergap is too large. Nonsense filter expected.\nNo powerloss compensation filter applied.");
          noplcflag=1;}
                 
         
         //Overwrite order 0 filter in output matrix
         if (noplcflag==0)
         {
             if (filtergap>(double)5.0f)
                {mexPrintf("\nWARNING: Filtergap is very large.");}
             
             if(fadeover==0) //Auto fadeover size
             {
                 fadeover=(long long int)(numberOfkrValues/100);
                 if (a_maxdB>0)
                 {
                     fadeover=fadeover/(long long int)ceil(a_max/4);
                 }

                if((fadeover>locatemindistance) || ((locatemindistance+fadeover)>numberOfkrValues))
                     {
                     if ((long long int)(locatemindistance-fadeover)<(long long int)(numberOfkrValues-(locatemindistance+fadeover)))
                      {
                         fadeover=locatemindistance;
                      }
                       else
                       {
                         fadeover=numberOfkrValues-locatemindistance;
                      }                  
                      }
                  mexPrintf("\nAuto filter size: %d Taps.", fadeover);
             }

             if((fadeover>locatemindistance) || ((locatemindistance+fadeover)>numberOfkrValues))
             {
               if ((long long int)(locatemindistance-fadeover)<(long long int)(numberOfkrValues-(locatemindistance+fadeover)))
                {
                    fadeover=locatemindistance;
                }
                 else
                 {
                    fadeover=numberOfkrValues-locatemindistance;
                }
                mexPrintf("\nWARNING: Filter fade size too high. Reduced to %d Taps.", fadeover);
             }

             mix=0;   
             for(ctr=0; ctr<=locatemindistance-fadeover; ctr++)
             {
                OutputArray[0][ctr]=xi[ctr];            
             }

             for(ctr=locatemindistance-fadeover+1; ctr<=locatemindistance+fadeover-1; ctr++)
                 {
                    mix+=1/(double)(2*fadeover);
                    OutputArray[0][ctr]=OutputArray[0][ctr]*complex<double>(mix,0)+xi[ctr]*complex<double>(1.0f-mix,0);
                }
             
             
             }//NoPLCFlagEnd
     }//PLC=1 low kr only
          
     if(plc==2) // full spectrum
     {
             for(ctr=0; ctr<numberOfkrValues; ctr++)
             {
                OutputArray[0][ctr]=xi[ctr];            
             }
     }  
     
     
     //Ouput kr-Response     
     beamresponse= new complex<double>[numberOfkrValues]; //Array for kr-response-Term

     normalizebeam=((N+1)*(N+1));
     
     for(ctr=0; ctr<numberOfkrValues; ctr++)
     {
         beamresponse[ctr]=0;                            //Initialize         
         for(ctrb=0; ctrb<=N; ctrb++)                   //ctrb=n;
         {
             for(ctrc=0; ctrc<=(2*ctrb); ctrc++)        //ctrc=m;
             {                 
                 beamresponse[ctr]=beamresponse[ctr]+bn(ctrb,krm[ctr],krs[ctr],ac)*OutputArray[ctrb][ctr];                 
             }        
         }   
         beamresponse[ctr]/=complex<double>((double)normalizebeam,0);        
     }

     
     
     //Declare return ARG Matrix 
     plhs[0] = mxCreateDoubleMatrix(N+1,numberOfkrValues,mxCOMPLEX);
     plhs[1] = mxCreateDoubleMatrix(1,numberOfkrValues,mxCOMPLEX);
     
     ReturnReal = mxGetPr(plhs[0]);
     ReturnImag = mxGetPi(plhs[0]);
     
     //Return ARG Matrix Fill 
     ctrc=0;     
     for(ctr=0;ctr<numberOfkrValues;ctr++)
        {
        for(ctrb=0;ctrb<=N;ctrb++)
             {
                ReturnReal[ctrc]=real(OutputArray[ctrb][ctr]);
                ReturnImag[ctrc]=imag(OutputArray[ctrb][ctr]);
                ctrc++;
              }
        }      
  
    //Return Beamresponse     
     ReturnRealB = mxGetPr(plhs[1]);
     ReturnImagB = mxGetPi(plhs[1]);
     
     for(ctr=0;ctr<numberOfkrValues;ctr++)
         {         
                ReturnRealB[ctr]=real(beamresponse[ctr]);
                ReturnImagB[ctr]=imag(beamresponse[ctr]);             
         }   
            
      if(plc!=0)
      {delete [] xi;}
      delete [] beamresponse;
      for(ctr = 0; ctr < (N+1); ctr++) delete OutputArray[ctr]; 
     
}


