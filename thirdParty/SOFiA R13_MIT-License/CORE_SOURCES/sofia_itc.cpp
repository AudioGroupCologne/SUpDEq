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
% I/T/C Fast Inverse spatial Fourier Transform Core R13-0306
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com                        
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
% 
% p = sofia_itc(Pnm, angles, [N])
% ------------------------------------------------------------------------
% p      sound pressures (complex data)  
%        Columns : FFT bins
%        Rows    : angles   
% ------------------------------------------------------------------------
% Pnm    spatial Fourier coefficients (e.g. from SOFiA S/T/C)
%        Columns : FFT bins 
%        Rows    : nm coeff  
% 
% angles target angles [AZ1 EL1; AZ2 EL2; ... AZn ELn]
%        Columns : Angle Number 1...n
%        Rows    : AZ EL    
%          
% [N]     *** Optional: Maximum transform order 
%            If not specified the highest order available included in
%            the Pnm coefficients will be taken.
% 
% This is a pure ISFT core that does not involve extrapolation. 
% (=The pressures are referred to the original radius)
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
   double *grid, *PnmDataReal, *PnmDataImag;   
   int    N, Nmax, PnmDataLength;
   int    n, m, j;
   long long unsigned int ctr, ctrb, ctrc, f;
   long unsigned int FFTBlocklength;
   int    numberOfAngles, numberOfAngleInfos;
   double *AzimuthAngles,*ElevationAngles;
   double *ReturnReal,*ReturnImag; 
   complex<double> **OutputArray;
   complex<double> SHresult;
   
   #ifndef DEBUG
   mexPrintf("SOFiA I/T/C - Inverse spatial Transform Core R13-0306\n");
   #endif
   if(nrhs<2) {mexErrMsgIdAndTxt("SOFiA:ITC:notEnoughInputs",
                 "Minimum 2 Inputs required: p = sofia_itc(Pnm coefficients, angles, [N])");}
   
   //Get ARGS
         
    PnmDataReal = mxGetPr(prhs[0]);
    PnmDataImag = mxGetPi(prhs[0]);
    
    if(PnmDataImag==NULL)
        {mexErrMsgIdAndTxt("SOFiA:ITC:InputArgError",
                 "PnmData: Complex Input Data expected.");}
    
    FFTBlocklength=mxGetN(prhs[0]);            
      
    grid  = mxGetPr(prhs[1]);   
    numberOfAngles = mxGetM(prhs[1]);
    numberOfAngleInfos  = mxGetN(prhs[1]);
    
    if(numberOfAngleInfos<2)
        {mexErrMsgIdAndTxt("SOFiA:ITC:InputArgError",
         "Error: Delivered angles are not valid. Must consist of [AZ1 EL1; AZ2 EL2; ...; AZn ELn] pairs.");}
    
    PnmDataLength = mxGetM(prhs[0]);
    
    if(nrhs>2){
        N    = (int)mxGetScalar(prhs[2]); 
        Nmax = (int)sqrt((float)PnmDataLength)-1;
        if (N>Nmax)
            {mexErrMsgIdAndTxt("SOFiA:ITC:MaximumOrderExceed",
                 "Requested order too high: Maximum available order in Pnm coefficients exceeded.");}        
    }
    else {                           
        N = (int)sqrt((float)PnmDataLength)-1; 
        Nmax = N; 
    }                   
   
//     mexPrintf("Delivered FFT Taps      %ld\n",FFTBlocklength);
//     mexPrintf("Spatial Positions       %d\n",numberOfAngles);
//     mexPrintf("Transform Order         %d",N);   
    
    
    //Create Dynamic Arrays (angles)
    AzimuthAngles   = new double [numberOfAngles]; 
    ElevationAngles = new double [numberOfAngles]; 

    
    //Fill Cells GRID
    for(ctr=0; ctr<numberOfAngles; ctr++)
        {
            AzimuthAngles[ctr]   = grid[ctr];
            ElevationAngles[ctr] = grid[ctr+numberOfAngles];                           
        }        
       
    try
     {
        //Allocate Output Array Size [numberOfAngles]*[FFTBins]
        OutputArray = new complex<double>*[numberOfAngles];

        for(ctr = 0; ctr < numberOfAngles; ctr++)		
           {OutputArray[ctr] = new complex<double>[FFTBlocklength];}


        //Initializate Output Array
        for(ctr = 0; ctr < numberOfAngles; ctr++)
           {
            for(ctrb = 0; ctrb < FFTBlocklength; ctrb++)
               {
                    OutputArray[ctr][ctrb] = complex<double>(0,0);                
               }
           }        
     }
     catch(...) 
     { mexErrMsgIdAndTxt("SOFiA:ITC:OutArrayAllocation",
                 "Not able to allocate memory for the output Matrix. Maybe to large?");
     }
             
   
    //SHT TRANSFORM CORE -------------------------------------------------
        
    ctr=0;
   
     for(n=0; n<=N; n++)
     {
         for(m=-n; m<=n; m++)
         {               
             for(j=0; j<numberOfAngles; j++)
             {
                 //SHresult=complex<double>(1/(4*M_PI),0);                 
                 SHresult = spherical_harmonic(n,m,ElevationAngles[j], AzimuthAngles[j]);
                 for(f=0; f<FFTBlocklength; f++)
                 {                   
                     OutputArray[j][f] = OutputArray[j][f] + SHresult * complex<double>(PnmDataReal[ctr+f*PnmDataLength],PnmDataImag[ctr+f*PnmDataLength]);                  
                 }                 
             }
             ctr++;                
         }
     }                
     
        
    
     //Declare return ARG Matrix 
     plhs[0] = mxCreateDoubleMatrix(numberOfAngles,FFTBlocklength,mxCOMPLEX);
     ReturnReal = mxGetPr(plhs[0]);
     ReturnImag = mxGetPi(plhs[0]);
     
     //Return ARG Matrix Fill 
     ctrc=0;     
     for(ctr=0;ctr<FFTBlocklength;ctr++)
        {
        for(ctrb=0;ctrb<numberOfAngles;ctrb++)
             {
                ReturnReal[ctrc]=real(OutputArray[ctrb][ctr]);
                ReturnImag[ctrc]=imag(OutputArray[ctrb][ctr]);
                ctrc++;
              }
        }         

      delete [] AzimuthAngles;
      delete [] ElevationAngles;
      for(ctr = 0; ctr < numberOfAngles; ctr++) delete OutputArray[ctr];      
}


