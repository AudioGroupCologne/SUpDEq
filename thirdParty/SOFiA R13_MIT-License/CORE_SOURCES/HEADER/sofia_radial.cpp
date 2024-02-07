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
% sofia_radial.cpp (Internal C++ Header)
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%                        and Nils Peters, nils 'AT' icsi.berkeley.edu 
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
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

#include "sofia_radial.h"

#ifndef M_PI
	#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#endif

// spherical functions     
std::complex<double> spbessel(int n, double kr) 
{
    std::complex<double> jn, prefactor; 
    double np;
    np=(double)n+0.5;
    prefactor=std::complex<double>(sqrt(M_PI/(2*kr)),0);
    jn = prefactor*std::complex<double>(boost::math::cyl_bessel_j(np,kr),0);  
    return jn;
}

std::complex<double> dspbessel(int n, double kr) 
{
    std::complex<double> djn;    
    djn = spbessel(n-1,kr)-spbessel(n,kr)*(std::complex<double>(n+1,0)/std::complex<double>(kr,0));   
    return djn;
}

std::complex<double> spneumann(int n, double kr) 
{
    std::complex<double> yn, prefactor;
    double np;
    np=(double)n+0.5;
    prefactor=std::complex<double>(sqrt(M_PI/(2*kr)),0);
    yn = prefactor*std::complex<double>(boost::math::cyl_neumann(np,kr),0);  
    return yn;
}   

std::complex<double> dspneumann(int n, double kr) 
{
    std::complex<double> dyn;      
    dyn = spneumann(n-1,kr)-spneumann(n,kr)*(std::complex<double>(n+1,0)/std::complex<double>(kr,0));  
    return dyn;
}   

std::complex<double> sphankel(int n, double kr) 
{
    std::complex<double> hn;  
    hn = spbessel(n,kr)-std::complex<double>(0,1)*spneumann(n, kr);
    return hn;
}   

std::complex<double> dsphankel(int n, double kr) 
{
    std::complex<double> dhn;   
    dhn = std::complex<double>(0.5,0)*(sphankel(n-1,kr)-sphankel(n+1,kr)-sphankel(n,kr)/kr);
    return dhn;
}  


//--- bn 

std::complex<double> bn_openP(int n, double krm)
{
    using namespace std;
    complex<double> bn;            
    bn  = spbessel(n,krm);
    return bn;
}

std::complex<double> bn_openPG(int n, double krm)
{
    using namespace std;
    complex<double> jn, diffjn, bn;     
    bn = complex<double>(0.5,0)*(spbessel(n,krm)-complex<double>(0,1)*dspbessel(n,krm));
    //                   ! 1/2
    return bn;
}

std::complex<double> bn_rigidP(int n, double krm, double krs)
{
    using namespace std;    
    complex<double> bn;
    bn = spbessel(n,krm)-(dspbessel(n,krs)/dsphankel(n,krs))*sphankel(n,krm);
    return bn;
}

std::complex<double> bn_rigidPG(int n, double krm, double krs) 

/* Rerence for Filter design for rigid sphere with cardioid microphones:
   P. Plessas, F. Zotter: Microphone arrays around rigid spheres for spatial recording and holography, DAGA 2010
   krm: for mic radius, krs: for sphere radius
   Implementation by Nils Peters, November 2011 */
{
    using namespace std;    
    complex<double> bn;
    bn = (spbessel(n,krm)-complex<double>(0,1)*dspbessel(n,krm) + (complex<double>(0,1)*dsphankel(n,krm)-sphankel(n,krm)) * (dspbessel(n,krs)/dsphankel(n,krs)));
    return bn;
}

std::complex<double> bn_dualOpenP(int n, double kr1, double kr2) 

/* Reference: Rafaely et al, 
   High-resolution plane-wave decomposition in an auditorium using a dual-radius scanning spherical microphone array
   JASA, 122(5), 2007

   kr1, kr2 are the kr values of the two different microphone spheres
   Implementation by Nils Peters, November 2011*/

{
    using namespace std;    
    complex<double> bn1, bn2;
    
    bn1 = bn_openP(n,kr1);
    bn2 = bn_openP(n,kr2);
    
    if (abs(bn1) >= abs(bn2))
        return bn1;
    else
        return bn2;	
}

std::complex<double> bn_npf(int n, double krm, double krs, int ac)
{
    std::complex<double> bnval;
    switch(ac)
                    {
                     case 0: bnval = bn_openP(n, krm);
                             break;

                     case 1: bnval = bn_openPG(n, krm);
                             break;

                     case 2: bnval = bn_rigidP(n, krm, krs);
                             break; 

					 case 3: bnval = bn_rigidPG(n, krm, krs);
		                     break;
                             
                     case 4: bnval = bn_dualOpenP(n, krm, krs);
		                     break;  	                    
                    }     
    
    return bnval;
}


std::complex<double> bn(int n, double krm, double krs, int ac)
{    
    std::complex<double> bnval;
    bnval = bn_npf(n, krm, krs, ac) * std::complex<double>(4*M_PI,0) * pow((std::complex<double>(0,1)),n);    
    return bnval;
}




