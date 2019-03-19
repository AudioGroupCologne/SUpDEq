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
% sofia_radial.h (Internal C++ Header)
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



#ifndef SOFIA_RADIAL_H
#define SOFIA_RADIAL_H

#include <boost/math/special_functions/bessel.hpp>

std::complex<double> spbessel(int n, double kr); 
std::complex<double> dspbessel(int n, double kr); 
std::complex<double> spneumann(int n, double kr); 
std::complex<double> dspneumann(int n, double kr); 
std::complex<double> sphankel(int n, double kr); 
std::complex<double> dsphankel(int n, double kr); 
std::complex<double> bn_openP(int n, double krm);
std::complex<double> bn_openPG(int n, double krm);
std::complex<double> bn_rigidP(int n, double krm, double krs);
std::complex<double> bn_rigidPG(int n, double krm, double krs); 
std::complex<double> bn_dualOpenP(int n, double kr1, double kr2);
std::complex<double> bn(int n, double krm, double krs, int ac);
std::complex<double> bn_npf(int n, double krm, double krs, int ac);

#endif