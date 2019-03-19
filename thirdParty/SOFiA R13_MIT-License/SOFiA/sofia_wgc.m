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
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
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
%          0  Open Sphere with p Transducers (NO plc!)
%          1  Open Sphere with pGrad Transducers
%          2  Rigid Sphere with p Transducers
%          3  Rigid Sphere with pGrad Transducers (Thx to Nils Peters!)
%          4  Dual Open Sphere with p Transducers (Thx to Nils Peters!)
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

