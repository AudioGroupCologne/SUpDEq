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
