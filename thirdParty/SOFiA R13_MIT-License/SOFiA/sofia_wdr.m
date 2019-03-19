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
% SOFiA W/D/R Wigner-D Rotation R13-0306
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% PnmRot = sofia_wdr(Pnm, xAngle, yAngle, zAngle, deg)
% ------------------------------------------------------------------------
% PnmRot   Output: Rotated spatial Fourier coefficients
%
% Pnm      Input: Spatial Fourier coefficients
% xAngle   Rotation angle around the x-Axis 
% yAngle   Rotation angle around the y-Axis 
% zAngle   Rotation angle around the z-Axis 
% deg      true  (1) - Angles in degree
%          false (0) - Angles in rad (#default if not set)
%
% This routine enables the 3D rotation of the spatial Fourier coefficients
% (Pnm) in the spherical wave spectrum domain.
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
% The external routine for the wigner-D matrix calculation: wignerd.p is
% taken from the free EasySpin Toolbox (http://www.easyspin.org/)
% -----------------------------------------------------------------------
% Stefan Stoll and Arthur Schweiger, "EasySpin, a comprehensive software
% package for spectral simulation and analysis in EPR," In: J. Magn.
% Reson. 178(1), 42-55 (2006)
% -----------------------------------------------------------------------
% The external routine is not part of the SOFiA MIT License.


function PnmRot = sofia_wdr(Pnm, xAngle, yAngle, zAngle, deg)

disp('SOFiA W/D/R Wigner-D Rotation R13-0306');

if nargin < 5
    deg = false;
end

if deg
    xAngle = xAngle*180/pi;
    yAngle = yAngle*180/pi;
    zAngle = zAngle*180/pi;
end

PnmRot = zeros(size(Pnm));

for i=0:sqrt(size(Pnm,1))-1
    wignerD = wignerd(i,[xAngle, yAngle, zAngle],'-');
    for bin = 1:size(Pnm,2)
        PnmRot(i^2+1:(i+1)^2,bin) = (Pnm(i^2+1:(i+1)^2, bin)' * wignerD')';
    end
end

