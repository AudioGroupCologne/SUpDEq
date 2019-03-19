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
% Make Matrix 3D Data (HD) R13-0306
%
%   Preparing pressure data for visual3d() or mtxToGixel()
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% mtxData = sofia_makeMtx(N, Pnm, dn, krIndex) 
% ------------------------------------------------------------------------     
% mtxData   SOFiA 3D-matrix-data in 1° steps
% ------------------------------------------------------------------------              
% N         Order of the spatial fourier transform     [default = 3]
% Pnm       Spatial Fourier Coefficients (from S/T/C)
% dn        Modal Radial Filters (from M/F)
% krindex   Index of kr Vector                         [default = 1]
% oversize  Integer Factor to increase the resolution. Set oversize = 1
%           (default) to use the mtxData matrix for visual3D(), map3D().
%
% Dependencies: SOFiA P/D/C
% 
% The file generates a SOFiA mtxData Matrix of 181x360 pixels for the
% visualisation with sofia_visual3d() in 1° Steps (65160 plane waves).
% The HD version generally allows to raise the resolution (oversize > 1).  
% (visual3D(), map3D() admit 1° data only, oversize = 1)
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

function mtxData = sofia_makeMTX(N, Pnm, dn, krIndex, oversize)  

disp('SOFiA Visual Matrix Generator Wrapper HD R13-0306');

if nargin < 5
   oversize = 1;
end

oversize = round(oversize);

if oversize < 1 
    oversize = 1;
end
   
angle = zeros(65160*oversize^2,2);
cnt=1;

Elev = linspace(0,180, 181*oversize);
Azim = linspace(0,359, 360*oversize);

for ElevCNT=1:size(Elev,2)
     for AzimCNT=1:size(Azim,2)         
         angle(cnt,:) = [Azim(AzimCNT) Elev(ElevCNT)];
         cnt = cnt+1;
     end
end

angle=angle*pi/180;
Y = sofia_pdc(N, angle, Pnm(:,krIndex), dn(:,krIndex));

cnt=1;
mtxData=zeros(181*oversize,360*oversize);

for ElevCNT=1:size(Elev,2)
    for AzimCNT=1:size(Azim,2)
        mtxData(ElevCNT, AzimCNT)=Y(cnt);
        cnt=cnt+1;
    end
end

