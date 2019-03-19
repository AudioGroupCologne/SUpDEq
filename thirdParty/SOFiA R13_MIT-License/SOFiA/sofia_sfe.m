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
% S/F/E Sound Field Extrapolation R13-0306
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% Pnm_krb = sofia_sfe(Pnm_kra, kra, krb, problem) 
% ------------------------------------------------------------------------     
% Pnm_krb Spatial Fourier Coefficients, extrapolated to rb
% ------------------------------------------------------------------------              
% 
% Pnm_kra Spatial Fourier Coefficients from SOFiA S/T/C
% 
% kra     k*ra Vector
% krb     k*rb Vector
% problem 1: Interior problem [default]  2: Exterior problem
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


function Pnm_krb = sofia_sfe(Pnm_kra, kra, krb, problem) 

disp('SOFiA S/F/E - Sound Field Extrapolation R13-0306');
disp(' ');

if nargin < 3
    error('Arguments missing. 3 Args required: Pnm_krb = sofia_sfe(Pnm_kra, kra, krb)')
end

if nargin < 4
   problem = 1;
end

if size(kra,2) ~= size(krb,2) ||  size(kra,2) ~= size(Pnm_kra,2) 
   error('FFT bin number or dimension mismatch. Pnm_kra, kra and krb must have the same M-dimension.') 
end

FCoeff  = size(Pnm_kra,1);
N       = sqrt(FCoeff)-1;

nvector=zeros(FCoeff,1);

index = 1;
for n=0:N
    for m=-n:n
        nvector(index) = n;
        index = index+1;
    end
end

nvector = repmat(nvector,1,size(Pnm_kra,2));
kra     = repmat(kra,FCoeff,1);
krb     = repmat(krb,FCoeff,1);

if problem == 2
    hn_kra  = sqrt(pi./(2*kra)).*besselh(nvector+.5,1,kra);
    hn_krb  = sqrt(pi./(2*krb)).*besselh(nvector+.5,1,krb);
    exp     = hn_krb./hn_kra;
else
    jn_kra  = sqrt(pi./(2*kra)).*besselj(n+.5,kra);
    jn_krb  = sqrt(pi./(2*krb)).*besselj(n+.5,krb);
    plot(abs(jn_kra'))
    exp     = jn_krb./jn_kra;
    if ~isempty(find(abs(exp)>1e2)) %40dB
       disp('WARNING: Extrapolation might be unstable for one or more frequencies/orders!');
    end
end

Pnm_krb = Pnm_kra.*exp;

