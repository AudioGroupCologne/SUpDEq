function [midfreq_out rms_out] = thirdoctrmsanalysis24(insig,fs)
%THIRDOCTRMSANALYSIS  XXX Description
%   Usage: [midfreq_out rms_out] = thirdoctrmsanalysis24(x,fs);
%
%   Input parameters:
%     insig  : Input signal
%     fs     : Sampling frequency
%
%   Output parameters:
%     midfreq_out : XXX DESC
%     rms_out     : XXX DESC
%
%   [midfreq_out,rms_out] = THIRDOCTRMSANALYSIS24(x,fs) computes XXX
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/thirdoctrmsanalysis24.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% insig = input;
% insig = real(ifft(ones(1,1000)));
% fs = 2000;
N = length(insig);
X = (fft(insig));
X_mag  = abs(X) ;
X_power = X_mag.^2/N ;% power spectrum.
X_power_pos = X_power(1:fix(N/2)+1) ;

%take positive frequencies only and mulitply by two-squared to get the same
%total energy(used since the integration is only performed for positive
%freqiencies)
X_power_pos(2:end) = X_power_pos(2:end).* (2);

freq= linspace(0,fs/2,length(X_power_pos));

% freq = linspace(0,fs/2,N/2 +1);
% X_pos = X(1:N/2);
f_axis = [-1*freq(end:-1:2)  freq(1:end-1)];

%resolution of data
resol=freq(2)-freq(1);

%band cross-over frequencies
%  not needed: 16 20 25 32 40 50 63 80 100 % % 6300 8000 10000 12500 16000
%  20000
midfreq=[63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 ...
         2500 3150 4000 5000 6300 8000 10000 12500 ];
crossfreq(1)=midfreq(1)/(2^(1/6));
crossfreq(2:length(midfreq)+1)=midfreq*(2^(1/6));

%cross-over indicies
y=crossfreq/resol;

%initialize output matrix
% output = zeros(length(midfreq));
% output_specs = zeros(length(midfreq));
%rounding up
crosselem=round(y);
for n=1:length(y)
    if crosselem(n)<y(n)
        crosselem(n)=crosselem(n)+1;
    end
end

nn=1;
while crossfreq(nn+1)<=freq(end)
    
    rms_out(nn) = sqrt(sum(X_power_pos(crosselem(nn):crosselem(nn+1)-1))/N);
   
    nn=nn+1;  
end
midfreq_out = midfreq(1:nn-1);

