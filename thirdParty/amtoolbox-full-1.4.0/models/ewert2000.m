function [fcs, powers, varargout] = ewert2000(data,fs)
%EWERT2000 Modulation filterbank (based on EPSM)
%   Usage: [fcs,powers] = ewert2000(data,fs)
%
%   Input parameters:
%     data   : Powerspectrum data to be filtered.
%     fs     : sampling frequency
%
%   Output parameters:
%     fcs   : centre frequencies of the filters
%
%     powers: filter coefficients
%
%   [fcs powers] = EWERT2000(data,fs) computes the
%   EPSM-filterbank as presented by Ewert & Dau 2000, without the additional
%   lowpass filter with a cut-off at 150 Hz. This implementation consists of
%   a lowpass filter with a cutoff at 1 Hz, in parallel with 6 bandpass
%   filters with octave spacing. the Center-frequencies of the bandpass
%   filters are lower than the original from Ewert & Dau (2000). In each of
%   the filters data is integrated within the pass-band of the filter to
%   give the power within that filter.
%
%   See also: demo_ewert2000
%
%   References:
%     S. Ewert and T. Dau. Characterizing frequency selectivity for envelope
%     fluctuations. J. Acoust. Soc. Am., 108(3):1181--1196, 2000.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/ewert2000.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Author: Søren Jørgensen
%   #Author: Clara Hollomey (2020): moved plot to demo

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if mod(length(data),2) == 0
    %number is even
    data = data(1:end-1);
else
    %number is odd
end


Q = 1;
N = length(data);
X = fft(data);
X_mag  = abs(X) ;
X_power = X_mag.^2/N ;% power spectrum.
X_power_pos = X_power(1:fix(N/2)+1) ;
%take positive frequencies only and mulitply by two-squared to get the same total energy
X_power_pos(2:end) = X_power_pos(2:end).* (2)  ; 
pos_freqs= linspace(0,fs/2,length(X_power_pos));
freqs = [pos_freqs -1*fliplr(pos_freqs(2:end))];

%band center frequencies
fcs=[ 2 4 8 16 32 64 ];


% Initialize transfer function
TFs = zeros(length(fcs),length(freqs));

% Calculating frequency-dmoain transferfunction for each center frequency:
for k = 1:length(fcs)
    TFs(k+1,2:end) = 1./(1+ (1j*Q*(freqs(2:end)./fcs(k) - fcs(k)./freqs(2:end)))); % p287 Hambley.
end


% squared filter magnitude transfer functions
Wcf = (abs(TFs)).^2;
% cutoff frequency of lowpassfilter:
fcut = 1;
% order:
n = 3;

% Lowpass filter squared transfer function: third order butterworth filter
% TF from: http://en.wikipedia.org/wiki/Butterworth_filter
Wcf(1,:) =  1./(1+((2*pi*freqs/(2*pi*fcut)).^(2*n))); 

TFs(1,:) = sqrt(Wcf(1,:));
% initialize output product:
Vout = zeros(length(fcs),length(pos_freqs));
powers = zeros(1,7);


% ------------ DC-power, --------------------------
% %here devided by two such that a fully modulated tone has an AC-power of 1.
DC_power =  X_power_pos(1)/ N /2;

% ------------------------------------------------

for k = 1:size(Wcf,1)
    Vout(k,:) = X_power_pos.*Wcf(k,1:floor(end/2)+1); 
    % Integration estimated as a sum from f > 0 
    % integrate envelope power in the
    % passband of the filter. Index goes from 2:end since
    powers(k) =  sum(Vout(k,2:end))/ N / DC_power;
    % integration is for f>0
end
fcs = [1 fcs];

if nargout >= 3
  varargout{1} = TFs;
  if nargout >= 4
    varargout{2} = freqs;
  end
end


