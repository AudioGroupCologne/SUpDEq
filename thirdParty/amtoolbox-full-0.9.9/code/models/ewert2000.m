function [fcs powers] = ewert2000(data,fig,fs)
%ewert2000  Modulation filterbank, EPSM version
%   Usage: [fcs,powers] = ewert2000(data,fig,fs)
%
%   Input parameters:
%    freqs  :  Frequencies of the data, starting from f = 0
%    data   :  Powerspectrum data to be filtered.
%    fcs    :  Centre frequencies of the filters
%    powers :  Integrated power at the output of each filter
%    Wcf    :  Matrix containing the frequency-domain squared transfer function
%    fig    :  flag for plotting fitler transferfunctions 1 = yes, 0 = no
%              plotting
%
%   [fcs powers] = EWERT2000(data,fig,fs) computes the
%   EPSM-filterbank as presented by Ewert & Dau 2000, without the additional
%   lowpass filter with a cut-off at 150 Hz. This implementation consists of
%   a lowpass filter with a cutoff at 1 Hz, in parallel with 6 bandpass
%   filters with octave spacing. the Center-frequencies of the bandpass
%   filters are lower than the original from Ewert & Dau (2000). In each of
%   the filters data is integrated within the pass-band of the filter to
%   give the power within that filter.
%
%   References:
%     S. Ewert and T. Dau. Characterizing frequency selectivity for envelope
%     fluctuations. J. Acoust. Soc. Am., 108(3):1181-1196, 2000.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/ewert2000.php

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

% AUTHOR: Søren Jørgensen

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
% fig = 1;
if fig == 1
    figure
    fnts = 14;
    lw = 2;
    % plot(freqs,10*log10(abs(TFs(1,:))),'linewidth',lw), hold on
    for k = 1:7
        plot(freqs,10*log10(TFs(k,:)),'linewidth',lw),hold on
    end
    title('Squared transfer functions fo the filterbank')
    xlabel('Frequency [Hz]','FontSize',fnts)
    ylabel('Filter attenuation [dB]','FontSize',fnts)
    % set(gca,'XScale','linear','Xtick',fcs,'FontSize',fnts,'FontWeight','b');
    xlim([0 79])
    ylim([-20 5])
    
    % signalpower = data(1);
    % 1/data(1)*sum(Vout(8,:))
    figure
    semilogx(pos_freqs(2:end),10*log10(Vout(1,:))), hold on
    for k = 2
        semilogx(pos_freqs(2:end),10*log10(Vout(k,:)),'r')
    end
    % xlim([0 32])
    % ylim([-90 10])
    
    figure
    % for k = 1:9
    plot(outTimeEPSM), hold on
    % end
    % plot(time_data,'r')
end

