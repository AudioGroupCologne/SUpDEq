function output = takanen2013_mso(ipsilateral, contralateral, fs, fc, printfigs)
%takanen2013_mso Model of the medial superior olive
%   Usage: output = takanen2013_mso(ipsilateral, contralateral, fs, fc, printfigs)
%
%   Input parameters:
%        ipsilateral   : The ipsilateral "where" stream output from the
%                        model of the periphery
%        contralateral : The contralateral "where" stream output from the
%                        model of the periphery
%        fs            : sampling rate
%        printFigs     : boolean value that defines whether several figures
%                        illustrating the processing steps in the model are
%                        plotted or not. As default, no figures are 
%                        plotted.
%
%   Output parameters:
%        output : Spatial cues for separate narrow bandwidths
%
%   This function models the medial superior olive (MSO) by processing the 
%   output of the periphery model with the following steps:
%
%   1) Delay and convolution is applied to the contralateral signal and
%      the ipsilateral signal is limited and normalized.
%
%   2) Coincidence detection between the ipsilateral and contralateral
%      signals.
%
%   3) Weighted and self-weighted moving average filters are applied to
%      the outputs of the coincidence detection and contralateral signal,
%      respectively, and the output is limited.
%
%   See also: takanen2013, takanen2013_periphery, takanen2013_weightedaveragefilter
%
%   References:
%     V. Pulkki and T. Hirvonen. Functional count-comparison model for
%     binaural decoding. Acta Acustica united with Acustica, 95(5):883 - 900,
%     Sept./Oct. 2009.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147-163, 2014. PMID: 24513586.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/takanen2013_mso.php

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

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland


nrows = size(contralateral,1);
t=(0:(nrows-1))./fs;
% If desired, the computations are illustrated at two characteristic
% frequencies
band=[8,10];

%% ------ The contralateral ear input is delayed by 0.4 ms ----------------
contraDelay = round(0.0004*fs);
contralateral = [zeros(contraDelay,size(contralateral,2));...
    contralateral(1:size(contralateral,1)-contraDelay,:)];

if(printfigs)
    figure(91);
    g(1)=subplot(2,1,1);plot(t,ipsilateral(:,band(1)),'-b',t,contralateral(:,band(1)),'--r');
    g(2)=subplot(2,1,2);plot(t,ipsilateral(:,band(2)),'-b',t,contralateral(:,band(2)),'--r');
    linkaxes(g,'x');title('Ipsi and contra inputs');
end

%% ------ Convolution with the contra response ----------------------------
limit2 = find(fc>=1500,1,'first');
for freqind = 1:length(fc)
    if freqind<limit2
        n = (0:1:(fs/fc(freqind)))';
        f =0.25*(cos(2*pi*(fc(freqind)*n/fs).^.25-pi)+1).^2;
    end
    convolved = conv(contralateral(:,freqind),f)/sum(f);
    contralateral(:,freqind) = convolved(1:nrows);
end

if(printfigs)
    figure(92);
    g(1)=subplot(2,1,1);plot(t,ipsilateral(:,band(1)),'-b',t,contralateral(:,band(1)),'--r');
    g(2)=subplot(2,1,2);plot(t,ipsilateral(:,band(2)),'-b',t,contralateral(:,band(2)),'--r');
    linkaxes(g,'x');title('Ipsi and contra after the contra response');
end

%% ------ Limiting of the ipsilateral input -------------------------------
% Limit values are determined for each frequency based on white noise at 30
% dB SPL
x=amt_load('takanen2013','msolimits.mat');
limits=x.limits;
limited = ipsilateral./(ones(nrows,1)*limits);
limited(limited>1) =1;

%% ------ Coincidence detection -------------------------------------------
output = (limited.*contralateral);

if(printfigs)
    figure(93);
    plot(t,output(:,band(2)),t,contralateral(:,band(2)),'-r');
    title('Coincidence and contralateral');
end

%% ------ Weighted and self-weighted moving averages of 4 ms --------------
output = takanen2013_weightedaveragefilter(output,contralateral,fs,0.004) ./ (takanen2013_weightedaveragefilter(contralateral,contralateral,fs,0.004)+1e-30);
% tau = 0.004; 
% B = 1-exp(-1/(tau*fs));
% A = [1 -exp(-1/(tau*fs))];
% self_weighted = filter(B,A,(contralateral.^3));
% weighted = filter(B,A,(output.*(contralateral.^2)));
% output = weighted./(self_weighted+1e-30);
output(output>1) = 1;
