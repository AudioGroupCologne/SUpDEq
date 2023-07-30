function [output energy] = takanen2013_wbmso(ipsilateral, contralateral, fs, widthinerbs, fc, printfigs)
%TAKANEN2013_WBMSO Wideband MSO model
%   Usage: [output energy] = takanen2013_wbmso(ipsilateral, contralateral, fs, widthinerbs, fc, printfigs)
%
%   Input parameters:
%        ipsilateral   : The ipsilateral "where" stream output from the
%                        model of the periphery
%        contralateral : The contralateral "where" stream output from the
%                        model of the periphery
%        fs            : Sampling rate
%        widthinerbs   : The number of adjacent ERB bands the information 
%                        is gathered over
%        fc            : Characteristic frequencies
%        printFigs     : Boolean value that defines whether several figures
%                        illustrating the processing steps in the model are
%                        plotted or not. As default, no figures are 
%                        plotted.
%
%   Output parameters:
%        output : Spatial cues for frequency bands that are summed together
%                 to form a width defined by widthinerbs
%        energy : "What" stream for the wideband MSO
%
%   The wideband MSO model simulates the ability of the human auditory 
%   system to extract localization cues based on interaural envelope time 
%   shifts. These time shifts are decoded into directional cues for the 
%   model. This is done by processing the output of the periphery model 
%   with the following steps:
%
%   1) Delay the contralateral signal.
%
%   2) Adjacent frequency bands are summed on both sides.
%
%   3) The average of the signal is removed to emphasize prominent peaks
%      by applying a self-weighted moving average filter and delay to the
%      signal and deducting this from the signal after summing over
%      frequency bands.
%
%   4) Both sides are convolved with a Hanning window and a phase-locked
%      impulse generator is applied.
%
%   5) The contralateral side is convolved and the ipsilateral side is
%      limited and normalized.
%
%   6) Coincidence detection between the ipsilateral and contralateral
%      signals.
%
%   7) Weighted and self-weighted moving average filters are applied to
%      the outputs of the coincidence detection and contralateral signal,
%      respectively, and the output is limited.
%
%   See also: takanen2013, takanen2013_periphery, takanen2013_weightedaveragefilter
%
%   References:
%     V. Pulkki and T. Hirvonen. Functional count-comparison model for
%     binaural decoding. Acta Acustica united with Acustica, 95(5):883 --
%     900, Sept./Oct. 2009.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147--163, 2014. PMID: 24513586.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/takanen2013_wbmso.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB M-Signal
%   #Author: Marko Takanen (2013)
%   #Author: Olli Santala 
%   #Author: Ville Pulkki

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%if desired, the computations are illustrated at two characteristic
%frequencies, namely around 500 Hz and 4.5 kHz
band=[8,25];
t=(0:(size(ipsilateral,1)-1))./fs;
[nrows nBands] = size(contralateral);

%% ------ The contralateral ear input is delayed by 0.2 ms ----------------
contradelay = round(0.0002*fs);
contralateral = [zeros(contradelay,nBands);...
    contralateral(1:size(contralateral,1)-contradelay,:)];

%% ------ Calculation of the sums across frequency bands ------------------
ipsi = ipsilateral;
contra = contralateral;
for freqind = 1:length(fc)
    columnRange = max(1,freqind-floor(widthinerbs/2)):min(length(fc),freqind+floor(widthinerbs/2));
    contralateral(:,freqind) = sum(contra(:,columnRange),2);
    ipsilateral(:,freqind) = sum(ipsi(:,columnRange),2);
end
tempn =ipsilateral;
tempm = contralateral;

%% ------ Post-processing of periphery output -----------------------------
mcoeff=2;
tempL = mcoeff*takanen2013_weightedaveragefilter(ipsilateral,ipsilateral,fs,0.01);
tempR = mcoeff*takanen2013_weightedaveragefilter(contralateral,contralateral,fs,0.01);

ipsilateral = ipsilateral-[zeros(floor(0.0005*fs),nBands);tempL(1:end-floor(0.0005*fs),:)];
contralateral = contralateral-[zeros(floor(0.0005*fs),nBands);tempR(1:end-floor(0.0005*fs),:)];

ipsilateral = ipsilateral.*(ipsilateral>0);
contralateral=contralateral.*(contralateral>0);

if(printfigs)
    figure(94);
    g(1)=subplot(2,1,1);plot(t,tempn(:,band(2)),'-b',t,tempm(:,band(2)),'--r');
    g(2)=subplot(2,1,2);plot(t,ipsilateral(:,band(2)),'-b',t,contralateral(:,band(2)),'--r');
    linkaxes(g,'x');
    title('Ipsi and contra after the short-time average');
end

%% ------ Convolution with a Hanning window -------------------------------
x = hanning(floor(0.002*fs));
maxWidth =11; % the length of the square-wave is at maximum of 0.2-ms long
for i=1:length(fc)
    %convolution with a gaussian window
    temp = conv(ipsilateral(:,i),x','same')./sum(x);
    %the original values are stored for plotting purposes into temporary variable 
    tempn(:,i) = temp;
    ipsilateral(:,i) = zeros(nrows,1);

    startingpoints = strfind((temp'>0),[0 1]);
    endpoints = [strfind((temp'>0),[1 0])+1 nrows];
    if(isempty(startingpoints))
        startingpoints = 1;
    end
    if(startingpoints(1)>=endpoints(1))
        startingpoints=[1 startingpoints];
    end
    %search for the mass centroid of a half-wave
    for indx=1:length(startingpoints)
        tempsum = sum(temp(startingpoints(indx):endpoints(indx)));
        temp2sum = cumsum(temp(startingpoints(indx):endpoints(indx)));
        location = startingpoints(indx)+find((temp2sum>=tempsum/2),1,'first');
        peak = max(temp(startingpoints(indx):endpoints(indx)));
        range = max(1,(location-floor(maxWidth/2))):min(nrows,(location+floor(maxWidth/2)));
        ipsilateral(range,i) = ones(length(range),1)*peak;
    end
    %convolution with a gaussian window
    temp = conv(contralateral(:,i),x','same')./sum(x);
    %the original values are stored for plotting purposes into temporary variable
    tempm(:,i) = temp;
    contralateral(:,i)= zeros(nrows,1);
    startingpoints = strfind((temp'>0),[0 1]);

    endpoints = [strfind((temp'>0),[1 0])+1 nrows];
    if(isempty(startingpoints))
        startingpoints = 1;
    end
    if(startingpoints(1)>=endpoints(1))
        startingpoints=[1 startingpoints];
    end
    %search for the mass centroid of a half-wave
    for indx=1:length(startingpoints)
        tempsum = sum(temp(startingpoints(indx):endpoints(indx)));
        temp2sum = cumsum(temp(startingpoints(indx):endpoints(indx)));
        location = startingpoints(indx)+find((temp2sum>=tempsum/2),1,'first');
        peak = max(temp(startingpoints(indx):endpoints(indx)));
        range = max(1,(location-floor(maxWidth/2))):min(nrows,(location+floor(maxWidth/2)));
        contralateral(range,i) = ones(length(range),1)*peak;
    end
end

if(printfigs)
    figure(95);
    g(1)=subplot(2,1,1);plot(t,tempn(:,band(2)),'-b',t,tempm(:,band(2)),'--r');
    g(2)=subplot(2,1,2);plot(t,ipsilateral(:,band(2)),'-b',t,contralateral(:,band(2)),'--r');
    linkaxes(g,'x');
    title('Ipsi and contra after the short-time average');
end

%% ------ Computing of the energy output of the wideband MSO  -------------
x=amt_load('takanen2013','wbmsomultp.mat');
multp=x.multp;
energy = contralateral.*(ones(nrows,1)*multp);
temp = takanen2013_weightedaveragefilter(energy,energy,fs,0.01);
energy = temp.*(ones(nrows,1)*(max(energy)./max(temp)));

%% ------ Convolution with the contra response ----------------------------
n = (0:1:(fs/fc(1)))';
f =0.25*(cos(2*pi*(fc(1)*n/fs).^.25-pi)+1).^3;
for freqind = 1:length(fc)
    convolved = conv(contralateral(:,freqind),f)./sum(f);
    contralateral(:,freqind) = convolved(1:nrows);
end

%% ------ Limiting of the ipsilateral input -------------------------------
limits = 1e-10;
limited=ipsilateral./limits;%
limited(limited>1) =1;

if(printfigs)
    figure(96);
    g(1)=subplot(2,1,1);plot(t,ipsilateral(:,band(1)),'-b',t,contralateral(:,band(1)),'--r');
    g(2)=subplot(2,1,2);plot(t,ipsilateral(:,band(2)),'-b',t,contralateral(:,band(2)),'--r');
    linkaxes(g,'x');
    title('Ipsi and contra after the contra response.');
end

%% ------ Coincidence detection -------------------------------------------
output = (limited.*contralateral);

if(printfigs)
    figure(97);
    plot(t,output(:,band(2)),t,contralateral(:,band(2)),'-r');
    title('Coincidence and contralateral');
end

%% ------ Weighted and self-weighted moving averages of 1 ms --------------
output = takanen2013_weightedaveragefilter(output,contralateral,fs,0.001) ./ (takanen2013_weightedaveragefilter(contralateral,contralateral,fs,0.001)+1e-30);
% tau = 0.001; 
% B = 1-exp(-1/(tau*fs));A = [1 -exp(-1/(tau*fs))];
% selfweighted = filter(B,A,(contralateral.^3));
% weighted = filter(B,A,(output.*(contralateral.^2)));
% output = (weighted)./(selfweighted+1e-80);
output(output>1) = 1;


