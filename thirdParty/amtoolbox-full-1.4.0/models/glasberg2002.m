function [results] = glasberg2002(inSig,fs, varargin)
%GLASBERG2002  Loudness model for time-variant signals
%   Usage: [results] = glasberg2002(inSig,fs);
%   
%   Input parameters:
%     inSig : monaural input signal
%     fs : sampling frequency [Hz]
%
%   Output parameters:
%     results : output struct containing the long-term and short-term loudness
%
%   Optional parameters:
%
%     'fs',fs     model-internal sampling frequency [Hz]; it is strongly
%                 recommended to use the default of 32000 Hz
%
%     'flow',flow   lowest frequency at which to evaluate the outer/middle
%                   ear transfer function
%
%     'fhigh',fhigh   highest frequency at which to evaluate the outer/middle
%                     ear transfer function
%
%     'order',order   order of the FIR filter used for deriving the outer/middle
%                     ear transfer function
%
%     'fftLen',fftLen     length of the fft used for preprocessing the inputSignal
%
%     'vLimitingIndices',vLimitingIndices     vector [1 x 5] containing boundaries [Hz] at which to
%                                             separate the inputSignal during frequency analysis
%
%     'hannLenMs',hannLenMs     vector [1 x 5] containing the size [ms] of the
%                               windows used for frequency analysis of the inputSignal
%
%     'timeStep',timeStep     integration duration for short-term loudness [s]
%
%     'compensationtype'     flag; specifies the outer/middle ear transfer
%                            function. Can be 'tfOuterMiddle1997','tfOuterMiddle2007', or 'specLoud'
%
%   This code calculates the excitation patterns, the long-term and
%   short-term loudness. 
%
%   Example:
%
%     fs = 32000; 
%     t = linspace(0,1,fs);
%     sig = sin(2*pi*1000*t).';
%     inSig = scaletodbspl(sig,100);  
%
%   Note that currently fs must be 32000 Hz.
%
%   See also: data_glasberg2002 exp_glasberg2002 exp_moore1997
%             glasberg2002 moore1997 moore2016
%
%   References:
%     B. R. Glasberg and B. C. J. Moore. A Model of Loudness Applicable to
%     Time-Varying Sounds. J. Audio Eng. Soc, 50(5):331--342, 2002.
%     
%     B. C. J. Moore, B. R. Glasberg, and T. Baer. A Model for the Prediction
%     of Thresholds, Loudness, and Partial Loudness. J. Audio Eng. Soc,
%     45(4):224--240, 1997.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/glasberg2002.php


%   #StatusDoc: Satisfactory
%   #StatusCode: Submitted
%   #Verification: Untrusted
%   #Requirements: M-Signal
%   #Author: Thomas Deppisch (2017)
%   #Author: Piotr Majdak (2017)
%   #Author: Clara Hollomey (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% todo: two-channel inSig, different fs (-> window sizes), inSig normalization

%% model
definput.import = {'glasberg2002'};    

[flags,kv]  = ltfatarghelper({},definput,varargin);

if size(kv.hannLenMs, 2) ~= 6 || size(kv.vLimitingIndices, 2) ~= 5
    error('vLimitingIndices needs to be of dim [1 x 5] and hannLenMs of dim [1 x 6]');
end

if fs ~= kv.fs
    inSig = resample(inSig, kv.fs, fs);
    fs = kv.fs;
end

%     compute ffts in 1ms intervals
fVec = kv.flow:kv.fhigh;
data = data_glasberg2002(flags.compensationtype,'fVec',fVec);

% create FIR filter
tfLinear = 10.^(data.tfOuterMiddle/10);
outerMiddleFilter = fir2(kv.order, linspace(0, 1, length(fVec)), tfLinear);
earSig = filtfilt(outerMiddleFilter,1,inSig);   % why does filter(..) not work?

updateRate = round(kv.timeStep*kv.fs);
numBlocks = ceil(length(earSig)./updateRate);
hannLenSmp = round(kv.hannLenMs./1000 * kv.fs); % windows size in samples
earSigPad = earSig;
earSigPad(end+1:end+hannLenSmp(6)) = zeros(hannLenSmp(6),1);  % zero padding

% calculate hann windows
hannWin2 = hann(hannLenSmp(1));
hannWin4 = hann(hannLenSmp(2));
hannWin8 = hann(hannLenSmp(3));
hannWin16 = hann(hannLenSmp(4));
hannWin32 = hann(hannLenSmp(5));
hannWin64 = hann(hannLenSmp(6));

spect1 = zeros(numBlocks, kv.fftLen);
spect2 = zeros(numBlocks, kv.fftLen);
spect3 = zeros(numBlocks, kv.fftLen);
spect4 = zeros(numBlocks, kv.fftLen);
spect5 = zeros(numBlocks, kv.fftLen);
spect6 = zeros(numBlocks, kv.fftLen);

for ii=1:numBlocks
    % fft for each window
    lower = (ii-1)*updateRate+1;
    upper = (ii-1)*updateRate+hannLenSmp(1);
    spect1(ii,:) = fft(earSigPad(lower:upper) .* hannWin2, kv.fftLen);

    upper = (ii-1)*updateRate+hannLenSmp(2);
    spect2(ii,:) = fft(earSigPad(lower:upper) .* hannWin4, kv.fftLen);

    upper = (ii-1)*updateRate+hannLenSmp(3);
    spect3(ii,:) = fft(earSigPad(lower:upper) .* hannWin8, kv.fftLen);

    upper = (ii-1)*updateRate+hannLenSmp(4);
    spect4(ii,:) = fft(earSigPad(lower:upper) .* hannWin16, kv.fftLen);

    upper = (ii-1)*updateRate+hannLenSmp(5);
    spect5(ii,:) = fft(earSigPad(lower:upper) .* hannWin32, kv.fftLen);

    upper = (ii-1)*updateRate+hannLenSmp(6);
    spect6(ii,:) = fft(earSigPad(lower:upper) .* hannWin64, kv.fftLen);
end

    oneHz = (kv.fftLen+2)/kv.fs;  % number of frequency bins representing 1Hz    
    % truncate ffts to match the frequency ranges specified in glasberg2002
    % and put PSD for each window and each time interval in matrix -> window
    % normalization
    spect = zeros(numBlocks,kv.fftLen/2+1);
    spect(:,round(kv.vLimitingIndices(1)*oneHz)+1:kv.fftLen/2 +1) = abs(spect1(:,round(kv.vLimitingIndices(1)*oneHz)+1:kv.fftLen/2+1)).^2/sum(hannWin2.^2); % 4050-fs/2
    spect(:,round(kv.vLimitingIndices(2)*oneHz)+1:round(kv.vLimitingIndices(1)*oneHz)) = abs(spect2(:,round(kv.vLimitingIndices(2)*oneHz)+1:round(kv.vLimitingIndices(1)*oneHz))).^2/sum(hannWin4.^2); % 2540-4050Hz
    spect(:,round(kv.vLimitingIndices(3)*oneHz)+1:round(kv.vLimitingIndices(2)*oneHz)) = abs(spect3(:,round(kv.vLimitingIndices(3)*oneHz)+1:round(kv.vLimitingIndices(2)*oneHz))).^2/sum(hannWin8.^2); % 1250-2540Hz
    spect(:,round(kv.vLimitingIndices(4)*oneHz)+1:round(kv.vLimitingIndices(3)*oneHz)) = abs(spect4(:,round(kv.vLimitingIndices(4)*oneHz)+1:round(kv.vLimitingIndices(3)*oneHz))).^2/sum(hannWin16.^2); % 500-1250Hz
    spect(:,round(kv.vLimitingIndices(5)*oneHz)+1:round(kv.vLimitingIndices(4)*oneHz)) = abs(spect5(:,round(kv.vLimitingIndices(5)*oneHz)+1:round(kv.vLimitingIndices(4)*oneHz))).^2/sum(hannWin32.^2); % 80-500Hz
    spect(:,1:round(kv.vLimitingIndices(5)*oneHz)) = abs(spect6(:,1:round(kv.vLimitingIndices(5)*oneHz))).^2/sum(hannWin64.^2); % 0-80Hz
    compInt = 2*spect./fs;  %   psd
    compFq = linspace(0,fs/2,kv.fftLen/2+1);
    
    %% calculating excitation patterns
    % calculate ERB numbers corresponding to ERB mid frequencies
    erbStep = 0.25;    % according to moore1997, glasberg2002
    erbFcMin = 50;
    erbFcMax = 15000;
    erbNMin = fc2erb(erbFcMin);
    erbNMax = fc2erb(erbFcMax);
    erbN = erbNMin:erbStep:erbNMax;    % numbers of erb bands
    erbFc = erb2fc(erbN);               % center frequency of erb bands

    erbLoFreq = erb2fc(erbN-0.5); % lower limit of each ERB filter
    erbHiFreq = erb2fc(erbN+0.5); % upper limit of each ERB filter

    % calculate intensity for each ERB (dB/ERB) and each time step
    erbInt = zeros( size(compInt, 1), length(erbFc));
    for ii=1:length(erbFc)
        range = round(erbLoFreq(ii)*oneHz):round(erbHiFreq(ii)*oneHz);
        erbInt(:,ii) = sum(compInt(:,range),2);   % intensity sum in each erb
    end
    erbdB = 10*log10(erbInt./(20e-6)^2);   % intensity level in each erb using reference SPL of 20 uPa

    % p determines bandwidth and slope of the erb filter and is generally
    % asymmetrical
    % p is roughly symmetrical for an excitation level of 51dB per ERB
    p51 = 4*erbFc./f2erb(erbFc);    % p for erb center frequencies and a level of 51dB
    p511 = 4*1000/f2erb(1000);    % p for fc=1kHz and a level of 51dB (at 1kHz filters are symmetrical)

    pU = p51;   % pU for all erbFc and all time steps
    g = abs(repmat(compFq,149,1) - repmat(erbFc.',1,length(compFq)))./repmat(erbFc.',1,length(compFq));    % normalized deviation of each f to erbFc for each erb band

    pL=zeros(numBlocks,length(compFq),length(erbFc));
    p=zeros(size(pL));
    w=p; e=p; 
    eL=zeros(numBlocks,length(erbFc));
    erbdB2f=zeros(numBlocks,length(compFq));
    for ii = 1:numBlocks % time steps
        erbdB2f(ii,:) = interp1([0 erbFc fs/2], [min(erbdB(ii,:)) erbdB(ii,:) min(erbdB(ii,:))], compFq);   % map erbFc to compFq
        pL(ii,:,:) = repmat(p51,length(compFq),1) - 0.35.*(repmat(p51,length(compFq),1)./p511).*(repmat(erbdB2f(ii,:).',1,length(erbN))-51); 
        p(ii,:,:) = pL(ii,:,:);
        for jj = 1:length(erbN)
            p(ii,round(erbFc(jj)*oneHz)+1:end,jj) = pU(jj); % p(f>erbFc) = pU
        end
        w(ii,:,:) = (1+squeeze(p(ii,:,:)).*g.').*exp(-squeeze(p(ii,:,:)).*g.');    % calculate weighting function
        e(ii,:,:) = squeeze(w(ii,:,:)).*repmat(compInt(ii,:).',1,length(erbN));
        eL(ii,:) = sum(squeeze(e(ii,:,:)),1);   % sum excitation level in each erb
    end
    results.eLdB = 10*log10(eL./(20e-6)^2);
    results.erbN = erbN;

    %% calculating specific loudness 

    dataSL = data_glasberg2002('specLoud','fVec',erbFc);
    tQdB = dataSL.tQ;
    tQ = 10.^(tQdB./10);
    tQdB500 = dataSL.tQ500;
    %gdB = dataSL.g;    % low level gain in cochlea amplifier
    g = 10.^((tQdB500-tQdB)/10);
    a = dataSL.a;    % parameter for linearization around absolute threshold
    alpha = dataSL.alpha;    % compressive exponent
    c = dataSL.c; % constant to get loudness scale to sone

    specLoud = zeros(size(eL));

    specLoud1 = c*(2*eL./(eL+repmat(tQ,numBlocks,1))).^1.5 .*((repmat(g,numBlocks,1).* ...
        eL+repmat(a,numBlocks,1)).^repmat(alpha,numBlocks,1)-repmat(a.^alpha,numBlocks,1));
    specLoud2 = c * ((repmat(g,numBlocks,1) .*eL+repmat(a,numBlocks,1)).^repmat(alpha,numBlocks,1) - ...
        repmat(a.^alpha,numBlocks,1));
    specLoud3 = c*(eL./1.04e6).^0.5;
    specLoud(eL<repmat(tQ,numBlocks,1)) = specLoud1(eL<repmat(tQ,numBlocks,1));
    specLoud(eL<=10^10 & eL>repmat(tQ,numBlocks,1)) = specLoud2(eL<=10^10 & eL>repmat(tQ,numBlocks,1));
    specLoud(eL>10^10) = specLoud3(eL>10^10);

    %% monaural/binaural loudness (= instantaneous loudness), short term loudness (STL), long term loudness (LTL)
    results.monauralLoudness = sum(specLoud,2) * erbStep;     % integrate over the erbs
    results.binauralLoudness = 2*results.monauralLoudness;  % use moore2016 (Modeling binaural loudness) for better results

    % STL and LTL:
    aSTL = 0.045;
    rSTL = 0.02;
    STL = zeros(length(results.binauralLoudness),1);
    aLTL = 0.01;
    rLTL = 0.0005;
    LTL = zeros(length(results.binauralLoudness),1);
    for ii = 2:length(results.binauralLoudness)
        if results.binauralLoudness(ii)>STL(ii-1)
            STL(ii) = aSTL*results.binauralLoudness(ii)+(1-aSTL)*STL(ii-1);
        else
            STL(ii) = rSTL*results.binauralLoudness(ii)+(1-rSTL)*STL(ii-1);
        end
        if STL(ii)>LTL(ii-1)
            LTL(ii) = aLTL*STL(ii)+(1-aLTL)*LTL(ii-1);
        else
            LTL(ii) = rLTL*STL(ii)+(1-rLTL)*LTL(ii-1);
        end
    end
    results.STL = STL;
    results.LTL = LTL;


