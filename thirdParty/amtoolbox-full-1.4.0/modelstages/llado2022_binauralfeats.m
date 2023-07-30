function [x_feat] = llado2022_binauralfeats(ir_input,stim,fs)
%LLADO2022_BINAURALFEATS calculates binaural features
%   Usage: [x_feat] = llado2022_binauralfeats(ir_input,stim,fs);
%
%   Input parameters:
%     ir_input       : Impulse response according to the following matrix
%                      dimensions: direction x time x channel/ear
%     stim           : stimulus
%     fs             : samplef frequency
%
%   Output parameters:
%     x_feat         : binaural features according to the following matrix
%                      dimensions: nOutputs x features. Features combine the
%                      itd values below 1.5 kHz (rows 1-18) and ild values
%                      over 1 kHz (rows 19-36).
%
%   Binaural auditory model based on the example 13.6.2 in: Ville Pulkki and 
%   Matti Karjalainen. Communication  acoustics:  an  introduction  to  
%   speech,  audio and psychoacoustics.  John Wiley & Sons, 2015.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/llado2022_binauralfeats.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB M - Communication Systems
%   #Author: Pedro Llado (2022)
%   #Author: Petteri Hyvärinen (2022)
%   #Author: Ville Pulkki (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



if(ndims(ir_input) == 3)
    nIRs = length(ir_input(:,1,1));
    ir = ir_input;
else
    if(ismatrix(ir_input) == 1)
        nIRs = 1;
        ir(1,:,:) = ir_input;
    else
    amt_disp('Wrong input');
    end
end

%% USING AMTOOLBOX FOR GAMMA TONE MODEL %%%%%
%create of a gammatone filter bank using the AMT function
%(http://amtoolbox.sourceforge.net)
%
for n_ir = 1:nIRs
    
    %convolution: stimulus - hrirs
    insig = [conv(stim,ir(n_ir,:,1)) conv(stim,ir(n_ir,:,2))];
    %some parameters for auditory model
    fLow = 50;%the lowest characteristic frequency of the filter bank
    fHigh = 8000;%and the highest
    fCut = 1000; % cut-off frequency of the low-pass filter
    maxLag= floor(0.00067*fs); % IACC-values are computed within-1...1 ms

    %FILTERS CHARACTERISTIC FREQUENCIES
    cfs = erbspacebw(fLow,fHigh);   %charact. frequencies of the filter bank

    [b,a] = gammatone(cfs,fs,'complex');
    %'complex'  Generate filter modulated by exponential functions
    % b,a: nominator,denominator coefficients

    %% Windowing init
    % Determine length of input signal
    nSamples = size(insig,1);

    % Block parameter 
    blockSec  = 20e-3; % Block size in seconds
    hopSec    = 10e-3; % Step size in seconds


    % Block processing parameters
    blockSamples = 2 * round(fs * blockSec / 2);
    hopSamples   = 2 * round(fs * hopSec / 2);
    overSamples  = blockSamples - hopSamples;

    %Haning Window
    w = hann(blockSamples);

    % Calculate number of frames
    nFrames = fix((nSamples-overSamples)/hopSamples);


    for i = 1 : nFrames
        %GENERATING FILTERS (UFILTERBANKZ) AND FILTERING THE INDIVIDUAL SIGNALS
        %processing the signal through the filter bank
        filterOut(i).left=2*real(ufilterbankz(b,a,insig((i-1)*hopSamples+1:...
            (i-1)*hopSamples+blockSamples,1)));
        filterOut(i).right=2*real(ufilterbankz(b,a,insig((i-1)*hopSamples+1:...
            (i-1)*hopSamples+blockSamples,2)));

        %%%%%EMULATION OF THE NEURAL TRANSDUCTION with :
        % 1) half-wave rectification and
        rectified(i).left = filterOut(i).left.*(filterOut(i).left>0);
        rectified(i).right = filterOut(i).right.*(filterOut(i).right>0);
        % 2) low-pass filtering of the filter bank output
        %a first-order IIR filter is used as the low-pass filter
        beta = exp(-fCut/fs);
        outSig(i).left = filter(1-beta,[1 -beta],rectified(i).left);
        outSig(i).right = filter(1-beta,[1 -beta],rectified(i).right);

        %%%%%COMPUTE INTERAURAL CROSS-CORRELATION AT EACH FREQ BAND
        iaccFuncts = zeros(2*maxLag+1,length(cfs));
        lagValues = (-maxLag:maxLag)./fs;% vector of lags at fs from min to max
        for freqInd=1:length(cfs)
            iaccFuncts(:,freqInd) = xcorr(outSig(i).left(:,freqInd),...
                outSig(i).right(:,freqInd),maxLag,'coeff'); %COEFF: 
                                        %Normalizes the sequence so that the
                                        %autocorrelations at zero lag equal 1
        end

        %%%%%COMPUTE ITD FOR EACH FREQ BAND
        [~,lag] = max(iaccFuncts);
        itdEst(i,:) = lagValues(lag); % Time corresponding to lag displacement

        %%%%%COMPUTE ILD FOR EACH FREQ BAND
        ildEst(i,:) = dbspl(outSig(i).right(:,:))-dbspl(outSig(i).left(:,:));
        
        %% Store binaural faetures for the next stage of the model 
            

    end
    itdild_feat(n_ir,4:35) = mean(itdEst(:,1:end));
    itdild_feat(n_ir,36:67) = mean(ildEst(:,1:end));
    %%%%%
end
x_feat = [itdild_feat(:,4:21)';itdild_feat(:,50:end)']; % Use ITD below 1.5 kHz and ILD over 1 kHz
end


