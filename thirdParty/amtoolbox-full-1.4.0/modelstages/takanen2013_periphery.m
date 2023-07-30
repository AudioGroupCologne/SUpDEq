function periphOutput = takanen2013_periphery(insig,fs,outputPlot)
%TAKANEN2013_PERIPHERY Process the input through the model of periphery         
%   Usage: periphOutput = takanen2013_periphery(insig,fs,outputPlot);
%          periphOutput = takanen2013_periphery(insig,fs);
%
%   Input parameters:
%        insig      : binaural input signal to be processed. Optionally,
%                     the output of the cochlear model by Verhulst et. al.
%                     2012 can be used as well
%        fs         : sampling rate
%        outputPlot : boolean value that defines whether the periphery
%                     model output at different frequency bands are plotted
%
%   Output parameters:
%        periphOutput : Structure consisting of the following elements
%
%                       periphOutput.left          Left ear "where" stream output
% 
%                       periphOutput.right         Right ear "where" stream output
%
%                       periphOutput.fc            Characteristic frequencies
%
%                       periphOutput.ventralLeft   Left hemisphere "what" stream output 
%
%                       periphOutput.ventralLeft   Right hemisphere "what" stream output
%
%   This function processes the binaural input signal through the model of
%   periphery presented by Takanen, Santala, Pulkki 2013, the model which
%   consists of a nonlinear time-domain model of cochlea and model of
%   cochlear nucleus. The processing contains the following steps:
%
%   1) the binaural input signal is processed with the nonlinear time-
%      domain model of cochlea by Verhulst et. al. 2013 to obtain the
%      velocity of basilar membrane movement at different positions
%
%   2) the obtained velocity information is half-wave rectified
%
%   3) the half-waves are replaced with Gaussian pulses centered around
%      the local maxima of the half-waves
%
%   4) the frequency-dependent delays of the cochlea model are
%      compensated
%
%   See also: takanen2013, verhulst2012
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
%     S. Verhulst, T. Dau, and C. A. Shera. Nonlinear time-domain cochlear
%     model for transient stimulation and human otoacoustic emission. J.
%     Acoust. Soc. Am., 132(6):3842 -- 3848, 2012.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/takanen2013_periphery.php


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


%% ------ Check the input arguments ---------------------------------------
if nargin<3
    outputPlot =0;
end

%specify the characteristic frequencies of the model
fc = erbtofreq(4:39);
spl = 60; % sound pressure level of the input signal 
%% ------ Process the binaural input with the cochlear model --------------

if isstruct(insig)
    % input to the periphery model can be also a output of the Verhulst
    % cochlear model
    cochlear = insig;
else
    norm_factor=max(abs(insig(:,1)));     %first channel as reference for rms normalization
    insig(:,1)=insig(:,1)./norm_factor;
    insig(:,2)=insig(:,2)./norm_factor;
    [V,Y,E,CF] = verhulst2012(insig,fs,fc,[spl spl]);
    cochlear.velocityLeft=V(:,:,1);
    cochlear.velocityRight=V(:,:,2);
end

%load a vector of frequency specific delays computed for the Verhulst model
%outputs at fc
x=amt_load('takanen2013','cochleardelays.mat');
cochlear.delaysV = x.velocitydelays;

%use the velocity of the basilar membrane movement as input to CN model
processed.left = cochlear.velocityLeft;processed.right = cochlear.velocityRight;
[nVals nBands] = size(processed.left);

%% ------ Half-wave rectification -----------------------------------------
processed.left= processed.left.*(processed.left >0);
processed.right= processed.right.*(processed.right >0);
periphOutput.ventralLeft = processed.right;
periphOutput.ventralRight = processed.left;

%% ------ Impulse generation ----------------------------------------------
%search of the local maximas of each half-wave separately for left and 
%right ear signals
for chanInd=1:nBands
    left = processed.left(:,chanInd);
    right = processed.right(:,chanInd);
    processed.left(:,chanInd) = zeros(nVals,1);
    processed.right(:,chanInd) = zeros(nVals,1);

    %start end end points for each half-wave
    startingPoints = strfind((left'>0),[0 1]);
    endpoints = [strfind((left'>0),[1 0])+1 length(left)];
    if (isempty(startingPoints) && ~isempty(endpoints))
        startingPoints = 1;
    else
        if startingPoints(1)>=endpoints(1)
            startingPoints = [1 startingPoints];
        end
    end
    
    rmsVals=zeros(size(startingPoints));
    locations = rmsVals;
    nsamples = endpoints(1:length(startingPoints))-startingPoints+1;
    %compute the rms values of each half-wave and position it at the local
    %maximum
    for locInd=1:length(startingPoints)
        rmsVals(locInd)= norm(left(startingPoints(locInd):endpoints(locInd)))/sqrt(nsamples(locInd));
        [unnecessary, locations(locInd)] = max(left(startingPoints(locInd):endpoints(locInd)));
    end
    processed.left(locations+startingPoints-1,chanInd) = rmsVals';

    %processing of the right channel
    
    %start end end points for each half-wave
    startingPoints = strfind((right'>0),[0 1]);
    endpoints = [strfind((right'>0),[1 0])+1 length(right)];
    if (isempty(startingPoints) && ~isempty(endpoints))
        startingPoints = 1;
    else
        if startingPoints(1)>=endpoints(1)
            startingPoints = [1 startingPoints];
        end
    end
    %compute the rms values of each half-wave and position it at the local
    %maximum
    rmsVals= zeros(size(startingPoints));locations = rmsVals;
    nsamples = endpoints(1:length(startingPoints))-startingPoints+1;
    for locInd=1:length(startingPoints)
        rmsVals(locInd)= norm(right(startingPoints(locInd):endpoints(locInd)))/sqrt(nsamples(locInd));
        [unnecessary, locations(locInd)] = max(right(startingPoints(locInd):endpoints(locInd)));
    end
    processed.right(locations+startingPoints-1,chanInd) = rmsVals';

end
%% ------ Convolution with Gaussian window-functions ----------------------
% the width of the Gaussian window in the peripheral hearing model depends
% on the center frequency
N = zeros(size(fc));
N(fc<800) = round((2./fc(fc<800))*fs);
indMid = find(((800<=fc).*(fc<=2800))==1);N(indMid) = round(0.0024*(0.6+0.4*fc(indMid)./800)*fs);
N(fc>2800) = round(0.0048*fs);
% the constant alpha is set to 20
alpha=20;

periphOutput.left = zeros(nVals,nBands);periphOutput.right = periphOutput.left;
for chanInd=1:nBands
    left = processed.left(:,chanInd);
    right = processed.right(:,chanInd);
    n = -N(chanInd)/2:1:N(chanInd)/2;
    winFunction = (exp(-0.5*(alpha*n/(N(chanInd)/2)).^2))';
    %convolution with the gaussian window function
    left = (1/sum(winFunction))*conv(left,winFunction,'same');
    right = (1/sum(winFunction))*conv(right,winFunction,'same');

    %% compensation for the cochlear model delays
    periphOutput.left(:,chanInd) = [left(cochlear.delaysV(chanInd)+1:end);zeros(cochlear.delaysV(chanInd),1)];
    periphOutput.right(:,chanInd) = [right(cochlear.delaysV(chanInd)+1:end);zeros(cochlear.delaysV(chanInd),1)];
    periphOutput.ventralLeft(:,chanInd) = [periphOutput.ventralLeft(cochlear.delaysV(chanInd)+1:end,chanInd);zeros(cochlear.delaysV(chanInd),1)];
    periphOutput.ventralRight(:,chanInd) = [periphOutput.ventralRight(cochlear.delaysV(chanInd)+1:end,chanInd);zeros(cochlear.delaysV(chanInd),1)];
end

periphOutput.fc = fc;
%% ------ Plot the periphery model output, if desired ---------------------
if outputPlot
    figure(90);
    g(1) = subplot(6,1,1);plot((0:length(insig(:,1))-1)./(fs/1000),insig(:,1));ylabel('Input');set(gca,'yTick',[]);
    g(2) = subplot(6,1,2);plot((0:size(periphOutput.left,1)-1)./(fs/1000),periphOutput.left(:,4));title([num2str(round(fc(4))) ' Hz']);
    g(3) = subplot(6,1,3);plot((0:size(periphOutput.left,1)-1)./(fs/1000),periphOutput.left(:,9));title([num2str(round(fc(9))) ' Hz']);
    g(4) = subplot(6,1,4);plot((0:size(periphOutput.left,1)-1)./(fs/1000),periphOutput.left(:,16));title([num2str(round(fc(16))) ' Hz']);
    g(5) = subplot(6,1,5);plot((0:size(periphOutput.left,1)-1)./(fs/1000),periphOutput.left(:,22));title([num2str(round(fc(22))) ' Hz']);
    g(6) = subplot(6,1,6);plot((0:size(periphOutput.left,1)-1)./(fs/1000),periphOutput.left(:,28));title([num2str(round(fc(28))) ' Hz']);
    xlabel('Time [ms]');
    linkaxes(g,'x');
    clear g
end


