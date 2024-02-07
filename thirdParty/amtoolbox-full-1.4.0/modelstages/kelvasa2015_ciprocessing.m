function [electrodogram, vTime] = ...
            kelvasa2015_ciprocessing(insig,fs,varargin)
%KELVASA2015_CIPROCESSING CI ACE processing strategy used in Kelvasa and Dietz 2015 binaural model
%   Usage: [electrodogram, vTime] = kelvasa2015aceprocessing(insig,fs);
%
%   Input parameters:
%       insig       : single channel column vector signal
%       fs          : sampling rate (Hz) 
%
%   Output parameters:
%       electrodogram : N x M matrix of electrode current values in 
%                       (microAmpere) with N being the number of CI 
%                       electrodes and M being a time vector with 
%                       1/pulse rate/maxima sampling period 
%
%       vTime         : time vector in seconds with M samples 
%
%   KELVASA2015_ciprocessing(insig,fs,varargin) simuluates the ACE
%   CI signal processing strategy for an input signal using user defined 
%   parameters. 
%
%   References:
%     S. Fredelake and V. Hohmann. Factors affecting predicted speech
%     intelligibility with cochlear implants in an auditory model for
%     electrical stimulation. Hearing Research, 287(1):76 -- 90, 2012.
%     [1]http ]
%     
%     B. A. Swanson. Pitch perception with cochlear implants. PhD thesis, The
%     University of Melbourne, 2008.
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S0378595512000639
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/kelvasa2015_ciprocessing.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal M-Stats
%   #Author: Daryl Kelvasa (2016)
%   #Author: Mathias Dietz (2016)
%   #Author: Clara Hollomey (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input paramters 
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig) || size(insig,2) > size(insig,1)
  error('%s: insig has to be a numeric column vector signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

definput.import={'kelvasa2015'}; % defaults from arg_kelvasa2015
[~,kv]  = ltfatarghelper({},definput,varargin);

%% Main CI processing

%Initialize parameters
insig = process_preemphasis_CI(insig, kv.FS_ACE); % Hochpass
lenSignal = length(insig);
W = hann(kv.NFFT);

%%Compare pulse rate with signal sampling rate
switch kv.pps
    % muss zyklisch abgetastet werden
    case 900
        vorschub = [18 18 17 18 18 18 17 18 18];  % wegen rundungsfehler (16000/900 = 17.7778; mean(vorschub) = 17.7778)
    otherwise
        vorschub = round(kv.FS_ACE/kv.pps); 
        amt_disp('Warning: Potential rounding inaccuraries, because the signal sampling rate is not an integer multiplier of PPS');
end

%%Calculate windowed indices of signal
numVorschub = (lenSignal - kv.NFFT)/(sum(vorschub));
numFFTwin = (floor(numVorschub)*numel(vorschub)) + ...
                    ceil(rem(numVorschub,1)*numel(vorschub));
vorschub = [1, repmat(vorschub,1,ceil(numVorschub))];
winInd = zeros(kv.NFFT,numFFTwin); winInd(:,1) = (1:kv.NFFT);
for ind =  2 : numFFTwin;
    winInd(:,ind) = (0:kv.NFFT-1)+sum(vorschub(1:ind)); end
insig(max(winInd)) = 0;    

%%Process signal
    %Window signal
    x = insig(winInd); %clear winInd insig
    x = repmat(W,1,size(x,2)).*x;

    % spectra
    S = abs(fft(x,kv.NFFT));

%%Main loop over signal
electrodogram = zeros(kv.NumOfChannels,numFFTwin*kv.maxima);
pulseInd = 1 : kv.maxima;
vTime   = (0:size(electrodogram,2)-1).*(1/kv.pps/kv.maxima);
for ind = 1 : numFFTwin
    
    % "envelope" per time unit
    F = process_ACE_filterbank_demo(S(:,ind)', kv);
    
    % choose n of m
    [A,iElectrode] = process_nOfm(F,kv.maxima);
    iElectrode = sort(iElectrode,1,'descend');

    % compression between Tlevels and Clevels
    C = process_compression_ci(A,kv.B,kv.M,kv.alpha_c);
    cl = zeros(length(kv.vActiveElectrodes),1);
    logicalNoPulse = C(iElectrode) == 0;
    
    cl(iElectrode) = round(kv.TCL(iElectrode,:) + ...
                             (kv.MCL(iElectrode,:)- ...
                             kv.TCL(iElectrode,:)).*C(iElectrode,:));
    
    % current amplitude in micro-Ampere
    I = zeros(length(kv.vActiveElectrodes),1);
    I(iElectrode) = conversion_CL2Iamp(cl(iElectrode),kv.devicename);
    I(iElectrode(logicalNoPulse)) = 0;  
%     stimulationsequence = 1:length(kv.vActiveElectrodes);
%     iElectrode = stimulationsequence(iElectrode);
    
    electrodogram(sub2ind(size(electrodogram),...
                    iElectrode',pulseInd)) = I(iElectrode);
    pulseInd = pulseInd + kv.maxima;
end

end

%% Model helperfunctions
%% Pre-emphasis function
function y = process_preemphasis_CI(x, FS)
% Author: Stefan Fredelake
% Date: 23-10-2008

fc = 1200; 
w = 2*fc/FS;

[b,a] = butter(1,w,'high');

y = filter(b,a,x);

end

%% ACE filterbank
function F = process_ACE_filterbank_demo(S, CIParams)
%Code from the dissertation of Brett Swanson (2008)
% Author: Stefan Fredelake
% Date: 23-10-2008
%
% usage:  F = process_ACE_filterbank(S)
% input:  S = complex spectrum
%         CIParams = Struct containing the CI Params. In this function only
%         NumofChannels, Q_Sum, G, indexOrg and vActiceElectrodes are
%         needed. See set_global_constant.m for a description of the
%         mentioned parameters
% output: F = f�ltered output spectrum
% 
% This filterbank is only applicable for a complex spectrum S with a length
% of 128 bins, and a sampling frequency of 16 kHz. 

S2 = abs(S).^2; % square the amplitude spectrum


F  = zeros(CIParams.NumOfChannels,1);
for n = 1:CIParams.NumOfChannels
    if n == 1
        index2 = CIParams.indexOrg;  % only for the first loop
    end
    index1 = index2 + 1;
    index2 = index1 + CIParams.Q_SUM(n)-1;
    index  = index1:index2;
    if CIParams.Q_SUM(n) == 1
        weight = CIParams.G(1);
    elseif CIParams.Q_SUM(n) == 2
        weight = CIParams.G(2);
    else
        weight = CIParams.G(3);
    end
    F(n,:) = weight*sqrt(sum(S2(index)));   % eq. 5.4
end

end

%% NofM
function [A,iElectrode] = process_nOfm(F,nChns)
% Reference: Laneau, 2005 PhD-thesis
%
% Author: Stefan Fredelake
% Date: 23-10-2008
%
% usage:  A     = process_nOfm(F,nChns)
% input:  F          = filterbank output (22 filters)
% output: nChns      = number of channels to simulate
%         iElectrode = index of the stimulating electrodes
% 
% This function chooses the channels to stimulate. The channels with the
% nChns maximal values are stimulated

A               = zeros(size(F));
[~,iElectrode] = sort(F,1,'descend');
iElectrode      = iElectrode(1:nChns);
A(iElectrode)   = F(iElectrode);      % eq 5.16

end

%% Process Compression
function C = process_compression_ci(A,B,M,alpha_c)

% usage:  C = process_compression_ci(A)
% input:  A = the amplitude values in each channel to simulate
%         B = Base level
%         M = saturation level
%         alpha = controls the steepness of the function
% output: C = clinical current unit
% 
% This function compresses the envelopes into the electrical dynamic range
% of CI users. 

%Code from the dissertation of Brett Swanson (2008)
% Author: Stefan Fredelake
% Date: 23-10-2008

C = zeros(size(A));
C(A<B) = 0;
index = find(A>=B & A<M); 
C(index) = log( 1 + alpha_c * ((A(index)-B)/(M-B) ))/log(1+alpha_c);
C(A>=M) = 1; 

end

function I = conversion_CL2Iamp(cl,szWhichDevice)
% Reference: Laneau, 2005 PhD-thesis
%
% Author: Stefan Fredelake
% Date: 23-10-2008

switch szWhichDevice
    case 'nucleus'
        I = 10.175.^(cl/255); % aus laneau
        % I = I/24.838;
    case 'CI22'
        I = 100e-6 * 1.0205.^(cl-100); % aus hearcom website
    case 'CI24M'
        I = 100e-6 * 1.0205.^(cl-100); % aus hearcom website
    case 'CI24R'
        I = 100e-6 * 1.0205.^(cl-100); % aus hearcom website
    case 'freedom' % also CI24RE
        I = 100e-6 * 1.0182.^(cl-100); % aus hearcom website
    otherwise
        error('this device is not known')
end
end


