function definput=arg_kelvasa2015(definput)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for ACE CI signal processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_kelvasa2015.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
definput.keyvals.devicename = 'nucleus';

definput.keyvals.NumOfChannels = 22;                                       % Number of electrodes

definput.keyvals.vActiveElectrodes =(definput.keyvals.NumOfChannels:-1:1)';% Vector which contains the numbers of the active electrodes after processing. Initialized with all possible Electrodenumbers.

definput.keyvals.NUMOFCHANNELS = ...                                       % BW:(???) Not sure why this is needed. Might allow "disabled" CI-Electrodes by placing "nan" in vActiveElectrodes?
                              length(definput.keyvals.vActiveElectrodes...
                             (~isnan(definput.keyvals.vActiveElectrodes)));

definput.keyvals.X_EL = 8.125:0.75:23.875;                                 % position of the electrodes (mm)

definput.keyvals.TCL    = [100, 100, 100, 100, 100, 100, 100, 100, ...     % MAP-T-Levels (in clinical units, denotes hearing threshold with CI)
                                100, 100, 100, 100, 100, 100, 100, 100,...
                                100, 100, 100, 100, 100, 100]';

definput.keyvals.MCL    = [200, 200, 200, 200, 200, 200, 200, 200, ...     % MAP-C-Levels (in clinical units, denotes maximum loudness threshold with CI)
                                200, 200, 200, 200, 200, 200, 200, 200, ...
                                    200, 200, 200, 200, 200, 200]';

definput.keyvals.Tph    = 25e-6;                                           % rectangular pulse phase duration [s]

definput.keyvals.ipg    =  8e-6;                                           % interphase gap between up- and down rectangular pulse path

definput.keyvals.T_SPL  = 25.0;                                            % The dB SPL Value, which is mapped to the T-Value [threshold db SPL]

definput.keyvals.C_SPL  = 65.0;                                            % The dB SPL Value, which is mapped to the C-Value [maximum loudness db SPL]

definput.keyvals.FS_ACE = 16000;                                           % Sampling rate of the ACE-Processing (Fs from digital sampling of the microphone or fs from the output of the filterbank?)

definput.keyvals.NFFT = 128;                                               % Sets the number of fft-bins for analysis window

definput.keyvals.pps = 900;                                                % Sets the pulses-per-second rate of the Cochlea Implant

definput.keyvals.maxima = 8;                                               % maximal number of electrodes in n-of-m strategy

definput.keyvals.indexOrg = 2;

%Filter bank
definput.keyvals.Q_SUM = [  ones(9,1); ...                                 %Information, how to sum up the frequency bins (num. f-Bins in table below)
                            2*ones(4,1); ...
                            3*ones(2,1); ...
                            4*ones(2,1); ...
                            5*ones(2,1);
                            6;7;8];

definput.keyvals.vLowerFreq = [ 188  313  438  563  688  813  938 1063 ... %Upper and lower frequencies for each electrode
                                1188 1313 1563 1813 2063 2313 2688 3063 ...
                                3563 4063 4688 5313 6063 6938]';

definput.keyvals.vUpperFreq = [ 313  438  563  688  813  938 1063 1188 ...
                                1313 1563 1813 2063 2313 2688 3063 3563 ...
                                4063 4688 5313 6063 6938 7938]';

% calculate the weights for the fft from the dissertation of Brett Swanson (2008)
W    = hann(definput.keyvals.NFFT);
W1 = abs(fft(W));                                                          % for G(2), which is derived from W(0.5)
W2 = abs(fft([W; zeros(definput.keyvals.NFFT,1)]));

definput.keyvals.G(1) = 2/W1(1);                                           % for channels with 1 bin

definput.keyvals.G(2) = 2/(sqrt(2)*W2(2));                                 % for channels with 2 bins

definput.keyvals.G(3) = 2/78.38;                                           % for channels with more than 2 bins.

%compression characteristics
definput.keyvals.B = 0.0156;                                               % base level

definput.keyvals.M = 0.5859;                                               % saturation level

definput.keyvals.alpha_c = 415.96;                                         % controls the steepness of the compression function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for electrode nerve interface and AN model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Set parameters for auditory nerve model
definput.keyvals.N_nervecells = 500;                                       % Number of simulated auditory nerve fibers

definput.keyvals.X_NZ = linspace(0,35,definput.keyvals.N_nervecells)';     % position of the nerve cell on the basilar membrane (mm)

%%Spatial spread function from Hamacher: p. 95

definput.keyvals.lambda = 9;                                               % constant of the spatial spread function (mm)

definput.keyvals.v0 = 1;                                                   % amplitude of the spatial spread function

%%Membrance capacitance
definput.keyvals.C = 1;                                                    % (???) Value doesnt't seem to be important DK

%%Spatial spread
definput.keyvals.V = calculate_AVF(definput.keyvals.X_EL,...
                                    definput.keyvals.X_NZ,...
                                    definput.keyvals.lambda,...
                                    definput.keyvals.v0);

%%Get current output based upon CI model
definput.keyvals.Ith0 = conversion_CL2Iamp(definput.keyvals.TCL,...
                                            definput.keyvals.devicename);

%%Refractory constants from Hamacher p. 79
[definput.keyvals.MT_ARP,definput.keyvals.MTAU_RRP] = ...
             generateRefractoryConstants(definput.keyvals.N_nervecells);

%%Neural latency constants from Hamacher p. 102
[definput.keyvals.ML50,definput.keyvals.SIGMAL50] =...
             generateNeuralLatencyConstants(definput.keyvals.N_nervecells);

%Chronaxie and rheobase constants according to Hamacher pp 97-100.
[definput.keyvals.TAUCHR,definput.keyvals.EFFIRHEO] = ...
             generateChronaxieRheobase(definput.keyvals.N_nervecells);

%%Relative spread constants(note 1) from Hamacher p. 100
definput.keyvals.RS0_ind = ...
            calculateRelativeSpread(definput.keyvals.N_nervecells);


%Reference current
definput.keyvals.IREF = ...                                                %Calculates the deteministic current for threshold for a given pulse
                          calculateDeteminsticThresholdCurrent_2(...       %with a defined phase kv.Tph
                                     mean(definput.keyvals.EFFIRHEO),...
                                     mean(definput.keyvals.TAUCHR/log(2)),...
                                     definput.keyvals.Tph);

for ii = 1:length(definput.keyvals.Ith0)                                   % (???) hier empirische anpassung des CIs ans Hörnervenmodell, dabei soll die
definput.keyvals.V(:,ii) = definput.keyvals.V(:,ii).*...                   % Höhe von V so variiert werden, dass SR, ca, den Wert von 30 annimmt.
                    definput.keyvals.IREF/definput.keyvals.Ith0(ii,:);     % Entspricht in etwa 30 APs.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for binning stages and localization models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Paramters for binning nervogram into time and AN rate bins
definput.keyvals.binPos  = [ 10, 13, 16, 19, 22];                          %AN frequency band positions in mm along cochlea from apex

definput.keyvals.numNeuronsInBin  = 20;                                    %number of neurons in each AN frequency band. Must be 20 for LSO model

definput.keyvals.timeWindowSec = 0.2;                                      %time binning of entire input wavfile in seconds

%Choosing and calibration of localization model
%Signal params
definput.keyvals.azis = 0:5:90;                                            %vector of azimuthal source locations in degrees over which to calibrate the mapping

definput.keyvals.HRTFfile = 'Kayser2009_Anechoic.sofa';                    %Name of *.sofa file containing HRTF data

definput.keyvals.HRTFelevation = 0;                                        %desired elevation used in HRTF filtering

definput.keyvals.HRTFsourceDistance = 300;                                 %desired source distance used in HRTF filtering

definput.keyvals.HRTFchannels = [1,2];                                     %desired microphone output from HRTF filtering


%Localization model params
definput.keyvals.localizationModelCalibWav = 'SSN.wav';                    %wavfile must be in MATLAB path

definput.keyvals.localizationModelCalibStimulusLevelDB =  55;              %in model dB

definput.keyvals.identifier = 'inEar';                                     %Use to identify calibrations processed with different paramters. ie BTE = behind the
                                                                           %ear microphones.

definput.keyvals.localizationModel =  'ResponseDifferenceAN';              %choose localization model from Kelvasa and Dietz 2015. Can be
                                                                           %'RateLevel','ResponseDifferenceAN','MaxLikelihood'

definput.keyvals.dBRange = [35:75];                                        %range in dB over which to process calibration signal to extract rate level slopes

definput.keyvals.RateLevelCrvfitRange =  [35,75];                          %range in dB over which to calculate rate level slopes

definput.keyvals.AziCrvfitRange =  45;                                     %range in degrees azimuth over which to calculate ILD/angle and spikeDiff/angle slopes

binPosInd = round(definput.keyvals.binPos.*(1/definput.keyvals.X_NZ(2)));  %compute indices of AN fibers
if mod(definput.keyvals.numNeuronsInBin,2) == 0
    binPosInd = [binPosInd-(definput.keyvals.numNeuronsInBin/2);
                    binPosInd+(definput.keyvals.numNeuronsInBin/2)-1];
else binPosInd = [binPosInd-floor(definput.keyvals.numNeuronsInBin/2);
                   binPosInd+floor(definput.keyvals.numNeuronsInBin/2)];
end

definput.keyvals.numBin = numel(definput.keyvals.binPos);

for bin = 1 : definput.keyvals.numBin
                                   definput.keyvals.binPosInd(bin,:) = ...
                                   [binPosInd(1,bin):binPosInd(2,bin)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check_definput(definput.keyvals);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Helper Functions
% Author: Daryl Kelvasa
% Date: 01-2-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = calculate_AVF(x_el,x_nz,lambda,v0)

% usage:  v = calculate_AVF(x_el,x_nz,szWhichConfig)
% input:  x_el = position of first electrode [mm]
%         x_nz = position of neuron cell [mm]
%         lambda = constant of the spatial spread function
%           v0 = constant factor of the spatial spread
%
% output: v = spreading function
%
% calculates the spreading function of an electrical field according for
% monopolar stimulation
%
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 17-10-2008
% Copyright (C) 2008 Stefan Fredelake, Oldenburg University

v = zeros(length(x_nz),length(x_el));

    for iElectrode = 1:length(x_el)
        v(:,iElectrode) = abs(exp(-abs(x_nz-x_el(iElectrode))/lambda));
    end

v  = v0 * v;
end

%%
function [mT_ARP,mTau_RRP] = generateRefractoryConstants(N)

% usage:  [mT_ARP,mTau_RRP] = generateRefractoryConstants(N)
% input:  N   = number of neurons
%
% output: T_ARP   = absolute refractory phase
%         tau_RRP = relative refractory phase time constant
%
% generates a random vector with absoulte and relative refractory phases.
% both vectors are correlation with rho = 0.75
%
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 02-12-2008
% Copyright (C) 2008 Stefan Fredelake, Oldenburg University

mean_mT_ARP   = 0.7e-3; % mittlere absolute Refraktärphase
mean_mtau_RRP = 1.6e-3; % mittlere relative Refraktärphase

std_mT_ARP   = 0.15 * mean_mT_ARP;
std_mtau_RRP = 0.15 * mean_mtau_RRP;

x0   = randn(N,1);
xnew = (randn(N,1)*std_mtau_RRP) + mean_mtau_RRP;

x1 = (x0*std_mT_ARP)   + mean_mT_ARP;
x2 = (x0*std_mtau_RRP) + mean_mtau_RRP;


rho  = my_corr(x1,x2);
while rho > 0.75
    index = 0;
    while index == 0
        index = round(N*rand(1));
    end
    x2(index) = xnew(index);
    rho = my_corr(x1,x2);
end

mT_ARP   = x1;
mTau_RRP = x2;
end


%%
function [mL50,sigmaL50] = generateNeuralLatencyConstants(N)

% usage:  [mL50,sigmaL50] = generateNeuralLatency(N)
% input:  N   = number of neurons
%
% output: mL50 =     random number for the calculation of the latency, each
%                    element stands for a neuron
%         sigmaL50 = random number for the calculation of the latency, each
%                    element stands for a neuron
%
%
% generates a random vector with numbers for the calculation of the
% propagation latency
%
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 17-09-2009
% Copyright (C) 2009 Stefan Fredelake, Oldenburg University

mean_mL50     = 607e-6;  % p.102
mean_sigmaL50 = 106e-6;  % p.102

std_mL50     = 142e-6;   % p.102
std_sigmaL50 =  87e-6;   % p.102


x0 = randn(N,1);
x1 = (x0*std_mL50)     + mean_mL50;
x2 = (x0*std_sigmaL50) + mean_sigmaL50;

xnew = (randn(N,1)*std_sigmaL50) + mean_sigmaL50;

rho  = my_corr(x1,x2);
while rho > 0.5
    index = 0;
    while index == 0
        index = round(N*rand(1));
    end
    x2(index) = xnew(index);
    rho = my_corr(x1,x2);
end

mL50     = x1;
sigmaL50 = xnew;
end

%%
function [tauChr,effIRheo] = generateChronaxieRheobase(N)

% usage:  [tauChr,IRheo] = generateChronaxieRheobase(N)
% input:  N = number of neurons
%
% output: tauChr = chronaxie values
%         IRheo  = effecitve rheobase values
%
% Generates random values for chronaxie and rheobase according to Hamacher
% pp 97-100.
%
%
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 02-12-2008
% Copyright (C) 2008 Stefan Fredelake, Oldenburg University

meanTauChr = 255e-6;    % seconds
stdTauChr  =  57e-6;    % seconds

meanIRheo  =   32e-6;    % ampere
stdIRheo   =    6.5e-6;  % ampere

tauChr     = (randn(N,1)*stdTauChr) + meanTauChr;
effIRheo   = (randn(N,1)*stdIRheo)  + meanIRheo;

end

%%
function Ith0 = calculateDeteminsticThresholdCurrent_2(IRheo,tauM,Tph)

% usage:  Ith0 = calculateDeteminsticThresholdCurrent(IRheo,C,Tph)
% input:  IRheo = the rheobase
%         tauM  = time constant of the membrane given with RC
%         Tph   = phase of the stimulus
%
% output: Ith0 = the deterministic current for threshold
%
% calculates the deteministic current for threshold for a given pulse
% with a defined phase Tph
%
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 01-12-2008
% Copyright (C) 2009 Stefan Fredelake, Oldenburg University

    Ith0 = IRheo./(1-exp(-Tph./tauM));
end

%%
function I = conversion_CL2Iamp(cl,szWhichDevice)
% Author: Stefan Fredelake
% Date: 02-12-2008
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

%%
function [RS0_ind] = calculateRelativeSpread(N)

% usage:  RS0 = calculateRelativeSpread(Tph,N)
% input:  Tph = phase duration in seconds
%         N   = number of neurons
%
% output: RS0 = relative spread accoring to Bruce et al
%
% calculates the relative spread accorind to Bruce et al. Note Bruce et al
% assumes at this function is only applicable to 100us<=Tph<=5000us and
% bipolar stimulation. (p. 38) But Hamacher used it for lower Tph values
% and also for monopolar stimulation. (p. 79)
%
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 02-12-2008
%
% RS0_ind = 0.12;   % phase duration independent part of the relative spread

meanRS0_ind = 0.095; % from Hamacher p. 100
stdRS0_ind  = 0.04;  % from Hamacher p. 100

RS0_ind     = (randn(N,1)*stdRS0_ind)+meanRS0_ind;

end

%%
function rho = my_corr(x,y)
    x = x - my_mean(x);
    y = y - my_mean(y);
    rho = sum(x.*y)/sqrt((sum(x.^2)*sum(y.^2)));
end

%%
function xm = my_mean(x)
    xm = sum(x)/length(x);
end


%%
%check ACE-Constants for sensible input
function check_definput(kv)
% MHL and THL should have the same length
if size(kv.NumOfChannels) ~= size(kv.X_EL)
    error('Electrode positions do not correspond to number of electrodes. Please check the Parameterfile!');
end

% MHL and THL should have the same length
if size(kv.MCL) ~= size(kv.TCL)
    error('There must be the same number of M-Levels (kv.MCL) and T-Levels (kv.TCL). Please check the Parameterfile!');
end

% T levels should be smaller than M Levels
if kv.TCL >= kv.MCL
    error('T-Levels (kv.TCL) must be smaller than M-Levels (kv.MCL). Please check the Parameterfile!');
end

% kv.maxima should be smaller than kv.NumOfChannels
if kv.maxima > kv.NumOfChannels
    error('In a n- of m-Strategie n (kv.maxima) must be smaller than m (kv.NumOfChannels). Please check the Parameterfile!');
end

% kv.T_SPL should be smaller than kv.C_SPL
if kv.T_SPL >= kv.C_SPL
    error('T-Levels (kv.T_SPL) must be smaller than C-Levels (kv.C_SPL). Please check the Parameterfile!');
end

% Currently only exactly 22 Electrodes are supported
if kv.NumOfChannels ~= 22
    error('Electrode Number (kv.NumOfChannels) must be 22. Other Electrode numbers are noch supported at this time. Please check the Parameterfile!');
end

%Check that AN frequency bands are bound by cochlea
if sum(kv.binPosInd(:)<0) > 0 || sum(kv.binPosInd(:)> kv.N_nervecells) > 0
    error('Specified bin positions are wider than cochlea.');
end
end

% Information about the filterbank.
% Chn.No. | num. f-Bins | fc, low bin | fc, upp bin
%      1  |           1 |         250 |         250
%      2  |           1 |         375 |         375
%      3  |           1 |         500 |         500
%      4  |           1 |         625 |         625
%      5  |           1 |         750 |         750
%      6  |           1 |         875 |         875
%      7  |           1 |        1000 |        1000
%      8  |           1 |        1125 |        1125
%      9  |           1 |        1250 |        1250
%     10  |           2 |        1375 |        1500
%     11  |           2 |        1625 |        1750
%     12  |           2 |        1875 |        2000
%     13  |           2 |        2125 |        2250
%     14  |           3 |        2375 |        2625
%     15  |           3 |        2750 |        3000
%     16  |           4 |        3125 |        3500
%     17  |           4 |        3625 |        4000
%     18  |           5 |        4125 |        4625
%     19  |           5 |        4750 |        5250
%     20  |           6 |        5375 |        6000
%     21  |           7 |        6125 |        6875
%     22  |           8 |        7000 |        8000

% Copyright (C) 2016   AG Medizinische Physik,
%                      Universitaet Oldenburg, Germany
%                      http://www.physik.uni-oldenburg.de/docs/medi
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


