function [APvec] = kelvasa2015_anprocessing(electrodogram, vTime, varargin)
%KELVASA2015_ANPROCESSING  AN model used in Kelvasa and Dietz 2015 binaural model
%   Usage: [APvec] = kelvasa2015_anprocessing(electrodogram, vTime);
%
%   Input parameters:
%       electrodogram : N x M matrix of electrode current values in (microA) 
%                       with N being the number of CI electrodes and M being
%                       a time vector with 1/pulse rate/maxima sampling frequency 
%
%       vTime         : time vector in seconds with M samples
%
%
%   Output parameters:
%       APvec         : N x 2 matrix of AN spikes with Nx1 holding indices 
%                       of the spiking neuron and Nx2 holding corresponding
%                       spike time in seconds. 
%
%   KELVASA2015_anprocessing(insig,fs,varargin) computes auditory nerve
%   spike times over a given population of AN fibers using a simulated
%   electrode nerve interface as detailed in (Fredelake & Hohmann (2012))
%
%   References:
%     S. Fredelake and V. Hohmann. Factors affecting predicted speech
%     intelligibility with cochlear implants in an auditory model for
%     electrical stimulation. Hearing Research, 287(1):76 - 90, 2012.
%     [1]http ]
%     
%     V. Hamacher. Signalverarbeitungsmodelle des elektrisch stimulierten
%     GehÃ¶rs; 1. Aufl. PhD thesis, RWTH Aachen, Aachen, 2004. Zugl.: Aachen,
%     Techn. Hochsch., Diss., 2003.
%     
%     D. Kelvasa and M. Dietz. Auditory model-based sound direction
%     estimation with bilateral cochlear implants. Trends in Hearing,
%     19:2331216515616378, 2015.
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S0378595512000639
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/kelvasa2015_anprocessing.php

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

%          
%   Authors: 
%            Daryl Kelvasa (daryl.kelvasa@uni-oldenburg.de) 2016
%            Mathias Dietz (mdietz@uwo.ca) 2016
%            Stefan Fredelake
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Model Parameters
definput.import={'kelvasa2015'};
[~,kv]  = ltfatarghelper({},definput,varargin);

%% Check input paramters

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if size(electrodogram,2)~= numel(vTime)
  error('%s: Unequal number of time samples in electrodogram and time vector. ',upper(mfilename));
end;

if size(electrodogram,1)~= numel(kv.X_EL)
  error('%s: Unequal number of electrodogram electrodes and specified electrode locations. ',upper(mfilename));
end;


%% Implement AN model
%Create vector of electrode pulse length in seconds
numberNeuron = kv.N_nervecells;
numberCycle = numel(vTime);
Tph = ones(numberNeuron,1).*kv.Tph;

%%Prepare Membrane Electrical properties
tauM = kv.TAUCHR/log(2);  % the membrane time constant     (eq. 6.10)
R    = tauM/kv.C;         % the resistance of the membrane (eq. 6.10)


%%Prepare voltage threshold
Uth  = kv.EFFIRHEO .* R;

%%Prepare matrix of random correlated refractory values
[T_ARP,tau_RRP] = calculate_refracConstants(kv.MT_ARP,kv.MTAU_RRP,numberCycle);

%%Prepare membrane noise
UN = membraneNoise(numberNeuron,1/vTime(2),numberCycle); %1/vTime(2) is the sampling rate of the stimulation pattern
RS0 = kv.RS0_ind * (1 + 792*Tph(1) - 65833*Tph(1)^2); % phase dependent relative spread 

%%Apply spatial spread function to electrodogram (eq. 1)
effIamp = kv.V * electrodogram;

%%Calculate depolarisation potential of the cell membrane (eq. 2).
UD = effIamp .* repmat(R .* (1-exp(-Tph./(kv.TAUCHR/log(2)))),...
                        1,size(effIamp,2));
                                        
%%Calculate action potentials
% Init values
numberCycle = numel(vTime);
tLAP   = ones(numberNeuron,1)*-99/1000;  % arbtrary init value for the last action potential
tAP   = ones(numberNeuron,1)*-99/1000;  % arbtrary init value for the last action potential


APvec = [];
% Berechnung der Aktionspotentiale
for iPulse = 1:numberCycle
                            
    r = calculateRefractoryFunction(vTime(iPulse)-tLAP,T_ARP(:,iPulse),...
                                                    tau_RRP(:,iPulse)); 

    % correct the relative spread with the refractory function. the higher r,
    % the less is the relative spread
    RS0     = RS0./r;      % (eq. 6.30)

    % calculate the std of the noise
    sigma_Tph = Uth.*RS0; 

    % adjust the noise samples with the standard deviation. 
    vUN        = UN(:,iPulse) .* sigma_Tph;

    % calculate the probability of firing (eq. 6.28)
    P_AP = 1/2 + 1/2 * erf((UD(:,iPulse) - r.*Uth) ./ (sqrt(2)*sigma_Tph));

    % find for all neurons the depolarization current greater than the 
    % threshold current, multiplied with the refractory function and noised the 
    % noise samples (eq. 6.27)
    indexNeuron = 1:kv.N_nervecells;
    bool_AP   = (UD(:,iPulse) +vUN >= (r .* Uth) );
    bool_noAP = ~bool_AP;

    index_AP   = indexNeuron(bool_AP);
    index_noAP = indexNeuron(bool_noAP);

    % the time of the new action potentials... 
    tp = vTime(iPulse)*ones(kv.N_nervecells,1);

    if ~isempty(index_AP)

        tAP(index_AP) =  tp(index_AP) + ...
                        kv.TAUCHR(index_AP) .* ...
                        log2(effIamp(index_AP,iPulse)./(effIamp(index_AP,iPulse) ...
                        - (1-vUN(index_AP)./Uth(index_AP)).*kv.EFFIRHEO(index_AP)));% eq. 6.17
          
        tAP(index_AP(~isreal(tAP(index_AP)))) =   vTime(iPulse); % case that log2 is negative
        tAP(index_AP(isinf(tAP(index_AP))))   =   vTime(iPulse); % case that neuron fires with no input

        [dj] = calculate_latency_dj(kv.ML50,kv.SIGMAL50,P_AP);
        tAP(index_AP)   = tAP(index_AP) + dj(index_AP);
        tLAP(index_AP)   = tAP(index_AP);
   
        APvec(end+1:end+length(index_AP),:)= [index_AP' tAP(index_AP)]; %line from Mathias [#nervecell spiketime]
    end
    
    % Keep refrac constants constant for neurons that didn't fire
    T_ARP(index_noAP,iPulse+1) = T_ARP(index_noAP,iPulse);
    tau_RRP(index_noAP,iPulse+1) = tau_RRP(index_noAP,iPulse);   
end
 
end

%% Helper Functions
%%
function n=membraneNoise(N_nervecells,fs,N_lenNoise)
% usage:  n=membraneNoise(N_nervecells,fs,N_lenNoise)
% input:  N_nervecells = number of AN fibers 
%         fs = sampling frequency
%         N_lenNoise = length of noise in samples
%
% output: n   = membrane noise
%
% Generate membrane noise of N elements at a rate of fs [1/s].
%
% The output n will be a row-vector.
% (Calculation is speeded up for fs=7200 and fs=9000)
%
% Author: Stefan Fredelake
% Date: 02-12-2008
% Copyright (C) 2008 Stefan Fredelake, Oldenburg University


% fc: Butterworth corner frequency in Hz
fc= 200;

% fo: Butterworth filter order
%fo= 1;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% generate Gaussian noise with given noise power and length
n = randn(N_nervecells,N_lenNoise);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if fs > 400 % only filter the Gaussian noise if its sampling frequency is more than twice the desired corner frequency
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % design the filter coefficients + speed up with default parameters
    if (fc == 200)
        if (fs == 7200)
            b=[0.0805,0.0805];
            a=[1.0000,-0.8391];
        elseif (fs == 9000)
            b=[0.0654,0.0654];
            a=[1.0000,-0.8693];
        else
            [b,a]= butter(1,fc/(fs/2));
        end;
    end;
    
    
    if (size(n,1)*size(n,2)) > 5e7 %this would be a too large matrix for the filter command
        %split it up into 10 packs of cells and filter them separately
        number_of_packs = 10;
        newvec = zeros(N_nervecells,N_lenNoise);
        for iCounter = 1:number_of_packs
            newvec(1+(iCounter-1)*size(n,1)/number_of_packs:iCounter*size(n,1)/number_of_packs,:) ...
                = filter(b,a,n(1+(iCounter-1)*size(n,1)/number_of_packs:iCounter*size(n,1)/number_of_packs,:),[],2);
        end
    else
        % filter noise
        n= filter(b,a,n,[],2);
    end
end

end

%%
function [T_ARP,tau_RRP] = calculate_refracConstants(MT_ARP,MTAU_RRP,numberCycle)

% usage:  [T_ARP,tau_RRP] = calculateT_ARP_tauRRP(N)
% input:  N   = number of neurons
% 
% output: T_ARP   = absolute refractory phase
%         tau_RRP = relative refractory phase
% 
% generates a random vector with absoulte and relative refractory phases.
% both vectors are correlated with rho = 0.75
% 
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 02-12-2008
% Copyright (C) 2008 Stefan Fredelake, Oldenburg University
%

N = size(MT_ARP,1);
corrCoef = 0.75;
R = [1 corrCoef; corrCoef 1];
L = chol(R);

T_ARP = zeros(N,numberCycle);
tau_RRP = zeros(N,numberCycle);

std_mT_ARP   = 0.15 * repmat(MT_ARP,1,numberCycle);
std_mtau_RRP = 0.15 * repmat(MTAU_RRP,1,numberCycle);

x0   = randn(N,numberCycle);
x1 = x0 .* std_mT_ARP;
x2 = x0 .* std_mtau_RRP;

for ind = 1 : numberCycle

    temp = [x1(:,ind),x2(:,ind)] * L;
    T_ARP(:,ind)   = MT_ARP   +  temp(:,1);
    tau_RRP(:,ind) = MTAU_RRP +  temp(:,2);
end

end


%%
function r = calculateRefractoryFunction(t,T_ARP,tau_RRP)

% usage:  r = calculateRefractoryFunction(t)
% input:  t       = time after the last action potential
%         T_ARP   = absolute refractory phase
%         tau_RRP = relative refractory phase
%
% output: r = refractory factor, which will be multiplied with Uth
% 
% calculates the refractroy factor, which is needed for the raise of the
% threshold Uth. 
% 
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 02-12-2008
% 
q = 0.1; % constant from Hamacher p. 79
p = 0.68;  % constant from Hamacher p. 79

t = ones(size(T_ARP)).*t;
r = (1-exp(-(t-T_ARP)./(q*tau_RRP))).*(1-p*exp(-(t-T_ARP)./(tau_RRP)));
r = 1./r;
r(t<=T_ARP) = Inf;
 
end


%% 
function [dj] = calculate_latency_dj(mL50,sigmaL50,P_AP)

% usage:  [dj] = calculate_latency_dj(mL50,sigmaL50,P_AP)
% input:  mL50     = 
%         sigmaL50 = 
%         P_AP     = 
% 
% output: dj = latency
%       
% calculates the latency
% 
% Reference: Hamacher 2003 Ph.D.-thesis
%
% Author: Stefan Fredelake
% Date: 17-09-2009

N = size(mL50,1);

s1 = -248e-6;
s2 = - 79e-6;

md = mL50     + s1*(P_AP-0.5);
sd = sigmaL50 + s2*(P_AP-0.5); 

dj = abs(randn(N,1) .* sd + md);

end







