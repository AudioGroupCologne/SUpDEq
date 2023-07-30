function [spikeRatePerNeuron,spikeRatePerBin] = ...
                        kelvasa2015_anbinning(APvec,sigLengthSec,...
                                                    varargin)
%KELVASA2015_ANBINNING  AN and time binning from Kelvasa and Dietz 2015 binaural model
%   Usage: [spikeRatePerNeuron,spikeRatePerBin] = ...
%                        kelvasa2015_anbinning(APvec,sigLengthSec);
%
%   Input parameters:
%       APvec         : N x 2 matrix of AN spikes with Nx1 holding indices 
%                       of the spiking neuron and Nx2 holding corresponding
%                       spike time in seconds. 
%
%       sigLengthSec  : length of input signal in seconds
%
%   Output parameters:
%       spikeRatePerNeuron: N x M matrix of AN spike rates in spikes/second
%                           with N being the number of user defined AN 
%                           fibers and M being the number of time windows.
%
%       spikeRatePerBin   : N x M matrix of AN spike rates in spikes/second
%                           with N being the number of user defined AN fibe
%                           bands and M being the number of time windows. 
%
%   KELVASA2015_anbinning(APvec,sigLengthSec,varargin) bins auditory nerve 
%   spike times over a given population of AN fibers into user defined AN 
%   frequency bands and time bins as detailed in (Kelvasa & Dietz (2015))
%
%   References:
%     D. Kelvasa and M. Dietz. Auditory model-based sound direction
%     estimation with bilateral cochlear implants. Trends in Hearing,
%     19:2331216515616378, 2015.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/kelvasa2015_anbinning.php


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

%Retrieve and compute model paramters
definput.import={'kelvasa2015'};
[~,kv]  = ltfatarghelper({},definput,varargin);
                    
%% Main code

%initialize variables
numWin = ceil(sigLengthSec/kv.timeWindowSec);
winEdges = linspace(0,sigLengthSec,numWin+1);
spikeRatePerNeuron = zeros(kv.N_nervecells,numWin);
spikeRatePerBin = zeros(kv.numBin,numWin);

if ~isempty(APvec)
    
    [~,ind] = histc(APvec(:,2),winEdges);
    ind(ind==numWin+1) = numWin;
    
    for win = 1 : numWin
   
        APwin = APvec(ind == win,:);
        [H, ~] = histc(APwin(:,1),1:kv.N_nervecells);
        clear APwin
        spkRate =  H./kv.timeWindowSec;
        spikeRatePerNeuron(:,win) = spkRate;
        spikeRatePerBin(:,win) = mean(spkRate(kv.binPosInd),2);
    
    end
end

end


