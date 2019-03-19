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
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/kelvasa2015_anbinning.php

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
         
%   Authors: 
%            Daryl Kelvasa (daryl.kelvasa@uni-oldenburg.de) 2016
%            Mathias Dietz (mdietz@uwo.ca) 2016

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

