function [ANbinPredictions, weightedPredictions, mappingData] = ...
            kelvasa2015_localize(mappingData, SpkDiffPerBin,...
                                    SpkSumPerBin, varargin)
%KELVASA2015_LOCALIZE uses calibration data to map bilateral spike rate differences to an azimuthal angle
%   Usage: [ANbinPredictions, weightedPredictions, mappingData] = ...
%           kelvasa2015_localize(mappingData, SpkDiffPerBin, SpkSumPerBin, varargin)
%
%   Input parameters:
%       mappingData    : Structure containing the necessary data to map 
%                        bilateral spike rate differences to a predicted 
%                        azimuthal angle. 
%
%       SpkDiffPerBin  : N x M matrix of chan2 - chan1 spike rate
%                        differences in spikes/sec where N are user defined 
%                        AN frequency bands and M are time bins over the 
%                        entire signal. 
%
%       SpkSumPerBin   : N x M matrix of chan2 + chan1 spike rate
%                        sums in spikes/sec where N are user defined 
%                        AN frequency bands and M are time bins over the 
%                        entire signal. 
%
%   Output parameters:
%       ANbinPredictions   : N x M matrix of bin azimuthal angle predictions 
%                            in degrees with N being the number of user 
%                            defined AN frequency bands and M being the 
%                            number of time windows.
%
%       weightedPredictions: 1 x M vector of bin weighted azimuthal angle 
%                            predictions in degrees with M being the number
%                            of time windows.
%
%   Fields added to "mappingData":
%         RateLevelPerBin    : 
%                              N x M matrix of monaural spike rates in 
%                              spikes/sec with N being the number of user 
%                              defined AN frequency bands and M being the 
%                              range of signal levels in dB SPL over which
%                              the calibration signal was processed.
% 
%         ILDcentFreq        : 
%                              1 x N vector of center frequencies in Hz used 
%                              in CI processing that correspond to the N user 
%                              defined AN frequency bands.
% 
%         ILD                : 
%                              N x M matrix of interaural level differences in
%                              dB SPL for each of the the N user defined AN 
%                              frequency bands over the M azimuthal angles for 
%                              which the signal was processed.
%                              
%         IldAngSlope       : 
%                             1 x N vector of linear regression slopes in 
%                             dB/degree  for each of the the N user defined AN 
%                             frequency bands.
%                             
%         RateLevelSlope    : 
%                             1 x N vector of linear regression slopes in 
%                             spike Rate/dB  for each of the the N user defined 
%                             AN frequency bands.
%                             
%         spkDiffAziSlopes  : 
%                             1 x N vector of linear regression slopes in 
%                             bilateral spike rate differences /dB  for each 
%                             of the the N user defined AN frequency bands.
%                             
%         Avg               : 
%                             M x N vector of spike rate difference averages
%                             over all time bins where N are user defined AN 
%                             frequency bands and M are the range of azimuthal 
%                             angles over which the signal was processed.
%                               
%         covariance        : 
%                             M x N vector of spike rate difference covariances
%                             over all time bins where N are user defined AN 
%                             frequency bands and M are the range of azimuthal 
%                             angles over which the signal was processed.
%
%   KELVASA2015_localize(APvec,sigLengthSec,varargin)implements a user 
%   defined localization model as detailed in (Kelvasa & Dietz (2015)) to
%   map bilateral spike rate differences to a predicted azimuithal 
%   source angle.
%
%   References:
%     D. Kelvasa and M. Dietz. Auditory model-based sound direction
%     estimation with bilateral cochlear implants. Trends in Hearing,
%     19:2331216515616378, 2015.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/kelvasa2015_localize.php

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
%
% Copyright (C) 2016   AG Medizinische Physik,
%                      Universitaet Oldenburg, Germany
%                      http://www.physik.uni-oldenburg.de/docs/medi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input paramters
if nargin<4
  error('%s: Too few input parameters.',upper(mfilename));
end;

%Retrieve and compute model paramters
definput.import={'kelvasa2015'};
[~,kv]  = ltfatarghelper({},definput,varargin);

%% Map binaural spike rate differences to an azimuthal angle using the 
%specified localization model as per Kelvasa and Dietz 2015

if strcmp(kv.localizationModel, 'RateLevel')
    mappingData = computeRateLevelMapping(mappingData,kv);
    
    %1/(Rate/level * ILD/angle) = angle/rateDiff 
    slopes =1./(mappingData.IldAngSlope.*...
                    mappingData.RateLevelSlope);
    slopes = repmat(slopes',1,size(SpkSumPerBin,2));
    
    %This is where the magic happens
    predictions = slopes .* SpkDiffPerBin;
     
    %Limit model predictions
    ANbinPredictions = round(max(min(predictions,90),-10)./5).*5;
    
    %Compute spike sume weighted prediction over all bins
    weights = SpkSumPerBin./...
                            repmat(sum(SpkSumPerBin,1),kv.numBin,1);
                        
    weights(isnan(weights)) = 0; %for when there are no spikes
    weightedPredictions =   sum(ANbinPredictions.*weights,1);
    weightedPredictions = round(weightedPredictions./5).*5;
    mappingData.spkDiffAziSlopes = slopes;
    


elseif strcmp(kv.localizationModel, 'ResponseDifferenceAN')  
    mappingData = computeResponseDiffMapping(mappingData,kv);
    %Prepare linear mapping
    slopes =1./(mappingData.ResponseDiffSlope);
    slopes = repmat(slopes',1,size(SpkSumPerBin,2));
    
    %This is where the magic happens
    predictions = slopes .* SpkDiffPerBin;
     
    %Limit model predictions
    ANbinPredictions = round(max(min(predictions,90),-10)./5).*5;
    
    %Compute spike sume weighted prediction over all bins
    weights = SpkSumPerBin./...
                            repmat(sum(SpkSumPerBin,1),kv.numBin,1);
    weights(isnan(weights)) = 0; %for when there are no spikes                    
    weightedPredictions =   sum(ANbinPredictions.*weights,1);
    weightedPredictions = round(weightedPredictions./5).*5;
   
    

else % 'MaxLikelihood'
    %Compute statistical paramters for N (num of bins) dimensional PDFs
    Avg = zeros(numel(kv.azis),numel(kv.binPos));
    covariance = zeros(numel(kv.azis),numel(kv.binPos));
    
    for bin = 1 : numel(kv.binPos)
      data1 = squeeze(mean(...
      mappingData.calSpikeDiffPerNeuronPerAzi(:,kv.binPosInd(bin,:),:),2));
      Avg(:,bin) = mean(data1,2);   
      covariance(:,bin) = max(std(data1,[],2).^2,0.001);
    end
    
    mappingData.Avg = Avg;
    mappingData.covariance = covariance;

    weightedPredictions = computeMaxLikelihood(mappingData, ...
                                               SpkDiffPerBin, kv);
    
    %No bin predictions in this model
    ANbinPredictions = [];
end                      
end


%% Helper functions
%% Linear Rate Level localization model
function mappingData = computeRateLevelMapping(mappingData,kv)
%Initialize variables
    ILDcentFreq = zeros(1,numel(kv.binPos));

    ILD = zeros(numel(ILDcentFreq),...
                        numel(kv.azis));
                    
    IldAngSlope = zeros(1,numel(ILDcentFreq));
    
    RateLevelSlope = zeros(1,numel(ILDcentFreq));
    
    RateLevelPerBin = zeros(numel(kv.binPos),numel(kv.dBRange));
    
    for bin = 1 : numel(kv.binPos)
        
        %Find center frequencies of electrodes
        [~,ind] = min(abs(kv.X_EL - kv.binPos(bin)));
        ILDcentFreq(bin) = ((kv.vUpperFreq(ind) - kv.vLowerFreq(ind))/2)...
                            + kv.vLowerFreq(ind);
   
        %Compute ILD using gammatone filter 
            %Preset gamma filter parameters
            attenuation_db = 3;
            filter_order = 4;

            for ang = 1 : numel(kv.azis)
                
                %Get HRTF signal at angle
                HRTFsig = squeeze(mappingData.calibHRTFsig(ang,:,:));
                %ERB bandwidth
                bandwidth_hz = 24.7*(4.37*ILDcentFreq(bin)*0.001 + 1);
     
                %Implement gammatone filter
                filter = hohmann2002_filter(kv.FS_ACE, ILDcentFreq(bin),...
                    bandwidth_hz, attenuation_db, filter_order);
                
                [filtSigA, ~] = hohmann2002_process(filter,  HRTFsig(:,1));
              
                [filtSigB, ~] = hohmann2002_process(filter,  HRTFsig(:,1));
    
                %Compute ILD
                ILD(bin,ang) = 20*log10(rms(real(filtSigB))) - ...
                                            20*log10(rms(real(filtSigA)));
    
            end
                                            
        %Compute linear approximation of ILD data (includes mirroring
        %of data to negative azimuth and setting middle to 0)
        [~, ind] = min(abs(kv.azis - kv.AziCrvfitRange));
        crvfitRange = (-kv.azis(ind) : 5 : kv.azis(ind));
        data1 = ILD(bin,1:ind); data1(1) = 0;
        data2 = fliplr(data1(2:end)); data = [-data2,data1];
        
        [a, ~] = polyfit(crvfitRange, data,1);
        IldAngSlope(bin) = a(1);
        
        %Bin Rate Level Data
        data = mappingData.calSpikeRatePerNeuronPerLevel;
        data = mean(data(:,kv.binPosInd(bin,:)),2);
        RateLevelPerBin(bin,:) = data; 
        
        %Compute linear approximation of rate level data
        [~, indLow] = min(abs(kv.dBRange - ...
                                kv.RateLevelCrvfitRange(1)));
        [~, indHigh] = min(abs(kv.dBRange - ...
                                kv.RateLevelCrvfitRange(2)));
        
        crvfitRange = (kv.dBRange(indLow : indHigh));                   
        [a, ~] = polyfit(crvfitRange, data(indLow : indHigh)',1);
        RateLevelSlope(bin) = a(1);
    end  
    
    mappingData.RateLevelPerBin = RateLevelPerBin;
    mappingData.ILDcentFreq = ILDcentFreq;
    mappingData.ILD = ILD;
    mappingData.IldAngSlope = IldAngSlope;
    mappingData.RateLevelSlope = RateLevelSlope;
    
end


%% Response Difference AN model
function mappingData = computeResponseDiffMapping(mappingData,kv)
    %Initialize variables    
    ResponseDiffSlope = zeros(1,numel(kv.binPos));
    SpkDiffPerAzi = zeros(numel(kv.binPos),numel(kv.azis)); 
    
    for bin = 1 : numel(kv.binPos)
        
        %Bin AN Response Diff Data
        data1 = mean(mappingData.calSpikeDiffPerNeuronPerAzi,3);
        data1 = mean(data1(:,kv.binPosInd(bin,:)),2);
        SpkDiffPerAzi(bin,:) = data1;
                                           
        %Compute linear approximation of response diff data 
        %(includes mirroring of data to negative azimuth and 
        %setting middle to 0)
        [~, ind] = min(abs(kv.azis - kv.AziCrvfitRange));
        crvfitRange = (-kv.azis(ind) : 5 : kv.azis(ind));
        data1 = data1(1 : ind); data1(1) = 0; 
        data2 = flipud(data1); data = [-data2;data1(2:end)];
 
        [a, ~] = polyfit(crvfitRange, data',1);
        ResponseDiffSlope(bin) = a(1);
        
    end
    mappingData.SpkDiffPerAzi = SpkDiffPerAzi;
    mappingData.ResponseDiffSlope= ResponseDiffSlope;
end

%% Probabilistic Model
function predictions = computeMaxLikelihood(mappingData, SpkDiffPerBin, kv)

predictions = zeros(1,size(SpkDiffPerBin,2));

for frame = 1 : size(SpkDiffPerBin,2)                            
    
    predictionsT = zeros(1,size(mappingData.Avg,1));
    
    for ang = 1 : size(mappingData.Avg,1)
                 
             MU = mappingData.Avg(ang,:);
             
             SIGMA = mappingData.covariance(ang,:);
     
             data = SpkDiffPerBin(:,frame);
             data = data';
             [prob] = mvnpdf(data,MU,SIGMA);
     
             predictionsT(ang) = prob;
   
     end 
            [~, indA] = max(predictionsT);
            predictions(frame) = kv.azis(indA); 
            
end
end



