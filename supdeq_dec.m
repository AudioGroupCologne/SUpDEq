%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function sparseHRTFdataset_dec = supdeq_dec(sparseHRTFdataset, eqDataset, targetDistance)
%
% This function performs a distance error compensation (dec) of a (sparse)
% HRTF set. It determines the onsets of the equalized datasets and 
% adjusts/shifts the (sparse) dataset baset on these onsets.
%
% Output:
% sparseHRTFdataset_dec - Sparse HRTF dataset with shifted initial delays
%
% Input:        
% sparseHRTFdataset     - Struct with sparse HRTF dataset
% eqDataset:            - Struct with equalization dataset
%                         (supdeq_getEqDataset)
% targetDistance:       - Correct distance to target in m
%                         Default: 2.00m
%
% Dependencies: AKTools
%
% Reference:
% C. Pörschmann and J. M. Arend, ?How positioning inaccuracies influence 
% the spatial upsampling of sparse head-related transfer function sets,? 
% in Proceedings of the International Conference on Spatial Audio - ICSA 2019
%   
% (C) 2019 by CP,  Christoph Pörschmann
%             JMA, Johannes M. Arend  
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function sparseHRTFdataset_dec = supdeq_dec(sparseHRTFdataset, eqDataset, targetDistance)

if nargin < 3 || isempty(targetDistance)
    targetDistance = 2;
end

%Get data from eqDataset
c = eqDataset.c;
r = eqDataset.radius;

%Get data from sparseHRTFdataset and perform equalization
fs = sparseHRTFdataset.f(end)*2;
sparseSamplingGrid = sparseHRTFdataset.samplingGrid;
Nsparse = sparseHRTFdataset.Nmax;
eqHRTFdataset = supdeq_eq(sparseHRTFdataset,eqDataset,Nsparse,sparseSamplingGrid,false);

% Determine onsets
for kk = 1:length(sparseHRTFdataset.samplingGrid(:,1)) 
    %Get HRIRs
    eqHRIR_L=ifft([eqHRTFdataset.HRTF_L(kk,:) zeros(1,length(eqHRTFdataset.HRTF_L(kk,:))-2)],'symmetric');
    eqHRIR_R=ifft([eqHRTFdataset.HRTF_R(kk,:) zeros(1,length(eqHRTFdataset.HRTF_R(kk,:))-2)],'symmetric');
    %Determine onsets with 10 times oversampling and write in sparseHRTFdataset struct
    sparseHRTFdataset.onset(kk,1) =  AKonsetDetect(eqHRIR_L',10,-1,'rel',[4000 fs]);
    sparseHRTFdataset.onset(kk,2) =  AKonsetDetect(eqHRIR_R',10,-1,'rel',[4000 fs]);
end
%Get mean onset
meanOnset = mean(mean(sparseHRTFdataset.onset));

%"Normalize" onsets by substracting mean and transform to seconds
sparseHRTFdataset.onset(:,1)=(sparseHRTFdataset.onset(:,1)-meanOnset)/fs;
sparseHRTFdataset.onset(:,2)=(sparseHRTFdataset.onset(:,2)-meanOnset)/fs;

for kk = 1:length(sparseHRTFdataset.samplingGrid(:,1)) 
    %Get "startingDistance" by adding (potential) offset/distance error to
    %target distance
    startingDistance(kk,:) = targetDistance+sparseHRTFdataset.onset(kk,:)*c;
    %Get targetDirection from respective sampling grid
    targetDirection = sparseSamplingGrid(kk,:);
    %Shift distance / compensate for distance errors
    %Use mean of left/right ear starting distance as starting distance for shift
    [sparseHRTFdataset.HRTF_L(kk,:),sparseHRTFdataset.HRTF_R(kk,:)] = supdeq_shiftDistance(sparseHRTFdataset.HRTF_L(kk,:),sparseHRTFdataset.HRTF_R(kk,:),mean(startingDistance(kk,:)),targetDistance,targetDirection,r,fs,c);
end

%Save in new variable and add starting/target distance info
sparseHRTFdataset_dec = sparseHRTFdataset;
sparseHRTFdataset_dec.decStartingDistance = startingDistance;
sparseHRTFdataset_dec.decTargetDistance = targetDistance;
fprintf('Performed distance error compensation to target distance of %d m\n',targetDistance);
    
end

