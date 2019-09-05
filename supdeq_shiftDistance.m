%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [shiftedHRTF_L,shiftedHRTF_R] = supdeq_shiftDistance(startingHRTF_L, startingHRTF_R, startingDistance, targetDistance, targetDirection, headRadius, fs, c)
%
% This function performs a shifts HRTFs in distance considering a spherical head model
%
% Output:
% shiftedHRTF_L / R     - HRTFs shifted / adapted in distance
%
% Input:        
% startingHRTF_L / R    - Input HRTFs to be shifted / adapted in distance
% startingDistance:     - Distance at which the datasets were measured
% targetDistance:       - Distance to which the datasets are shifted / adapted
% targetDirection:      - Direction of the target (spatial sampling grid, 
%                         Qx2 matrix with azimuth and elevation) in degree.
% headRadius:           - Radius of the spherical head model. For details
%                         refer to AKsphericalHead.
% fs                    - Sampling rate in Hz
%                         Default: 48000
% c                     - Speed of sound in m/s
%                         Default: 343
%
% Dependencies: AKTools
%   
% (C) 2019 by CP,  Christoph Pörschmann
%             JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [shiftedHRTF_L,shiftedHRTF_R] = supdeq_shiftDistance(startingHRTF_L, startingHRTF_R, startingDistance, targetDistance, targetDirection, headRadius, fs, c)

if nargin < 7 || isempty(fs)
    fs = 48000;
end

if nargin < 8 || isempty(c)
    c = 343;
end

%Generate sampling grid with distance information based on targetDirection
%and startingDistance/targetDistance
sg=[targetDirection(1),90-targetDirection(2),startingDistance;targetDirection(1),90-targetDirection(2),targetDistance];
%Define symmetrical ear position of spherical head model
ear=[90 0];

%Get respective spherical head impulse responses for starting and target distance
hDistances=AKsphericalHead(sg,ear,false,headRadius,sg(1,3),100,length(startingHRTF_L)*2-2,fs,c);

%Get distance variation functions for left and right ear
H_L=fft(squeeze(hDistances(:,2,1)))./fft(squeeze(hDistances(:,1,1)));
H_R=fft(squeeze(hDistances(:,2,2)))./fft(squeeze(hDistances(:,1,2)));

%Get rid of mirror spectrum
H_L=H_L(1:length(startingHRTF_L)).';
H_R=H_R(1:length(startingHRTF_R)).';

%Apply distance variation functions to HRTFs
shiftedHRTF_L=H_L.*startingHRTF_L;
shiftedHRTF_R=H_R.*startingHRTF_R;
        
end

