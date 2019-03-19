%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [Hl_nm, Hr_nm] = supdeq_getEqDataset(earDistance, NFFT, fs)
%
% This function returns the sound pressure distribution of an incident 
% plane wave for different sound incidence directions on the left ear 
% position and on the right ear position of a sphere which models the 
% human head.
%
% Output:
% eqDataset     - Struct with SH-coefficients describing the sound 
%                 incidence at the letf/right ear position on the sphere 
%                 (Hl_nm/Hr_nm). Length of the SH-coefficients is NFFT/2+1 
%                 (single sided spectrum).         
%
% Input:        
% earDistance   - Distance between both ear positions in m 
%                 Default: 0.165
% NFFT          - FFT size of the SH-coefficients to be returned
%                 Default: 512
% fs            - Sampling rate
%                 Default: 48000
%
% Dependencies: SOFiA toolbox
%
% References:
% Benjamin Bernschütz: Microphone Arrays and Sound Field Decomposition 
% for Dynamic Binaural Recording. Ph.D. dissertation, Technical University
% Berlin (2016).
%   
% (C) 2018 by CP,  Christoph Pörschmann
%             JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function eqDataset = supdeq_getEqDataset(earDistance, NFFT, fs)

if nargin < 1
    earDistance = 0.165;
end

if nargin < 2
    NFFT = 512;
end

if nargin < 3
    fs = 48000;
end
    
%% Get SH-Coefficients with sofia_wgc

%Define required parameters
N       = 35; %Transform order N  
radius  = earDistance/2;
ac      = 2; %Array configuration, 2 - Rigid sphere with pressure transducers
delay   = 0; %Time delay in s - Set to 0 here

%Calculate SH-Coefficients with defaults c = 343 and wavetype = 0
eqDataset.Hl_nm = sofia_wgc(N, radius, ac, fs, NFFT, pi/2, pi/2, delay);
eqDataset.Hr_nm = sofia_wgc(N, radius, ac, fs, NFFT, -pi/2, pi/2, delay);
eqDataset.f = linspace(0,fs/2,NFFT/2+1);

end

