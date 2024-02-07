%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function output = supdeq_listen(HRIRs_sfd,testSignal,AzEl,playback)
%
% This function convolves any arbitrary test signal with an HRIR for the 
% given azimuth (az) and elevation (el). The output will be played back
% automatically by default.
%
% Output:
% output        - Convolution result of test signal and HRIR in time domain 
%
% Input:        
% HRIRs_sfd     - Struct with Spherical Harmonic coefficients 
%                 (SH-Coefficients) for the left (Hl_nm) and right (Hr_nm) 
%                 channel/ear, absolute frequency scale f, 
%                 transform order N, and FFToversize
% testSignal    - Mono test signal
% AzEl          - Spatial sampling point (Q x 2 matrix), where the value 
%                 is the azimuth and the second the
%                 elevation (both in degree).
% playback      - Play the result with soundsc - true/false
%                 Default: true
% fs            - Sampling rate
%                 Default: 48000
%
% Dependencies: -
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

%%
function output = supdeq_listen(HRIRs_sfd,testSignal,AzEl,playback)

if nargin < 4 || isempty(playback)
    playback = true;
end

if nargin < 5
    fs = 48000;
end

%Get HRIR according to AzEl
[hrir(:,1),hrir(:,2)] = supdeq_getArbHRIR(HRIRs_sfd,AzEl,'DEG');

%Convolve HRIR with test signal (FFT convolution)
NFFT = length(hrir) + length(testSignal) - 1;
testSignal = [testSignal,testSignal];
output = ifft(fft(hrir,NFFT) .* fft(testSignal,NFFT));

%Play back if intended...
if playback
    soundsc(output,fs);
end
