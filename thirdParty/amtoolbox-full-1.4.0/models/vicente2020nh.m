function [predicted_SNR, BE, BU] = vicente2020nh(target_in,int_in,fs)
%VICENTE2020NH Compute the effective SNR taking into account better ear and bmld advantages
%   Usage: [predicted_SNR, BE, BU] = vicente2020nh(target_in,int_in,fs)
%
%   Input parameters:
%     target_in     : target
%     int_in        : interferer
%     fs            : sampling frequency [Hz]
%
%   Output parameters:
%     predicted_SNR : SNR predicted by the model
%     BE            : better-ear advantage
%     BU            : binaural masking level difference advantage
%
%   VICENTE2020NH computes the effective SNR taking into account BU and BE 
%   by respective time frames, taking the target and interferer signals 
%   (sampled at fs) as inputs
%
%   See also: lavandier2022 vicente2020nh vicente2020 prudhomme2020 leclere2015 jelfs2011 exp_lavandier2022
%
%   References:
%     M. Lavandier, T. Vicente, and L. Prud'homme. A series of snr-based
%     speech intelligibility models in the auditory modeling toolbox. Acta
%     Acustica, 2022.
%     
%     T. Vicente and M. Lavandier. Further validation of a binaural model
%     predicting speech intelligibility against envelope-modulated noises.
%     Hearing Research, 390(107937), 2020.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/vicente2020nh.php


%   #StatusDoc: Perfect
%   #StatusCode: Good
%   #Verification: Verified
%   #Requirements: MATLAB
%   #Author: Matthieu Lavandier
%   #Author: Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



%MODEL PARAMETERS
ceiling=20;   %maximum SNR allowed 
window_size_BE=round(0.024*fs);     %set the BE/monaural time frame size to 24 ms (effective duration 12 ms for Hann windows) 
window_size_BU=round(0.300*fs);     %set the BU/binaural time frame size to 300 ms (effective duration 150 ms for Hann windows) 

%compute target info on the long term signal
[targ_left, targ_right, target_phase, fc_target] = local_get_target_stats(target_in,fs);

ll=length(int_in);   %length of interferer determines the number of time frames

%COMPUTE BE SNR
hh_BE=hann(window_size_BE)*[1 1]; % create a binaural Hann window
weighted_better_ear=zeros(floor(ll/(window_size_BE/2))-1,1);
jj=1;
for kk=1:window_size_BE/2:ll-window_size_BE
    weighted_better_ear(jj) = local_betterear_timeframe(int_in(kk:kk+window_size_BE-1,:).*hh_BE,fs,ceiling,targ_left,targ_right,fc_target);
    jj=jj+1;
end
BE=mean(weighted_better_ear);

%COMPUTE BU
hh_BU=hann(window_size_BU)*[1 1]; % create an binaural Hann windows
weighted_bmld=zeros(floor(ll/(window_size_BU/2))-1,1);
jj=1;
for kk=1:window_size_BU/2:ll-window_size_BU
    weighted_bmld(jj) = local_bmld_timeframe(int_in(kk:kk+window_size_BU-1,:).*hh_BU,fs,target_phase,fc_target);
    jj=jj+1;
end
BU=mean(weighted_bmld);

%EFFECTIVE SNR INVOLVING BE AND BU
predicted_SNR=BU+BE;

end

function [left_spectrum, right_spectrum, interaural_phase, fc] = local_get_target_stats(sig,fs)
%Compute the (left and right) spectrum and interaural phase of the input signal sig (stereo files=2-colum matrix) sampled at fs for
%each (gammatone) frequency band with center frequency given by fc
%Computations are similar to those used in lavandier2022.m

nerbs = 1:0.5:round(f2erbrate(fs/2));
fc = zeros(size(nerbs));
interaural_phase = zeros(size(nerbs));
left_spectrum = zeros(size(nerbs));
right_spectrum = zeros(size(nerbs));

for n = 1:length(nerbs)
    % get filter center frequency
    fc(n) = round(erbrate2f(nerbs(n)));         
    % filter target
    sig_left = auditoryfilterbank(sig(:,1),fs,fc(n), 'lavandier2022');    
    sig_right = auditoryfilterbank(sig(:,2),fs,fc(n), 'lavandier2022');   
    [interaural_phase(n), ~] = local_do_xcorr(sig_left,sig_right,fs,fc(n)); % cross-correlate
    % spectrum in dB based on rms of the signals (independent of signal length but not of 0 padding) rms=10*Log10(mean(sig.*sig))
    left_spectrum(n) = 10*log10(mean(sig_left.^2));
    right_spectrum(n) = 10*log10(mean(sig_right.^2));     
end
end

function [weighted_better_ear] = local_betterear_timeframe(int_in,fs,ceiling,targ_left,targ_right,fc_target)
%compute BE (better ear) SNR in a time frame taking the interferer signal and
%target spectrum at the ears as inputs
%if there is no interferer energy in the time frame, BE is set to ceiling

nerbs = 1:0.5:round(f2erbrate(fs/2));
fc = zeros(size(nerbs));
if length(fc)~= length(fc_target)                %check that fc for target and masker are the same
    disp('Target and masker stats should be computed at the same frequency') 
end
better_ear_prediction = zeros(size(nerbs));

for n = 1:length(nerbs)
    % get filter cf
    fc(n) = round(erbrate2f(nerbs(n)));         
    % filter interferer
    int_left = auditoryfilterbank(int_in(:,1),fs,fc(n), 'lavandier2022');       
    int_right = auditoryfilterbank(int_in(:,2),fs,fc(n), 'lavandier2022');
    if sum(int_left.^2)==0 ||  sum(int_right.^2)==0         %if there is no interferer energy in the time frame, set BE to ceiling
        better_ear_prediction(n)=ceiling;
    else
    left_SNR = targ_left(n) - 10*log10(mean(int_left.^2));
    right_SNR = targ_right(n) - 10*log10(mean(int_right.^2));
    better_ear_prediction(n) = min(ceiling,max(left_SNR,right_SNR));
    end  
end

%integration accross frequency using SII weightings
weightings = f2siiweightings(fc);
weighted_better_ear = sum(better_ear_prediction.*weightings');

end

function [weighted_bmld] = local_bmld_timeframe(int_in,fs,target_phase,fc_target)
%compute BU advantage (bmld) in a time frame taking the interferer signal and target
%interaural phase as inputs

nerbs = 1:0.5:round(f2erbrate(fs/2));
fc = zeros(size(nerbs));
if length(fc)~= length(fc_target)                %check that fc for target and masker are the same
    amt_disp('Target and masker stats should be computed at the same frequencies') 
end
bmld_prediction = zeros(size(nerbs));

for n = 1:length(nerbs)
    % get filter cf
    fc(n) = round(erbrate2f(nerbs(n)));         
    % filter interferer
    int_left = auditoryfilterbank(int_in(:,1),fs,fc(n), 'lavandier2022');       
    int_right = auditoryfilterbank(int_in(:,2),fs,fc(n), 'lavandier2022');
    if sum(int_left.^2)==0 ||  sum(int_right.^2)==0         %if there is no interferer energy in the time frame, set bmld to 0 
        bmld_prediction(n)=0;
    else
    % BMLD
    [int_phase, int_coherence] = local_do_xcorr(int_left,int_right,fs,fc(n)); % cross-correlate
    bmld_prediction(n) = bmld(int_coherence,target_phase(n),int_phase,fc(n));
    end   
end

%integration accross frequency using SII weightings
weightings = f2siiweightings(fc);
weighted_bmld = sum(bmld_prediction.*weightings');

end

function [phase, coherence] = local_do_xcorr(left, right, fs, fc)
    [iacc, lags] = xcorr(left,right,round(fs/(fc*2)),'coeff'); %round(fs/(fc*2)) is for conformity with Durlach's 1972 formulation which allows time delays up to 
                                                               %+/- half the period of the channel centre frequency.
    [coherence, delay_samp] = max(iacc);
    phase = fc*2*pi*lags(delay_samp)/fs;
end


