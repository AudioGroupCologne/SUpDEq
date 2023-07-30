function [twoears_benefit, weighted_bmld, weighted_better_ear] = leclere2015(target_in,int_in,fs)
%LECLERE2015 Compute the binaural useful-to-detrimental ratio for a reverberated target
%   Usage: [twoears_benefit, weighted_bmld, weighted_better_ear] = leclere2015(target_in,int_in,fs)
%
%   Input parameters:
%     target_in     : target
%     int_in        : interferer
%     fs            : sampling frequency [Hz]
%
%   Output parameters:
%     twoears_benefit     : useful-to-detrimental ratio
%     weighted_bmld       : weighted binaural masking level difference
%     weighted_better_ear : weighted better ear advantage
%
%   LECLERE2015 computed the binaural useful-to-detrimental ratio for a reverberated 
%   target and multiple stationary noise interferers. The early and late 
%   parts of the target BRIR are separated, the early part constitutes the useful
%   signal while the late part is concatenated with the interferer BRIRs to 
%   constitute the detrimental signal.
%
%
%   See also: lavandier2022 vicente2020nh vicente2020 prudhomme2020 leclere2015 exp_lavandier2022
%   jelfs2011
%
%   References:
%     M. Lavandier, T. Vicente, and L. Prud'homme. A series of snr-based
%     speech intelligibility models in the auditory modeling toolbox. Acta
%     Acustica, 2022.
%     
%     M. Lavandier and J. Culling. Speech segregation in rooms: Monaural,
%     binaural and interacting effects of reverberation on target and
%     interferer. J. Acoust. Soc. Am., 123(4):2237--2248, 2008.
%     
%     T. Leclère, M. Lavandier, and J. Culling. Speech intelligibility
%     prediction in reverberation: Towards an integrated model of speech
%     transmission, spatial unmasking and binaural de-reverberation. J.
%     Acoust. Soc. Am., 137(6):3335--3345, 2015.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/leclere2015.php


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


%%
%%%%%-----------------Windowing Settings------------------%%%%%%%%%%%%%%%
%default values of the room-independent version of the model in Lecl�re et al. (2015):
    %earlyLateLimit=30 ms, decayDuration=25 ms, windowShape=�linear�, lateWindow=�opposite�
    earlyLateLimit = 30;
    decayDuration=25;
    windowShape= 'linear'; % 'gate', 'linear', 'exp', 'cumulative'
    lateWindow='opposite';

padding = zeros(1024,2);
%early/late separation
[early, late] = local_windowing(target_in, earlyLateLimit, decayDuration, windowShape, lateWindow, fs);
useful = [early; padding];
detrimental = [late; padding; int_in];

%application of the binaural model on the useful and detrimental signals
nerbs = 1:0.5:round(f2erbrate(fs/2));
fc = zeros(size(nerbs));
bmld_prediction = zeros(size(nerbs));
better_ear_prediction = zeros(size(nerbs));

for n = 1:length(nerbs)
    % get filter cf
    fc(n) = round(erbrate2f(nerbs(n)));         
    % filter target and interferer separately
    targ_left = auditoryfilterbank(useful(:,1),fs,fc(n), 'lavandier2022');    
    targ_right = auditoryfilterbank(useful(:,2),fs,fc(n), 'lavandier2022');   
    int_left = auditoryfilterbank(detrimental(:,1),fs,fc(n), 'lavandier2022');       
    int_right = auditoryfilterbank(detrimental(:,2),fs,fc(n), 'lavandier2022');
    % BMLD
    [int_phase, int_coherence] = local_do_xcorr(int_left,int_right,fs,fc(n)); % cross-correlate
    [target_phase] = local_do_xcorr(targ_left,targ_right,fs,fc(n));    
    bmld_prediction(n) = bmld(int_coherence,target_phase,int_phase,fc(n));
    % better-ear SNR in dB based on energy (independent of 0 padding of
    % BRIRs) energy=10*Log10(sum(sig.*sig))
    left_SNR = 10*log10(sum(targ_left.^2)/sum(int_left.^2));
    right_SNR = 10*log10(sum(targ_right.^2)/sum(int_right.^2));
    better_ear_prediction(n) = max(left_SNR,right_SNR);  
end

%integration accross frequency using SII weightings
weightings = f2siiweightings(fc);
weighted_bmld = sum(bmld_prediction.*weightings');
weighted_better_ear = sum(better_ear_prediction.*weightings');

twoears_benefit = weighted_better_ear + weighted_bmld;

end


function [early, late] = local_windowing(signal, earlyLateLimit, decayDuration,  windowShape, lateWindow, fs)
% Time windowing function to split the impulse response into an early
% and late parts.
% signal : impulse response vector/matrix (can be multichannel). Each column
% represents an impulse response
% earlyLateLimit (in ms): time for which the window stands to 1
% decayDuration (in ms): duration of the transition for the window to decrease from 1 to 0
% shape of the window : the value 1 is hold until the earlyLateLimit,
% after it's going down to 0 according to the specified window shape :
%     - 'gate' :  from 1 to 0 in only one sample
%     - 'linear' : linear decrease
%     - 'exp' : exponential decrease
%     - 'cumulative' : decreasing according to the cumulative function of a 
%       normal distribution, the slope is settled with decayDuration
% lateWindow : defines the shape of the late window:
%     - 'opposite': the late window is the complement of the early
%     window, ie, earlyWindow + lateWindow = 1 (for every t)
%     - 'gate': the late window is set to 1 right after the early
%     window reached 0.
% fs: sampling rate (in Hz)
% 
% Details can be found in Lecl�re et al. (2015), J. Acoust. Soc. Am.,
% Vol 135, p 3335-3345


plotFigure = false;

[lines, columns] = size(signal);
t = (0:lines-1)';         
early = zeros(lines, columns);
late = zeros(lines, columns);

% Determine direct sound
dirSound = zeros(1,columns);
for m = 1:columns
    dirSound(m) = local_get_direct_lag(signal(:,m));
end
dirSound = min(dirSound);

window = zeros(lines,1);
if earlyLateLimit <0 || fs*earlyLateLimit/1000 > lines
        disp 'The Early/Late limit must be a positive value inferior to the length of the signal'
else
        % Convert temporal parameters in samples
        earlyLateLimit = round(fs * earlyLateLimit / 1000);
        decayDuration = round(fs * decayDuration / 1000);

        % Determine t1 and t2
        t1 = dirSound + earlyLateLimit;
        t2 = t1 + decayDuration;

        switch windowShape
            case 'gate'
                window = window;
                
            case 'cumulative'             
                
                sigma = (t1-t2) / (sqrt(2)*(erfinv(0.998)-erfinv(-0.998)));
                mu = t1-sigma*sqrt(2)*erfinv(0.998);
                
                window = 0.5*(1+erf((t-mu)/(sigma*sqrt(2))));
                
            case 'linear'                
                % linear equation knowing that the window is 1 at t = t1                
                window(:,1) = (t1 -t)/decayDuration + 1;

                
            case 'exp'
                
                window = [window(1:t1); exp(-t/slope)];
                window = window(1:lines);
                
            otherwise
                error('Wrong window shape')
        end
        
        % Impose 1 and 0 before t1 and after t2, respectively
        window(1:t1,1) = 1;
        window(t2:end, 1) = 0;
        
        switch lateWindow
            case 'opposite'
                lateWindowing = 1-window;
            case 'gate'
                lateWindowing = [zeros(t2,1); ones(length(signal) - t2,1)];
            otherwise
                error('Wrong late gate chosen')
        end
        
        for k = 1:columns % For each signal
        early(:,k) = signal(:,k).*window;
        late(:,k) = signal(:,k).*lateWindowing;
        end
end

if plotFigure
    figure, hold on
    plot(t/fs, early)
    plot(t/fs,window, 'LineWidth', 2)
end

end

function [direct] = local_get_direct_lag(ir)

%to be used with leclere2015.m
[~, column] = size(ir);
for k = 1:column
    [max_val, delay] = max(abs(ir(:,k)));
    if (max(abs(ir(1:delay-1,k))) > 0.8*max_val)
        delay(k) = local_get_direct_lag(ir(1:delay-1,k));
    end
end
direct = min(delay);

end

function [phase, coherence] = local_do_xcorr(left, right, fs, fc)
    [iacc, lags] = xcorr(left,right,round(fs/(fc*2)),'coeff'); %round(fs/(fc*2)) is for conformity with Durlach's 1972 formulation which allows time delays up to 
                                                               %+/- half the period of the channel centre frequency.
    [coherence, delay_samp] = max(iacc);
    phase = fc*2*pi*lags(delay_samp)/fs;
end


