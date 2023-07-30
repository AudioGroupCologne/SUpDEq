%DEMO_MAY2011 Demo of the model estimating the azimuths of concurrent speakers
%
%
%   DEMO_MAY2011 generates figures showing the result of the model estimating
%   the azimuth position of three concurrent speakers. Also, it returns the
%   estimated azimuths.
% 
%   Set demo to the following flags to shows other conditions:
%   
%   1R
%      one speaker in reverberant room
%
%   2
%      two speakers in free field
%
%   3
%      three speakers in free field (default)
%
%   5
%      five speakers in free field
%
%   Figure 1: Time-frequency-based azimuth estimates
%
%      This figure shows the azimuth estimates in the time-frequency
%      domain for three speakers.
%
%   Figure 2: Interaural time differences (ITDs)
%
%      This figure shows the ITDs in the time-frequency domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   Figure 3: Interaural level differences (ILDs)
%
%      This figure shows the ILDs in the time-frequency domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   Figure 4: Interaural coherence
%
%      This figure shows the interaural coherence in the time-frequency domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   Figure 5: Frame-based azimuth estimates
%
%      This figure shows the azimuth directions in the time domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   Figure 6: GMM pattern
%
%      This figure shows the pattern and the histogram obtained from the
%      GMM-estimator for the mixed signal of three concurrent speakers.
%
%   See also: may2011
%
%   References:
%     T. May, S. van de Par, and A. Kohlrausch. A probabilistic model for
%     robust localization based on a binaural auditory front-end. IEEE Trans
%     Audio Speech Lang Proc, 19:1--13, 2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_may2011.php


%   #Author: Tobias May (2009)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Select binaural recordings
% 
% 
% Select one demo
if ~exist('demo','var')
  demo='3';
end
 
% Create signals
switch lower(demo)
    case '1r'
        % Create input signal
        [signal,fs] = sig_competingtalkers('one_speaker_reverb');
        
        % Find all active sources
        nSources = 1;
        
    case '2'
        % Create input signal
        [signal,fs] = sig_competingtalkers('two_of_three');
        
        % Find all active sources
        nSources = 2;
        
    case '3'
        % Create input signal
        [signal,fs] = sig_competingtalkers('three_of_three');
        
        % Find all active sources
        nSources = 3;
        
    case '5'
        % Create input signal
        [signal,fs] = sig_competingtalkers('five_speakers');
        
        % Find all active sources
        nSources = 5;
end


%% Perform GMM-based sound source localization
% 
% 
% Perform localization
out = may2011(signal,fs);


%% Plot results
% 
% 
% Plot time-frequency-based azimuth estimates
figure;
imagesc(out.azimuth,[-90 90]);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Time-frequency-based azimuth estimates')
colorbar;
may2011_cbarlabel('Azimuth (deg)')
axis xy;

% Plot binaural cues
figure;
imagesc(out.itd,[-1e-3 1e-3]);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural time difference (ITD)')
colorbar;
may2011_cbarlabel('ITD (ms)')
axis xy;

figure;
imagesc(out.ild);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural level difference (ILD)')
colorbar;
may2011_cbarlabel('ILD (dB)')
axis xy;

figure;
imagesc(out.ic);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural coherence (IC)')
colorbar;
may2011_cbarlabel('IC')
axis xy;

% Plot frame-based azimuth estimates
%figure;
%plot(out.azFrames,'k.','linewidth',2);
%xlabel('Number of frames')
%ylabel('Azimuth (deg)')
%title('Frame-based azimuth estimates')
%xlim([-inf inf])
%ylim([-90 90])
%grid on;
%axis xy;

% Histogram analysis of frame-based localization estimates
azEst=may2011_estazimuthgmm(out,'HIST',nSources,1)

