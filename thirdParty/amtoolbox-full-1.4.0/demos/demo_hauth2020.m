%DEMO_HAUTH2020 framework to test the binaural speech intelligibility model hauth2020
%
%   hauth2020 predicts the SRTs from the Anechoic Experiment of Beutelmann and Brand 
%   (2006). In this experiment, SRTs are simulated for speech located at 0 deg in
%   the horizontal plane, and noise located at different positions [angles].
%   The binaural processing in this model works blindly and only requires the
%   mixture of speech an noise [required signal]. The back-end used here is the SII. 
%   SII values are produced for different signal-to-noise ratios. Depending
%   on the back-end you want to use, optional signals have to be used. They
%   are processed in the same way as the required signal. In this Demo, the
%   optional signals are clean speech and noise. 
%
%   Figure 1: Predicted speech intelligibility index as a function of SNR
%
%   See also: hauth2020
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_hauth2020.php


%   #Author: Christopher Hauth
%   #Author: Dr. Thomas Brand

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%-------------------------------------------------------------------------%
fast = 1; % set to 1 to get faster results, set to 0 to get more accurate results
display_level = 'no_debug'; % set to 'debug' to see more information, set to 'no_debug' to have less mess on your display
list_clean(1).name = '0_05413.wav';
list_clean(2).name = '0_31601.wav';
list_clean(3).name = '0_46393.wav';
%---------------Experimental Conditions-----------------------------------%
% Define your experimental conditions, number of Monte Carlo simulations
% (binaural processing) and number of sentences (statistics across sentences)
% For Matrix type sentences it is recommended to use 10 sentences, where each word
% of the test appears once.
if fast
  amt_disp('Fast mode activated - results may be inaccurate.'); 
  vSNR_test = -20:5:0;
  iNumofSentences = 1;
  iNumofMonteCarlo = 3;
else
  amt_disp('Accurate mode activated - please be patient.'); 
  vSNR_test = -20:0;
  iNumofSentences = 10;
  iNumofMonteCarlo = 10;
end

fs_model = 44100;

sentences_clean = {list_clean.name};
sentence_choose_clean = sentences_clean;
randomizer = randperm(length(sentence_choose_clean));

vangles_test  = 45;
% Vector of input SNRs. For each SNR, (iNumofSentences x iNumofMonteCarlo) SII values are obtained.
% Make sure to test different different SNRs in order to be able to map the
% SII to an SRT, e.g.:

% If you want to know the SII of a single SNR, use only one value:
 %vSNR_test = -18;                                
%% Calibration
% the calibration factor is mean level between the ears
% The calibration can be adjusted for your needs. For the SII, which is used here, 65
% dB FS (relative to full scale) is assumed to be 65 dB SPL.
% However, if you only aim for the
% resynthesized output, please adjust this level to avoid clipping. 
lev_desired = 65;   % the value 65 is required here for correct use of the SII;

% Use the co-located noise condition to calibrate the input (For the speeech signal, also the noise is used)
[calibnoise, fs] = amt_load('hauth2020', '0_speaker_reference_olnoise.wav');
calibnoise = calibnoise(1:end-round(1.5*fs),:);
    
lev_Speech = 20*log10(rms(calibnoise));  % actual rms-level of speech
lev_S = mean(lev_Speech);                % the mean between both ears is considered
Delta_L_speech = lev_desired - lev_S;    % Calibration Gain is the difference between the desired level and the actual level
Delta_L_speech_lin = 10.^((Delta_L_speech)/20); % Convert to linear gain

% Similar calibration of the noise
lev_Noise = 20*log10(rms(calibnoise));  % actual rms-level of the noise
lev_Noise = mean(lev_Noise);                % reference is MEAN level between the two ears
Delta_L_noise = lev_desired - lev_Noise;    % Calibration Gain is the difference between the desired level and the actual level
Delta_L_noise_lin = 10.^((Delta_L_noise)/20);% Convert to linear gain
%-------------------------------------------------------------------------%
clear sii_min_all sii_max_all sii_syn_all sii_L_all sii_R_all
    % Iterate through all SNRs
    for kk = 1:length(vSNR_test)
      amt_disp(['Processing SNR ' num2str(kk) ' out of ' num2str(length(vSNR_test)) '.']);
        % Iterate through the different sentences
        for ll = 1:iNumofSentences
            amt_disp(['Processing sentence ' num2str(ll) ' out of ' num2str(iNumofSentences) '.']);
            sentence_clean = sentence_choose_clean{randomizer(ll)};         
            % Read sentences and noise from wav files
            [speech_clean, fs_s]= amt_load('hauth2020', sentence_clean);
            [noise, fs_n] = amt_load('hauth2020', sprintf('%d_speaker_reference_olnoise.wav', vangles_test));
            % resample signals if necessary 
            if fs ~= fs_s
                speech_clean = resample(speech_clean,fs_model,fs_s);
                noise = resample(noise,fs_model,fs_n);
            end
            % Get length of the speech signal
            lenSpeech = length(speech_clean);
            lenNoise = length(noise);
            
            % Truncate noise to have the same lenght as speech
            noise = noise(1:lenSpeech,:);
            speech_clean = speech_clean(1:lenSpeech,:);
            % adjust level of speech:
            speech_clean = Delta_L_speech_lin.*speech_clean;
            
            % adjust level of noise:
            noise = Delta_L_noise_lin.*noise;
            
            % adjust SNR of mixed input signal (speech + noise)
            % This is a required signal:
            mixed_input = 10.^((vSNR_test(kk))/20).*speech_clean+noise;
            inputLen = length(mixed_input);
            
            % Adjust level of the clean speech if you want to 
            % use it as an optional input:
            speech_clean_proc = 10.^((vSNR_test(kk))/20).*speech_clean;
            
            % All optional signals are arranged in a matrix.
            % Here: [S_l(:) S_r(:) N_l(:) N_r(:)]
            OptionalSignals = [speech_clean_proc noise];
            
            % Apply the binaural model to the mixed signal 
            % (and optionally to the clean speech and noise)
            % Monte Carlo simulations are used to model the binaural uncerntainty
            for oo=1:iNumofMonteCarlo
                % Do binaural processing:
                % out_struct contains the processed mixed signal as well as
                % the processed optional signals: Moreover, as the SII is
                % used as back-end, it also contains the frequency-specific 
                % levels
                amt_disp(['Processing Monte Carlo ' num2str(oo) ' out of ' num2str(iNumofMonteCarlo) '...'],'volatile');
                 out_struct = hauth2020(mixed_input, fs,'OptSigs',OptionalSignals,display_level);

                % Speech Intelligibility back-end: (SII in this example)
                % Use your speech intelligibility back-end here
                [sii_min_temp(oo),A,Z] = hauth2020_sii(out_struct.levels.LevelOptSig1min,out_struct.levels.LevelOptSig2min,-Inf*ones(30,1),2,0); 
                [sii_max_temp(oo),A,Z] = hauth2020_sii(out_struct.levels.LevelOptSig1max,out_struct.levels.LevelOptSig2max,-Inf*ones(30,1),2,0); 
                [sii_syn_temp(oo),A,Z] = hauth2020_sii(out_struct.levels.LevelOptSig1syn,out_struct.levels.LevelOptSig2syn,-Inf*ones(30,1),2,0); 
                [sii_L_temp(oo),A,Z]   = hauth2020_sii(out_struct.levels.LevelOptSig1L,out_struct.levels.LevelOptSig2L,-Inf*ones(30,1),2,0); 
                [sii_R_temp(oo),A,Z]   = hauth2020_sii(out_struct.levels.LevelOptSig1R,out_struct.levels.LevelOptSig2R,-Inf*ones(30,1),2,0); 
                          
                sii_min_all(ll,kk,oo) = sii_min_temp(oo);
                sii_max_all(ll,kk,oo) = sii_max_temp(oo);
                sii_syn_all(ll,kk,oo) = sii_syn_temp(oo);
                sii_L_all(ll,kk,oo)   = sii_L_temp(oo);
                sii_R_all(ll,kk,oo)   = sii_R_temp(oo);     
            end
        end
    end
    amt_disp();
%% Plotting

sii_min_all_squeezed = squeeze(mean(sii_min_all,1));
sii_max_all_squeezed = squeeze(mean(sii_max_all,1));

sii_syn_all_squeezed = squeeze(mean(sii_syn_all,1));
sii_L_all_squeezed = squeeze(mean(sii_L_all,1));
sii_R_all_squeezed = squeeze(mean(sii_R_all,1));

for ii = 1:numel(vSNR_test)
  xlab{ii} = vSNR_test(ii);
end


subplot(3,1,1)
plot(max(sii_min_all_squeezed.' ),'k', 'linewidth', 2)
hold on
plot(min(sii_min_all_squeezed.' ),'k', 'linewidth', 2)
plot(sii_min_all_squeezed)
xlim([1 numel(vSNR_test)])
grid on
xlabel('SNR [dB]')
ylabel('Speech Intelligibility Index')
set(gca,'xtick',[1:numel(vSNR_test)],'xticklabels',xlab);
title('Minimum Speech Intelligibility Index')

subplot(3,1,2)
plot(max(sii_max_all_squeezed.' ),'k', 'linewidth', 2)
hold on
plot(min(sii_max_all_squeezed.' ),'k', 'linewidth', 2)
plot(sii_max_all_squeezed)
xlim([1 numel(vSNR_test)])
grid on
xlabel('SNR [dB]')
ylabel('Speech Intelligibility Index')
set(gca,'xtick',[1:numel(vSNR_test)],'xticklabels',xlab);
title('Maximum Speech Intelligibility Index')

subplot(3,1,3)
plot(max(sii_L_all_squeezed.' ), 'linewidth', 2)
hold on
plot(min(sii_R_all_squeezed.' ),'r', 'linewidth', 2)
xlim([1 numel(vSNR_test)])
grid on
xlabel('SNR [dB]')
ylabel('Speech Intelligibility Index')
set(gca,'xtick',[1:numel(vSNR_test)],'xticklabels',xlab);
legend('Left ear SII', 'Right ear SII','location', 'southeast')
title('L/R Speech Intelligibility Index')


