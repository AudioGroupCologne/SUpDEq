%DEMO_JOERGENSEN2013 Demo for the multi-resolution speech-based envelope spectrum model 
%
%   DEMO_JOERGENSEN2013 computes the signal-to-noise envelope-power ratio
%   (SNRenv) for a speech sentence in noise at an SNR of -3 dB.
%   The SNRenv can be then used to predict speech intelligibility. 
%
%   Pcorrect indicates the probability of correctly understanding the
%   sentence based on results from Joergensen et al. (2013). 
%
%   See also: joergensen2013 joergensen2011
%
%   References:
%     S. Joergensen and T. Dau. Predicting speech intelligibility based on
%     the signal-to-noise envelope power ratio after modulation-frequency
%     selective processing. J. Acoust. Soc. Am., 130(3):1475--1487, 2011.
%     
%     S. Jørgensen, S. D. Ewert, and T. Dau. A multi-resolution envelope
%     power based model for speech intelligibility. J. Acoust. Soc. Am.,
%     134(1):436--446, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_joergensen2013.php


%   #Author : Piotr Majdak

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

  % Load speech (from the CLUE corpus)
x=amt_load('joergensen2011','Danish_CLUE_10sentence_samples_22kHz.mat');

speech  = x.sentenceArray{1};
sentenceFileLevel = -26.00; % The RMS level of all CLUE sentence files corresponds to ...
SPL = 65; % ... this sound pressure level
SNR = -3;
fs = 22050;
speech = speech*10^((SPL-sentenceFileLevel)/20);
N = length(speech);

  % Load speech-shaped noise from the CLUE corpus
noise = amt_load('joergensen2011','SSN_CLUE_22kHz.wav');

% pick a random segment from the noise file
Nsegments = floor(length(noise)/N);
startIdx = randi(Nsegments-2 ,1)*N;
noise = noise(startIdx:startIdx+N -1)';
noise = noise./rms(noise)*10^((SPL-SNR)/20);
if size(noise) ~= size(speech), noise = noise'; end

  % Create a mixed signal
test = noise + speech;

  % Run the model without any priors
tmp = joergensen2013(test,noise,fs);
SNRenvs_noIOparameters = tmp.SNRenv

  % Run the model with parameters for the CLUE material from Joergensen et al., (2013).
IOparameters = [0.61 0.5 8000 0.6]; 
tmp = joergensen2013(test,noise,fs, IOparameters);
SNRenvs_withIOparameters = tmp.SNRenv
Pcorrect = tmp.P_correct


