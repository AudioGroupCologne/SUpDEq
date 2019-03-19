function plot_roenne2012(stim_level,waveVamp,waveVlat, simpot, ANout,varargin)
%plot_roenne2012 Plot the output from the Roenne 2012 model
%   Usage: plot_roenne2012(waveVamp,waveVlat,...);
%
%   Input parameters:
%     waveVamp   : Amplitude of simulated ABR wave V.
%     waveVlat   : Latency of simulated ABR wave V peak.
%
%   PLOT_ROENNE2012(stim_level,waveVamp,waveVlat) plots the output from
%   ROENNE2012.
%
%   The flag may be one of:
%
%     'fsmod',fsmod     Auditory nerve model frequency.
%                       Default value is 200000.
%
%     'flow',flow       Auditory nerve model lowest center frequency.
%                       Default value is 100 Hz.
%
%     'fhigh',fhigh     Auditory nerve model highest center frequency.
%                       Default value is 16000 Hz.
%      
%     'min_modellength',mn 
%                       Minimum length of modelling measured in ms.
%                       Default value is 40.
%
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%   References:
%     C. Elberling, J. Calloe, and M. Don. Evaluating auditory brainstem
%     responses to different chirp stimuli at three levels of stimulation. J.
%     Acoust. Soc. Am., 128(1):215-223, 2010.
%     
%     F. M. RÃ¸nne, T. Dau, J. Harte, and C. Elberling. Modeling auditory
%     evoked brainstem responses to transient stimuli. The Journal of the
%     Acoustical Society of America, 131(5):3903-3913, 2012. [1]http ]
%     
%     M. S. A. Zilany and I. C. Bruce. Representation of the vowel (epsilon)
%     in normal and impaired auditory nerve fibers: Model predictions of
%     responses in cats. J. Acoust. Soc. Am., 122(1):402-417, jul 2007.
%     
%     References
%     
%     1. http://scitation.aip.org/content/asa/journal/jasa/131/5/10.1121/1.3699171
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/plot/plot_roenne2012.php

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

% Define input flags
definput.keyvals.fsmod=200000;
definput.keyvals.flow = 100;
definput.keyvals.fhigh = 16000;
[flags,kv]      = ltfatarghelper({},definput,varargin);

%% PLOTS, extra plots created for all conditions used, i.e. three
% plots for each stimulus level x each chirp sweeping rate. If this
% is switched on and all other variables are set to default,
% 3 (levels) x 6 (chirps / clicks) x 3 (different figures) = 48
% figures will be created.

% load Unitary Response and its samlping frequency
[ur,fs]=data_roenne2012;

% Plot simulated ABR
figure;
t = 1e3.*[0:length(simpot)-1]./fs;
plot(t,simpot,'k','linewidth',2)
xlabel('time [ms]'), title(['Simulated ABR at ' num2str(stim_level) 'dB']),
set(gca, 'fontsize',12)

% create frequency vector indentical to the CFs simulated
% location of lowest and highest centre frequency 
xlo     = (1.0/0.06)*log10((kv.flow/165.4)+0.88);
xhi     = (1.0/0.06)*log10((kv.fhigh/165.4)+0.88);	

% equal spaced distances on the BM
vX      = linspace(xlo,xhi,500);

% and the resulting frequency vector
vFreq   = 165.4*(10.^(0.06*vX)-0.88); 


% Plot "AN-gram" - spectrogram-like representation of the discharge
% rate after the IHC-AN synapse 
figure;
set(gca, 'fontsize',12);
imagesc(ANout');
title(['ANgram at ' num2str(stim_level) 'dB'])
set(gca,'YTick',[1 100 200 300 400 500]),
set(gca,'YTicklabel',round(vFreq([1 100 200 300 400 500]))), 
set(gca,'XTick',(0:1000:8000)),
set(gca,'XTicklabel',(0:1000:8000)/kv.fsmod*1000);
ylabel('model CF'), 
xlabel('time [ms]'),
colorbar;

% Plot "AN-UR-gram" - spectrogram-like representation of the discharge
% rate convolved line by line with the UR. 
figure
ANUR = resample(ANout,fs,kv.fsmod);
ANUR = filter(ur,1,ANUR);
imagesc(ANUR');
set(gca,'YTick',[1 100 200 300 400 500]);
set(gca,'YTicklabel',round(vFreq([1 100 200 300 400 500])));
set(gca,'XTick',(0:150:1500));
set(gca,'XTicklabel',(0:150:1500)/fs*1000);
ylabel('model CF');
xlabel('time [ms]');
colorbar;
title(['AN-URgram at ' num2str(stim_level) 'dB']);



