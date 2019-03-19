function [waveVamp, waveVlat]  = roenne2012_chirp(levels,chirps,varargin)
%roenne2012_chirp Simulate chirp evoked ABRs
%   Usage: [waveVamp, waveVlat]  = ronne2012_chirp(flag)
%
%   Output parameters:
%     waveVamp   : Amplitude of simulated ABR wave V.
%     waveVlat   : Latency of simulated ABR wave V peak.
%
%   ronne2012_chirp(levels,chirps) returns simulations from Rønne et
%   al. (2012) for a range of input levels levels. The values of
%   chirps selects which stimuli from Elberling et al. (2010) stimuli to
%   simulate, 1 = click, 2 to 6 = chirp 1 to 5.
%
%   Simulates ABR responses to five chirps and one click using the ABR model
%   of Rønne et al. (2012). Fig. 6 or 7 of Rønne et al. (2012) can be
%   reproduced based on these data - use plot_ROENNE2012_CHIRP.  Simulations
%   can be compared to data from Elberling et al. (2010). Stimuli are
%   defined similar to Elberling et al. (2010).
%
%   This function takes the following optional parameters:
%
%     'plot2'   Plot extra figures for all individual simulations (3
%               figures x stimulus levels x number of chirps).
%
%     'noplot'  Do not plot main figure (fig 6 or 7).
% 
%     'plot'    Plot main figure (fig 6 or 7). See PLOT_ROENNE2012_CHIRP.
%
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%   See also: roenne2012, plot_roenne2012_chirp, exp_roenne2012
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
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/roenne2012_chirp.php

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
definput.flags.plot                 = {'noplot','plot','plot2'};
[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Init
fsmod               = 200000;       % AN model fs

% load Unitary Response and its samlping frequency
[ur,fs]=data_roenne2012;

% load Elberling et al (2010) click and chirp stimuli
[Elberlingstim,fsstim] = data_elberling2010('stim');

% length of modelling [ms]
modellength         = 40;           

% Output filter corresponding to recording settings (see Elberling et al
% 2010 section 3)
b=fir1(200,[100/(fs/2),3000/(fs/2)]);
a=1;

%% ABR model

% Loop over stimulus levels
for L = 1:length(levels)                                                    
  lvl  = levels(L);
  
  % Loop over chosen Elberling et al. (2010) stimuli
  for CE = chirps
    
    % Define length of stimulus, uses variable modellength
    stim  = zeros(modellength/1000*fsstim,1);
    
    % Use one of Elberling et al. (2010) stimuli
    stimCE = Elberlingstim.(['c' num2str(CE-1)]);
    
    % Create stimulus with chirp stimulus and concatenated zeros => combined 
    % length = "modellength"          
    stim(1:length(stimCE)) = stimCE;            
    
    % call AN model, note that lots of extra outputs are possible
    [ANout,vFreq] = zilany2007(lvl, stim', fsstim, fsmod);
    
    % Subtract 50 due to spontaneous rate
    ANout = ANout'-50;                                    
    
    % Sum in time across fibers = summed activity pattern        
    ANsum = sum(ANout,2);
    
    % Downsample ANsum to get fs = fs_UR = 30kHz
    ANsum = resample(ANsum ,fs,fsmod);                   
    
    % Simulated potential = UR * ANsum (* = convolved)
    simpot = filter(ur,1,ANsum);                           
    
    % apply output filter similar to the recording conditions in
    % Elberling et al. (2010)
    simpot = filtfilt(b,a,simpot);                         
    
    % Find max peak value (wave V)
    maxpeak(CE,L) = max(simpot);                             
    
    % Find corresponding time of max peak value (latency of wave V)
    waveVlat(CE,L) = find(simpot == maxpeak(CE,L));                      
    
    % find minimum in the interval from "max peak" to 6.7 ms later
    minpeak(CE,L) = min(simpot(find(simpot == maxpeak(CE,L)):find(simpot == maxpeak(CE,L))+200)); 
    
    %% PLOTS, extra plots created for all conditions used, i.e. three
    % plots for each stimulus level x each chirp sweeping rate. If this
    % is switched on and all other variables are set to default,
    % 3 (levels) x 6 (chirps / clicks) x 3 (different figures) = 48
    % figures will be created.
    if flags.do_plot2
      % Plot simulated ABR
      figure, t = 1e3.*[0:length(simpot)-1]./fs-15;
      plot(t,simpot,'k','linewidth',2),xlabel('time [ms]'), 
      title(['Simulated ABR (Elberling et al stimulus c' ...
             num2str(CE-1) ' at ' num2str(lvl) 'dB)']),
      set(gca, 'fontsize',12), axis([0 10 -.2 .3])
      
      % Plot "AN-gram" - spectrogram-like representation of the
      % discharge rate after the IHC-AN synapse 
      figure,  set(gca, 'fontsize',12),imagesc(ANout'), 
      title(['ANgram - (Elberling et al stimulus c' num2str(CE-1) ...
             ' at ' num2str(lvl) 'dB)'])
      set(gca,'YTick',[1 100 200 300 400 500]),
      set(gca,'YTicklabel',round(vFreq([1 100 200 300 400 500]))), 
      ylabel('model CF'), xlabel('time [ms]'),set(gca,'XTick',(0:500:8000)),
      set(gca,'XTicklabel',(0:500:8000)/fsmod*1000-15), xlabel('time [ms]'),
      xlim([12/1000*fsmod 26/1000*fsmod]),colorbar
      
      % Plot "AN-UR-gram" - spectrogram-like representation of the
      % discharge rate convolved line by line with the UR. 
      figure, ANUR = resample(ANout,fs,fsmod);
      ANUR = filter(ur,1,ANUR); imagesc(ANUR'),
      set(gca,'YTick',[1 100 200 300 400 500]),
      set(gca,'YTicklabel',round(vFreq([1 100 200 300 400 500]))), 
      ylabel('model CF'), xlabel('time [ms]'),
      set(gca,'XTick',(0:150:1500)), xlabel('time [ms]'),
      set(gca,'XTicklabel',(0:150:1500)/fs*1000-15)
      colorbar, xlim([450 1000]), 
      title(['AN-URgram (Elberling et al stimulus c' num2str(CE-1) ' stimulus at ' num2str(lvl) 'dB)'])
    end
  end
  
end

% Calculate wave V amplitude, as the difference between the peak and
% the dip.
waveVamp = (maxpeak-minpeak);

% Subtract 15 ms as click stimulus peaks 15 ms into the time series
waveVlat = waveVlat*1000/fs-15;


if flags.do_plot
  figure;
  plot_roenne2012_chirp(waveVamp, waveVlat);
end

