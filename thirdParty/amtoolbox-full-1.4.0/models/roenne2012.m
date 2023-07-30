function [waveVamp, waveVlat, varargout]  = roenne2012(stim,fsstim,stim_level,varargin)
%ROENNE2012 Simulate auditory brainstem responses (ABRs)
%   Usage: [waveVamp, waveVlat]  = roenne2012(stim,fsstim,stim_level)
%
%   Input parameters:
%     stim         : stimuli
%     fsstim       : sampling frequency
%     stim_level   : level
%
%   Output parameters:
%     waveVamp   : Amplitude of simulated ABR wave V.
%     waveVlat   : Latency of simulated ABR wave V peak.
%
%
%   ROENNE2012(stim,fsstim,stim_level) returns simulated ABR wave V
%   latency and amplitude. The stimulus stim must be defined in pascals
%   and calibrated so a pure tone stimulus has an RMS value of 1. Transient
%   stimuli (which this model is designed to simulate) has to be calibrated
%   in peSPL acoustically. This is *not* the same as "just" having a
%   numerical peak to peak value of the same level as the pure tone. For
%   calibrated click, chirps and tone bursts, see ROENNE2012_CLICK,
%   ROENNE2012_TONEBURSTS and ROENNE2012_CHIRP.
%
%   The parameter fsstim gives the sampling frequency of the input
%   stimulus, and stim_level the level. As input is calibrated to an
%   RMS-value of 1, a stimulus level in (pe)SPL has to be set.
%
%   Flags:
%
%     'plot'                : Plot the output. See PLOT_ROENNE2012.
%  
%     'no_plot'             : Do not plot. This is the default.
%
%     'fsmod',fsmod         : Auditory nerve model sampling frequency.
%                             Default value is 200000.
%      
%     'flow',flow           : Auditory nerve model lowest center frequency.
%                             Default value is 100 Hz.
%
%     'fhigh',fhigh         : Auditory nerve model highest center frequency.
%                             Default value is 16000 Hz.
%
%     'min_modellength',mn  : Minimum length of modelling measured in ms.
%                             Default value is 40.
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%   See also: data_roenne2012 plot_roenne2012 plot_roenne2012_chirp
%             plot_roenne2012_tonebursts demo_roenne2012 roenne2012_click
%             roenne2012_chirp roenne2012_tonebursts exp_roenne2012 zilany2014
%
%   References:
%     C. Elberling, J. Calloe, and M. Don. Evaluating auditory brainstem
%     responses to different chirp stimuli at three levels of stimulation. J.
%     Acoust. Soc. Am., 128(1):215--223, 2010.
%     
%     F. M. Rønne, T. Dau, J. Harte, and C. Elberling. Modeling auditory
%     evoked brainstem responses to transient stimuli. The Journal of the
%     Acoustical Society of America, 131(5):3903--3913, 2012. [1]http ]
%     
%     M. S. A. Zilany and I. C. Bruce. Representation of the vowel (epsilon)
%     in normal and impaired auditory nerve fibers: Model predictions of
%     responses in cats. J. Acoust. Soc. Am., 122(1):402--417, jul 2007.
%     
%     References
%     
%     1. http://scitation.aip.org/content/asa/journal/jasa/131/5/10.1121/1.3699171
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/roenne2012.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Peter L. Sondergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% Define input flags
definput.flags.plot     = {'no_plot', 'plot'};
definput.keyvals.fsmod=200000;
definput.keyvals.flow = 100;
definput.keyvals.fhigh = 16000;
definput.keyvals.min_modellength=40;
[flags,kv]      = ltfatarghelper({},definput,varargin);

%% Init
[ur,fs] = data_roenne2012;

% Assure minimum model length of 40ms
if length(stim)/fsstim < kv.min_modellength/1000                               
    stim_temp = zeros(1, fsstim*kv.min_modellength/1000);
    stim_temp(1:length(stim)) = stim;
    stim = stim_temp;
end

%% ABR model
% call AN model, note that lots of extra outputs are possible
[ANout,vFreq] = zilany2007(stim_level, stim, fsstim, kv.fsmod, 'flow',kv.flow, 'fhigh',kv.fhigh);   

% subtract 50 due to spontaneous rate
ANout = ANout-50;                                            

% Sum in time across fibers, summed activity pattern
ANsum1 = sum(ANout,2);                 

% Downsample ANsum to get fs = fs_UR = 32kHz
ANsum = resample(ANsum1,fs,kv.fsmod); 

% Simulated potential = UR * ANsum (* = convolution)
simpot = filter(ur,1,ANsum);        

% Find max peak value (wave V)
maxpeak = max(simpot);                                          

% Find corresponding time of max peak value (latency of wave V). The unit
% is [ms]. 
waveVlat = find(simpot == maxpeak)/fs*1000; 

% find minimum in the interval from "max peak" to 6.7 ms later
minpeak = min(simpot(find(simpot == max(simpot)):...
                     find(simpot == max(simpot))+200)); 

% Calculate wave V amplitude, as the difference between the peak and the
% dip, in [\mu p] (micro pascals).
waveVamp = (maxpeak-minpeak);                               

if nargout >= 3
    varargout{1} = simpot;
    if nargout >= 4
        varargout{2} = ANout;
    end
end
end


