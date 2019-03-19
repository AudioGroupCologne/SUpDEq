function [waveVamp, waveVlat] = exp_roenne2012(varargin)
%EXP_ROENNE2012 Figures from Rønne et al. (2012)
%   Usage: output = exp_roenne2012(flag)
%
%   EXP_ROENNE2012(flag) reproduces the results for the figure given
%   by flag from the Rønne et al. (2012) paper. Outputs are the ABR wave V
%   amplitude and latency of all datapoints in that given figure.
%   
%   The following flags can be specified;
%
%     'plot'     Plot the specified figure from Rønne et al. (2012). This is
%                the default. 
%
%     'noplot'   Don't plot, only return data.
%
%     'plot2'    Plot extra figures for all individual simulated points.
%                Note that this creates lots of extra figures (3 for each
%                simulated data point)
%
%     'redo'     Recalculate the experiment results.
%
%     'cached'   Use cached results. Default.
%
%     'fig5'     Plot Fig. 5 (Rønne et al., 2012). Latency of simulated ABR
%                wave V's compared to Neely et al. (1988) and Harte et al.
%                (2009) reference data.
%
%     'fig6'     Plot Fig. 6 (Rønne et al., 2012). Amplitude of simulated
%                wave V compared to Elberling et al. (2010) reference data.
%
%     'fig7'     Plot Fig. 7 (Rønne et al., 2012). Latency of simulated wave
%                V compared to Elberling et al. (2010) reference data.
%
%   Examples:
%   ---------
%
%   To display Figure 5 use :
%
%     exp_roenne2012('fig5');
%
%   To display Figure 6 use :
%
%     exp_roenne2012('fig6');
%
%   To display Figure 7 use :
%
%     exp_roenne2012('fig7');
%
%   References:
%     C. Elberling, J. Calloe, and M. Don. Evaluating auditory brainstem
%     responses to different chirp stimuli at three levels of stimulation. J.
%     Acoust. Soc. Am., 128(1):215-223, 2010.
%     
%     J. Harte, G. Pigasse, and T. Dau. Comparison of cochlear delay
%     estimates using otoacoustic emissions and auditory brainstem responses.
%     J. Acoust. Soc. Am., 126(3):1291-1301, 2009.
%     
%     S. Neely, S. Norton, M. Gorga, and J. W. Latency of auditory brain-stem
%     responses and otoacoustic emissions using tone-burst stimuli. J.
%     Acoust. Soc. Am., 83(2):652-656, feb 1988.
%     
%     F. M. RÃ¸nne, T. Dau, J. Harte, and C. Elberling. Modeling auditory
%     evoked brainstem responses to transient stimuli. The Journal of the
%     Acoustical Society of America, 131(5):3903-3913, 2012. [1]http ]
%     
%     References
%     
%     1. http://scitation.aip.org/content/asa/journal/jasa/131/5/10.1121/1.3699171
%     
%
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_roenne2012.php

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
definput.import={'amt_cache'};
definput.flags.type={'missingflag','fig5','fig6','fig7'};
definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
		flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}), ...
				sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
		error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5
  
  waveVamp = 0;

  stim_level = 40:10:100; % Default stimulus levels

  [waveVlat, click_latency] = amt_cache('get','fig5',flags.cachemode);  
  if isempty(waveVlat);
    [click_amplitude, click_latency]    = roenne2012_click(stim_level);     
    waveVlat = roenne2012_tonebursts(stim_level);    
    amt_cache('set','fig5',waveVlat,click_latency);
  end;
  
  if flags.do_plot;
    plot_roenne2012_tonebursts(waveVlat,click_latency);
  end  ;  

end;

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6;

  stim_level    = (20:20:60)+35.2;
  
  % Default chirp numbers. 1 = click, 2 to 6 = chirp 1 to 5.
  chirp_number  = 1:6;              

  [waveVamp, waveVlat] = amt_cache('get','fig6',flags.cachemode);
  if isempty(waveVamp)
    [waveVamp, waveVlat] = roenne2012_chirp(stim_level, chirp_number);
    amt_cache('set','fig6',waveVamp,waveVlat);
  end;
        
  if flags.do_plot
    plot_roenne2012_chirp(waveVamp, waveVlat,'amponly');
  end
end;

%% ------ FIG 7 -----------------------------------------------------------
if flags.do_fig7;

  stim_level    = (20:20:60)+35.2;
  
  % Default chirp numbers. 1 = click, 2 to 6 = chirp 1 to 5.
  chirp_number  = 1:6;              

  [waveVamp, waveVlat] = amt_cache('get','fig7',flags.cachemode);
  if isempty(waveVamp)
    [waveVamp, waveVlat] = roenne2012_chirp(stim_level, chirp_number);
    amt_cache('set','fig7',waveVamp,waveVlat);
  end;
        
  if flags.do_plot
    plot_roenne2012_chirp(waveVamp, waveVlat,'latonly');
  end
end;


