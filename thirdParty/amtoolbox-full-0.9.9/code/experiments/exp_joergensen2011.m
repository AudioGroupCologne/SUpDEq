function exp_joergensen2011(varargin)
%EXP_JOERGENSEN2011 Figures from Joergensen and Dau (2011)
%   Usage: output = exp_joergensen2011(flag)
%
%   EXP_JOERGENSEN2011(flag) reproduces the results for the figure given
%   by flag from the Joergensen and Dau (2011) paper. 
%    
%   The following flags can be specified;
%
%     'plot'     Plot the specified figure from Jørgensen and Dau (2011). This is
%                the default. 
%
%     'noplot'   Don't plot, only return data.
%
%
%     'redo'     Always recalculate the experiment results.
%
%     'cached'   Always use the cached version. Default.
%
%     'fig5'     Plot Fig. 5 (Joergensen and Dau, 2011).
%
%     'fig6'     Plot Fig. 6 (Joergensen and Dau, 2011). 
% 
%
%   Examples:
%   ---------
%
%   To display Figure 5 use :
%
%     exp_joergensen2011('fig5');
%
%   To display Figure 6 use :
%
%     exp_joergensen2011('fig6');
%
%   ---------
%
%   Please cite Joergensen and Dau (2011) if you use this model.
%
%   See also: joergensen2011, plot_joergensen2011, exp_joergensen2011
%
%   References:
%     S. Joergensen and T. Dau. Predicting speech intelligibility based on
%     the signal-to-noise envelope power ratio after modulation-frequency
%     selective processing. J. Acoust. Soc. Am., 130(3):1475-1487, 2011.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_joergensen2011.php

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
definput.flags.type={'missingflag','fig5','fig6'};
definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

NSpeechsamples = 10; % specify the number of speech samples to be used fro the simulations. The simulation takes longer the more samples are used. A minimum of 50 should be used for final validation.  

%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5;

  dSRT = amt_cache('get', ['fig5_' num2str(NSpeechsamples) 'sntcs'], flags.cachemode);

  if isempty(dSRT)
    dSRT = joergensen2011_sim(NSpeechsamples,'fig5');
    amt_cache('set',['fig5_' num2str(NSpeechsamples) 'sntcs'],dSRT);
  end;
        
  if flags.do_plot
    plot_joergensen2011(dSRT,'fig5');
  end
end;

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6;

  dSRT = amt_cache('get', ['fig6_' num2str(NSpeechsamples) 'sntcs'], flags.cachemode);
  if isempty(dSRT)

    dSRT = joergensen2011_sim(NSpeechsamples,'fig6');
		amt_cache('set',['fig6_' num2str(NSpeechsamples) 'sntcs'],dSRT);
  end;
        
  if flags.do_plot
    plot_joergensen2011(dSRT,'fig6');
  end
end;

