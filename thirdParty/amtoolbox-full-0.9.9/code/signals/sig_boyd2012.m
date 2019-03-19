function out = sig_boyd2012(in,source,varargin)
%sig_boyd2012 - Stimulus from Boyd et al. (2012) 
%   Usage: out = sig_boyd2012(in,source,mix,azi)
%
%   Time-aligned mixture between individualized binaural stimulus and 
%   head-absent simulation (ITD only).
%
%   Input parameters:
%     in      : binaural input signal.
%     source	: monaural source signal.
%     mix     : mixing ratio ranging from 0 to 1. Default is 1 and means 
%               out is same as in.
%     lp      : low-pass cut off frequency. Default is NAN and means
%               broadband.
%     fs      : sampling rate in Hz (only required for low-pass filtering).
%               Default is 48 kHz.
%     azi     : azimuth (positive to the left).
%
%   Output parameters:
%     out     : binaural output signal mixture.
%
%   References:
%     A. W. Boyd, W. M. Whitmer, J. J. Soraghan, and M. A. Akeroyd. Auditory
%     externalization in hearing-impaired listeners: The effect of pinna cues
%     and number of talkers. J. Acoust. Soc. Am., 131(3):EL268-EL274, 2012.
%     [1]arXiv | [2]www: ]
%     
%     References
%     
%     1. http://arxiv.org/abs/%20http://dx.doi.org/10.1121/1.3687015
%     2. http://dx.doi.org/10.1121/1.3687015
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_boyd2012.php

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

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

definput.keyvals.mix = 2;
definput.keyvals.lp = nan;
definput.keyvals.fs = 48e3;
definput.keyvals.azi = -30;

[flags,kv]  = ltfatarghelper({'mix','lp','fs','azi'},definput,varargin);

% create time-aligned head-absent simulation
[scor,lag] = xcorr(in(:,1),source(:));
[~,iL] = max(abs(scor));
iL = lag(iL);
[~,iR] = max(abs(xcorr(in(:,2),source(:))));
iR = lag(iR);
ha = zeros(length(in),2);
ha(iL+(1:length(source)),1) = source;
ha(iR+(1:length(source)),2) = source;
if kv.azi > 0
  ha = fliplr(ha);
end
SPL = mean(dbspl(in));
ha = setdbspl(ha,SPL);
   
% mixing
out = kv.mix*in + (1-kv.mix)*ha;

% low-pass filtering
if not(isnan(kv.lp) || isempty(kv.lp))
  [b,a]=butter(10,2*kv.lp/kv.fs,'low');
  out = filter(b,a,out);
end

end
