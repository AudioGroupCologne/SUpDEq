function out = sig_boyd2012(in,source,varargin)
%SIG_BOYD2012 - Stimulus from Boyd et al. (2012) 
%   Usage: out = sig_boyd2012(in,source,mix,azi)
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
%   Time-aligned mixture between individualized binaural stimulus and 
%   head-absent simulation (ITD only).
%
%
%   References:
%     A. W. Boyd, W. M. Whitmer, J. J. Soraghan, and M. A. Akeroyd. Auditory
%     externalization in hearing-impaired listeners: The effect of pinna cues
%     and number of talkers. J. Acoust. Soc. Am., 131(3):EL268--EL274, 2012.
%     [1]www: ]
%     
%     References
%     
%     1. http://dx.doi.org/10.1121/1.3687015
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_boyd2012.php


%   #Author: Robert Baumgartner (2021), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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
ha = scaletodbspl(ha,SPL, 100);
   
% mixing
out = kv.mix*in + (1-kv.mix)*ha;

% low-pass filtering
if not(isnan(kv.lp) || isempty(kv.lp))
  [b,a]=butter(10,2*kv.lp/kv.fs,'low');
  out = filter(b,a,out);
end

end


