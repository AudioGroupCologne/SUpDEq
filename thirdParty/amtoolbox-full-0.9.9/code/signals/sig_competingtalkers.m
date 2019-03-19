function [s,fs]=sig_competingtalkers(varargin)
%sig_competingtalkers  Load one of several test signals
%   Usage:  s=sig_competingtalkers(signame);
%           [s,fs]=sig_competingtalkers(signame);
%
%   SIG_COMPETINGTALKERS(signame) loads one of several test signals consisting
%   of competing talkers. All the talkers are taken from the TIMIT speech
%   corpus:
%   http://www.ldc.upenn.edu/Catalog/CatalogEntry.jsp?catalogId=LDC93S1.
%
%   The signals have 2 channels and are all recorded with a sampling rate of
%   16 kHz.
%
%   [sig,fs]=SIG_COMPETINGTALKERS(signame) additionally returns the sampling
%   frequency fs.
%
%   The value of signame may be one of:
%
%     'one_of_three'    XXX Description missing
%
%     'two_of_three'    XXX Description missing
%
%     'three_of_three'  XXX Description missing
%
%     'one_speaker_reverb'
%                       XXX Description missing
%
%     'two_speakers'    XXX Description missing
%
%     'five_speakers'   XXX Description missing
%
%     'bnoise'          Speech shaped noise
%
%   Examples:
%   ---------
%
%   The following plot shows an estimate of the power spectral density of
%   the first channels of the speech shaped noise:
%
%      s=sig_competingtalkers('bnoise');
%      pwelch(s(:,1),hamming(150));
%
%   See also: exp_dietz2011
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_competingtalkers.php

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

%   AUTHOR : Peter L. Søndergaard


definput.flags.sigtype={'missingflag','one_of_three','two_of_three',...
                    'three_of_three','one_speaker_reverb',...
                    'two_speakers','five_speakers','bnoise'};

[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.sigtype{2:end-2}),...
             sprintf('%s or %s',definput.flags.sigtype{end-1},definput.flags.sigtype{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


% fs = 16000;
[s,fs]=amt_load('sig_competingtalkers',[flags.sigtype '.wav']);
