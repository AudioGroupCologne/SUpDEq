function [outsig, fc, mfc] = dau1997_preproc(insig, fs, varargin);
%dau1997_preproc   Auditory model from Dau et. al. 1997
%   Usage: [outsig, fc] = dau1997_preproc(insig,fs);
%          [outsig, fc] = dau1997_preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
%   DAU1997_PREPROC(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz.
%  
%   [outsig,fc,mfc]=DAU1997_PREPROC(...) additionally returns the center
%   frequencies of the filter bank and the center frequencies of the
%   modulation filterbank.
%  
%   The model consists of the following stages:
%   
%   1) a gammatone filter bank with 1-erb spaced filtes.
%
%   2) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 1000 Hz.
%
%   3) an adaptation stage modelling nerve adaptation by a cascade of 5
%      loops.
%
%   4) a modulation filterbank
%
%   Any of the optinal parameters for AUDITORYFILTERBANK,
%   IHCENVELOPE and ADAPTLOOP may be optionally specified for this
%   function. They will be passed to the corresponding functions.
%
%   See also: auditoryfilterbank, ihcenvelope, adaptloop, modfilterbank
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/dau1997_preproc.php

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

%   References: dau1997mapI dau1997mapII

%   AUTHOR : Torsten Dau, Morten Løve Jepsen, Peter L. Søndergaard
  
% ------ Checking of input parameters ------------

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop'};
definput.importdefaults={'gtf_dau', 'ihc_dau', 'adt_dau'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);


