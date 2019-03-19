function [outsig, fc] = auditoryfilterbank(insig, fs, varargin);
%AUDITORYFILTERBANK   Linear auditory filterbank
%   Usage: [outsig, fc] = auditoryfilterbank(insig,fs);
%          [outsig, fc] = auditoryfilterbank(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
%   AUDITORYFILTERBANK(insig,fs) applies an auditory filterbank to the
%   imput signal insig sampled with a frequency of fs Hz. The filterbank
%   is composed of gammatone filters with 1 ERB wide filters.
%  
%   [outsig,fc]=AUDITORYFILTERBANK(...) additionally returns the center
%   frequencies of the filter bank.
%
%   The following parameters may be passed at the end of the line of
%   input arguments:
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 80 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 8000 Hz.
%
%     'basef',basef  Ensure that the frequency basef is a center frequency
%                    in the filterbank. The default value of [] means
%                    no default.
%
%   AUTHOR : Peter L. Søndergaard
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/auditoryfilterbank.php

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

definput.import={'auditoryfilterbank'};

[flags,keyvals,flow,fhigh]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(flow, fhigh, keyvals.bwmul, keyvals.basef);

% Calculate filter coefficients for the gammatone filter bank.
[gt_b, gt_a]=gammatone(fc, fs, 'complex');

% Apply the Gammatone filterbank
outsig = 2*real(ufilterbankz(gt_b,gt_a,insig));


