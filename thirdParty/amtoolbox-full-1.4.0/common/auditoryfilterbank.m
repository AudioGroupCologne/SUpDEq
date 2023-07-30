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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/auditoryfilterbank.php


%   #Author : Peter L. Søndergaard

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
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
if flags.do_default
    % find the center frequencies used in the filterbank, 1 ERB spacing
    fc = erbspacebw(flow, fhigh, keyvals.bwmul, keyvals.basef);

    % Calculate filter coefficients for the gammatone filter bank.
    [gt_b, gt_a]=gammatone(fc, fs, 'complex');

    % Apply the Gammatone filterbank
    outsig = 2*real(ufilterbankz(gt_b,gt_a,insig));
else    
    fc = varargin{1};
    outsig = comp_auditoryfilterbank_singlefc(insig, fs, fc);
end


