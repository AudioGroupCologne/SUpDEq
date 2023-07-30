function [outsig, fc] = dau1996(insig, fs, varargin)
%DAU1996   Linear filtering for monaural masking (basic)
%   Usage: [outsig, fc] = dau1996(insig,fs);
%          [outsig, fc] = dau1996(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%
%   Output parameters:
%     outsig : output signal.
%     fc     : center frequencies
%
%
%   DAU1996(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz as described in Dau,
%   Puschel and Kohlrausch (1996a).
%  
%   [outsig,fc]=DAU1996(...) additionally returns the center frequencies of
%   the filter bank.
%
%   Any of the optional parameters for AUDITORYFILTERBANK, IHCENVELOPE
%   and ADAPTLOOP may be specified for this function. They will be passed
%   to the corresponding functions.
%
%   The model implemented in this file is not identical to the model
%   published in Dau et. al. (1996a). An overshoot limit has been added to
%   the adaptation stage to fix a problem where abrupt changes in the
%   input signal would cause unnaturally big responses. This is described
%   in Dau et. al. (1997a).
%
%
%   The Dau 1996 model consists of the following stages:
%   
%   1) a gammatone filter bank with 1-erb spaced filtes.
%
%   2) an envelope extraction stage done by half-wave rectification followed by low-pass filtering to 1000 Hz.
%
%   3) an adaptation stage modelling nerve adaptation by a cascade of 5 loops.
%
%   4) a modulation low pass filter liming modulations to below 50 Hz.
%
%
%   *Warning:* This code is incorrect, the Dau 1996 models uses a
%   transmission-line model from Strube 1985, while this code erroneously
%   uses the gammatone filters. If/when the Strube model is included in
%   AMToolbox, this function will be fixed.
%
%   References:
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     J. Acoust. Soc. Am., 99(6):3615--3622, 1996a.
%     
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     "effective" signal processing in the auditory system. II. Simulations
%     and measurements. J. Acoust. Soc. Am., 99:3623--3631, 1996b.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%
%   See also: auditoryfilterbank, ihcenvelope, adaptloop, dau1997
%             plot_audspecgram demo_lopezpoveda2001
%             baumgartner2016_spectralanalysis
%             exp_osses2022 baumgartner2013
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/dau1996.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Torsten Dau 
%   #Author: Morten Løve Jepsen 
%   #Author: Peter L. Søndergaard (2011)

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

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop'};
definput.importdefaults={'ihc_dau1996','adt_dau1996'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

% Calculate filter coefficients for the 20 ms (approx.eq to 8 Hz) modulation
% lowpass filter.
% This filter places a pole /very/ close to the unit circle.
mlp_a = exp(-(1/0.02)/fs);
mlp_b = 1 - mlp_a;
mlp_a = [1, -mlp_a];

% Apply the low-pass modulation filter.
outsig = filter(mlp_b,mlp_a,outsig);

% Apply final resampling to avoid excessive data
if ~isempty(keyvals.subfs)
  outsig = fftresample(outsig,round(length(outsig)/fs*subfs));
end



