function [outsig, fc, mfc] = dau1997(insig, fs, varargin);
%DAU1997   Linear filtering for monaural masking (improved)
%   Usage: [outsig, fc] = dau1997(insig,fs);
%          [outsig, fc] = dau1997(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%
%   Output parameters:
%     outsig : output signal
%     fc     : center frequencies [Hz]
%  
%   DAU1997(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz.
%  
%   [outsig,fc,mfc]=DAU1997(...) additionally returns the center
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
%             plot_audspecgram relanoiborra2019 exp_osses2021
%             exp_osses2022 lopezpoveda2001 dau1996
%
%   References:
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. II. Spectral and temporal integration. J.
%     Acoust. Soc. Am., 102:2906--2919, 1997b.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/dau1997.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author : Torsten Dau, 
%   #Author: Morten Løve Jepsen 
%   #Author: Peter L. Søndergaard (2011)
%   #Author: Alejandro Osses (2021)

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

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop','modfilterbank'};
definput.importdefaults={'afb_dau1997', 'ihc_dau1996', 'adt_dau1997','mfb_jepsen2008'}; 
                    % The flag 'mfb_dau1997' would be exactly as in day1997mapI
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------
insig = gaindb(insig,keyvals.dboffset-100); % from here on, the input signal is
                               % assumed to be at a dboffset of 100 dB (default AMT)
if flags.do_outerear
  hp_fir = headphonefilter(fs);% Getting the filter coefficients at fs
  N = ceil(length(hp_fir)/2);  % group delay for a FIR filter of order length(hp_fir)
  M = 1; % assumes insig is monaural
  insig = [insig; zeros(N,M)]; % group delay compensation: step 1 of 2. 
  insig = filter(hp_fir,1,insig); % filtering
  insig = insig(N+1:end,1:M); % group delay compensation: step 2 of 2
end

if flags.do_middleear || flags.do_jepsen2008
  if flags.do_middleear
      filtertype = 'lopezpoveda2001';
  elseif flags.do_jepsen2008
      filtertype = 'jepsen2008';
  end
  me_fir = middleearfilter(fs,filtertype);
  N = ceil(length(me_fir)/2); % group delay for a FIR filter of order length(me_fir)
  M = 1; % assumes insig is monaural
  insig = [insig; zeros(N,M)]; % group delay compensation: step 1 of 2.
  insig = filter(me_fir,1,insig); % filtering
  insig = insig(N+1:end,1:M); % group delay compensation: step 2 of 2. 
  me_gain_TF = max( 20*log10(abs(freqz(me_fir,1,8192))) ); % max of the filter response
  insig = gaindb(insig,-me_gain_TF); % if me_fir is a non-unit gain filter, 
                                     % the gain of the FIR filter is compensated.
end

% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

if flags.do_ihc
    % 'haircell' envelope extraction
    outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
end

if flags.do_adt || flags.do_mfb
    % non-linear adaptation loops
    outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);
end

if flags.do_mfb
    %% Downsampling (of the internal representations)
    % Apply final resampling to avoid excessive data
    if ~isempty(keyvals.subfs)
        % In case of downsampling:
        outsig = fftresample(outsig,round(length(outsig)/fs*keyvals.subfs));
        subfs = keyvals.subfs;
    else
        % In case of no-resampling:
        subfs = fs;
    end

    % Modulation filterbank
    [outsig,mfc] = modfilterbank(outsig,subfs,fc,'argimport',flags,keyvals);
end


