function [outsig, fc, mfc, params] = osses2021(insig, fs, varargin);
%OSSES2021   Monaural perceptual similarity
%
%   Usage: [outsig, fc] = osses2021(insig,fs);
%          [outsig, fc] = osses2021(insig,fs,...);
%          [outsig, fc, params] = osses2021(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%
%   Output parameters:
%     outsig  : output acoustic signal.
%     fc      : filter bank center frequencies.
%  
%   OSSES2021(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz.
%  
%   [outsig,fc,mfc]=OSSES2021(...) additionally returns the center
%   frequencies of the filter bank and the center frequencies of the
%   modulation filterbank.
%  
%   The model consists of the following stages:
%   
%   1) an outer- and middle-ear filtering as used by Jepsen et al. 2008
%
%   2) a gammatone filter bank with 1-erb spaced filters.
%
%   3) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 770 Hz as used by Breebaart et al. 2001
%
%   4) an adaptation stage modelling nerve adaptation by a cascade of 5
%      loops using a limiter factor of 5 (Osses and Kohlrausch, 2021).
%
%   5) a modulation filterbank
%
%   Any of the optional parameters for AUDITORYFILTERBANK,
%   IHCENVELOPE and ADAPTLOOP may be optionally specified for this
%   function. They will be passed to the corresponding functions.
%
%   See also: auditoryfilterbank, ihcenvelope, adaptloop, modfilterbank
%             exp_osses2021 exp_osses2022 breebaart2001 lopezpoveda2001 dau1997
%
%   References:
%     A. Osses Vecchi and A. Kohlrausch. Perceptual similarity between piano
%     notes: Simulations with a template-based perception model. J. Acoust.
%     Soc. Am., 141(4), 2020.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. II. Spectral and temporal integration. J.
%     Acoust. Soc. Am., 102:2906--2919, 1997b.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/osses2021.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author : Alejandro Osses (2022)

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
% load defaults from arg_auditoryfilterbank, arg_ihcenvelope, arg_adaptloop, arg_modfilterbank and arg_osses2021
definput.import={'auditoryfilterbank','ihcenvelope','adaptloop','modfilterbank','osses2021'}; 
definput.importdefaults={'afb_osses2021','ihc_breebaart2001', 'adt_osses2021','mfb_jepsen2008'}; 
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

fc  = [];
mfc = [];
params = [];
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

if flags.do_no_outerear
    if keyvals.silent_mode == 0
        %fprintf('\t%s: outer-ear filtering by-passed\n',upper(mfilename));
        amt_disp([upper(mfilename),': outer-ear filtering by-passed']);
        amt_disp();
    end
end 

if flags.do_middleear || flags.do_no_middleear
    filtertype = 'lopezpoveda2001';
elseif flags.do_jepsen2008
    filtertype = 'jepsen2008';
end
me_fir     = middleearfilter(fs,filtertype);
me_gain_TF = max( 20*log10(abs(freqz(me_fir,1,8192))) ); % max of the filter response

if flags.do_middleear || flags.do_jepsen2008
    N = ceil(length(me_fir)/2); % group delay for a FIR filter of order length(me_fir)
    M = 1; % assumes insig is monaural
    insig = [insig; zeros(N,M)]; % group delay compensation: step 1 of 2.
    insig = filter(me_fir,1,insig); % filtering
    insig = insig(N+1:end,1:M); % group delay compensation: step 2 of 2. 
    insig = gaindb(insig,-me_gain_TF); % if me_fir is a non-unit gain filter, 
                                       % the gain of the FIR filter is compensated.
    if keyvals.silent_mode == 0
        %fprintf('\t%s: middle-ear filter was adjusted to have a 0 dB bandpass gain (gain applied=%.1f dB)\n',upper(mfilename),-me_gain_TF);
        amt_disp();
        amt_disp([upper(mfilename),': middle-ear filter was adjusted to have a 0 dB bandpass gain (gain applied=',-me_gain_TF,' dB']);
        amt_disp();
    end
end                                    

if flags.do_no_middleear
    if keyvals.silent_mode == 0
        %fprintf('\t%s: middle-ear filtering by-passed\n',upper(mfilename));
        amt_disp();
        amt_disp([upper(mfilename),': middle-ear filtering by-passed']);
        amt_disp();
    end
end

if flags.do_afb
    % Apply the auditory filterbank
    [outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);
end

if flags.do_no_afb
    outsig = insig;
    if keyvals.silent_mode == 0
        %fprintf('\t%s: Gammatone filter bank by-passed\n',upper(mfilename));
        amt_disp();
        amt_disp([upper(mfilename),': Gammatone filter bank by-passed']);
        amt_disp();
    end
end

if flags.do_ihc || flags.do_adt || flags.do_mfb && ~flags.do_no_ihc
    % 'haircell' envelope extraction
    outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
else
    if keyvals.silent_mode == 0
        %fprintf('\t%s: ihcenvelope processing by-passed\n',upper(mfilename));
        amt_disp();
        amt_disp([upper(mfilename),': ihcenvelope processing by-passed']);
        amt_disp();
    end
end

if flags.do_adt || flags.do_mfb && ~flags.do_no_adt
    % non-linear adaptation loops
    outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);
else
    if keyvals.silent_mode == 0
        %fprintf('\t%s: adaptation loops by-passed\n',upper(mfilename));
        amt_disp();
        amt_disp([upper(mfilename),': adaptation loops by-passed']);
        amt_disp();
    end
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
    if nargout >= 4
        [outsig,mfc,params] = modfilterbank(outsig,subfs,fc,'argimport',flags,keyvals);
    else
        [outsig,mfc] = modfilterbank(outsig,subfs,fc,'argimport',flags,keyvals);
    end
end

if flags.do_no_mfb
    if keyvals.silent_mode == 0
        %fprintf('\t%s: modulation filter bank processing by-passed\n',upper(mfilename));
        amt_disp();
        amt_disp([upper(mfilename),': modulation filter bank processing by-passed']);
        amt_disp();
    end
end


