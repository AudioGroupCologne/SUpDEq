function [outsig, fc, mfc, step] = king2019(insig,fs,varargin)
%KING2019 Modulation filterbank (based on nonlinear processing)
%   
%   Usage:    outsig = king2019(insig,fs,basef)
%             [outsig, fc, mfc, step] = king2019(insig,fs,varargin)
%              [outsig, fc, mfc, step]= king2019(insig,fs,flow,fhigh,parameters)
% 
%   Input parameter:
%     insig  : input acoustic signal.
%     fs     : sampling rate (Hz).
%     basef  : base frequency of the analysis (Hz).
%
%   Output parameters:
%     outsig : output signal
%     fc     : center frequencies filterbank
%     mfc    : center frequencies modulation filterbank
%     step   : struct containing intermediate model outputs
%
%   KING2019(insig,fs,'basef',basef) computes the internal representation 
%   around the frequency basef of the signal insig sampled with 
%   a frequency of fs Hz. outsig is a matrix of size [length fc mfc]. 
%  
%   [outsig,fc,mfc,step]=KING2019(...) additionally returns the 
%   centre frequencies of the filter bank and the center frequencies of the
%   modulation filterbank, and 'steps' is a structure containing the 
%   intermediate model outputs.
%  
%   The model consists of the following stages:
%
%   1) a gammatone filter bank with 1-erb spaced filtes.
%
%   2) a compression (exponential or broken-stick) stage applied to each 
%      individual gammatone filter. In these stages is extremely relevant
%      to use the correct calibration level (default is dboffset = 100), given
%      that the keyval 'compression_knee_dB' is referenced to the dboffset.
%      Note that if other dboffset values are used (e.g., dboffset=94 dB), the
%      knee point will in fact start compressing at higher levels.
%
%   3) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 1500 Hz.
%
%   4) a simplified adaptation stage simulated as a first-order high-pass 
%      filter
%
%   5) a modulation filterbank
%
%   Many parameters (keyvals) and flags can be optionally changed and/or 
%   switched 
%
%
%   References:
%     A. King, L. Varnet, and C. Lorenzi. Accounting for masking of frequency
%     modulation by amplitude modulation with the modulation filter-bank
%     concept. J. Acoust. Soc. Am., 145(2277), 2019.
%     
%
%   See also: auditoryfilterbank, ihcenvelope, king2019_modfilterbank
%             demo_king2019 exp_osses2022 dau1996
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/king2019.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Leo Varnet and Andrew King (2020) as king2019_preproc
%   #Author: Alejandro Osses (2020) Further adaptations
%   #Author: Clara Hollomey (2020) Adapted for AMT as king2019
%   #Author: Piotr Majdak (2021) Further adaptations to AMT 1.0
%   #Author: Alejandro Osses (2023) Updated documentation

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% ------ Checking of input parameters ------------
if nargin<3
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import={'auditoryfilterbank','ihcenvelope','king2019'}; % load from arg_king2019

[flags,kv]  = ltfatarghelper({'basef'},definput,varargin);

fc = [];
mfc = [];
if isempty(kv.flow) && isempty(kv.fhigh)
    if isempty(kv.basef)
        error('%s: Please specify the centre frequency of your analysis, type ''help %s''.',upper(mfilename),mfilename);
    end
end

step_erb = 2;
if isempty(kv.flow)
    kv.flow = floor(audtofreq(freqtoaud(kv.basef)-step_erb));
end
if isempty(kv.fhigh)
    kv.fhigh = ceil(audtofreq(freqtoaud(kv.basef)+step_erb));
end

% ------ do the computation -------------------------
if nargout >= 4
    step = [];
end

amt_disp('KING2019:',flags.disp);

Nsamples = length(insig);
% -------------------------------------------------------------------------
% --- Filter bank stages:
if flags.do_afb
    amt_disp('  Calculating the auditory filterbank...',flags.disp);
    % 'Oldenburg' gammatone filterbank
    [outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,kv);
    outsig    = squeeze(outsig);
    Nchannels = length(fc); % Number of auditory bands
    Nsamples  = size(outsig,1); % Number of samples per band
end
if flags.do_no_afb
    amt_disp('  Auditory filterbank skipped',flags.disp);
    outsig = insig;
end

if nargout >= 4
    step.gtone_response = outsig;
end

% -------------------------------------------------------------------------
% --- "Power" compression 
% WARNING: below the kneepoint this option actually performs an *expansion*
% of the signal

if flags.do_compression_power
    amt_disp('  Power compression...',flags.disp);
    
    comp_n = kv.compression_n;
    if length(comp_n) == 1
        comp_n = comp_n*ones(1,Nchannels);
    elseif length(comp_n) ~= Nchannels
        error('%s: compression_n does not have the same number of channels as the output of the Gammatone Filterbank',upper(mfilename));
    end
    dBFS = kv.dboffset; % dB full scale
    comp_knee = gaindb(1,kv.compression_knee_dB-dBFS);
        
    outsig = sign(outsig).*abs(outsig/comp_knee).^(ones(Nsamples,1)*comp_n)*comp_knee;
end

% -------------------------------------------------------------------------
% --- "Brockenstick" compression
if flags.do_compression_brokenstick
    amt_disp('  Broken-stick compression...',flags.disp); 
        
    comp_n = kv.compression_n;
    if length(comp_n)==1
    	comp_n = comp_n*ones(1,Nchannels);
    elseif length(comp_n) ~= Nchannels
        error('%s: compression_n does not have the same number of channels as the output of the Gammatone Filterbank',upper(mfilename));
    end
    dBFS = kv.dboffset; % dB full scale
    comp_knee = gaindb(1,kv.compression_knee_dB-dBFS);
    
    outsig = local_broken_stick(outsig,comp_knee,comp_n);
end

if nargout >= 4
    step.compressed_response = outsig;
end

% -------------------------------------------------------------------------
% --- 'haircell' envelope extraction
if (flags.do_ihc || flags.do_adt || flags.do_mfb) && ~flags.do_no_ihc
    % do_ihc is forced to be done if do_ihc, do_adaptation, or do_mfb are
    %   active modules, EXCEPT that the user explicitly deactivated the module
    %   (do_no_ihc or do_noihc).
    amt_disp('  Hair-cell envelope extraction',flags.disp); 
    
    outsig = ihcenvelope(outsig,fs,'argimport',flags,kv);
    
    if nargout >= 4
        step.ihc = outsig;
    end
end

% -------------------------------------------------------------------------
% --- Adaptation by high-pass filtering
if flags.do_adt || flags.do_mfb && ~flags.do_no_adt
    % do_adaptation is forced to be done if do_adaptation or do_mfb are active
    %   modules, EXCEPT that the user explicitly deactivated the module (do_no_adt).
    amt_disp('  Adaptation by high-pass filtering',flags.disp);
    adt_HP_fc = kv.adt_HP_fc; 
    adt_HP_order = kv.adt_HP_order;
    [b,a] = butter(adt_HP_order,2*(adt_HP_fc/fs),'high');

    outsig=filter(b,a,outsig);
    
    if nargout >= 4
        step.a_adapt_HP = a;
        step.b_adapt_HP = b;
        step.adapted_response = outsig;
    end
end

% -------------------------------------------------------------------------
% --- Modulation filterbank
%%% Modulation filterbank
if flags.do_mfb
    amt_disp('  Modulation filter bank',flags.disp);
    
    % % Parameters modulation filter:
    [outsig, mfc, step_mod] = king2019_modfilterbank(outsig, fs,'argimport',flags,kv);
        
    if nargout >= 4
        step.a_mfb = step_mod.a;
        step.b_mfb = step_mod.b;
        step.E_mod = step_mod.E_mod;
        step.fmc = mfc;
    end
end

% -------------------------------------------------------------------------
% --- Downsampling (of the internal representations)
%  Apply final resampling to avoid excessive data
if ~isempty(kv.subfs) && flags.do_mfb
    
    subfs = kv.subfs;
    if subfs ~= fs
        amt_disp(['  Downsampling to ' num2str(subfs) ' Hz'],flags.disp); 
        % In case of downsampling:
        outsig = fftresample(outsig,round(length(outsig)/fs*subfs));
    end
    
else
    % In case of no-resampling:
    subfs = fs;
end
step.subfs = subfs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = local_broken_stick(insig,comp_knee,comp_n)
% This function does the same as nltrans3.m from Stephan Ewert with no smoothing

Nchannels = size(insig,2);

outsig = insig;
for j = 1:Nchannels 
    if comp_n(j) ~= 1
        idxs = find(abs(insig(:,j))>comp_knee);
        outsig(idxs,j) = sign(insig(idxs,j)).*(abs(insig(idxs,j)).^comp_n(j) / comp_knee^comp_n(j) * comp_knee);
    end
end


