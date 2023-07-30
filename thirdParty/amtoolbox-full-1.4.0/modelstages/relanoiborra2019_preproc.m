function [outsig,fc,mfc] = relanoiborra2019_preproc(insig, fs, varargin)
%RELANOIBORRA2019_PREPROC the internal representations of the template and target signals
%   according to the CASP model (see references).
%
%   The code is based on previous versions of authors: Torsten Dau, Morten
%   Loeve Jepsen, Boris Kowalesky and Peter L. Soendergaard
%
%   Notes:
%
%   The model has been optimized to work with speech signals, and the
%   preprocesing and variable names follow this principle. The model is
%   also designed to work with broadband signals. In order to avoid undesired
%   onset enhancements in the adaptation loops, the model expects to recive a
%   prepaned signal to initialize them.
%
%   Inputs:
%
%      insig      :  input signal
%      fs         :  sampling frequency (Hz)
%      flow       :  lowest center frequency of auditory filterbank
%      fhigh      :  highest center frequency of auditory filterbank
%      sbj        :  subject profile for drnl definition
%
%   Outputs:
%      out           : correlation metric structure inlcuding:
%          .dint      : correlation values for each modulation band
%          .dsegments : correlation values from each time window and mod. band.
%          .dfinal    : final (averaged) correlation
%
%   References:
%     H. Relaño-Iborra, J. Zaar, and T. Dau. A speech-based computational
%     auditory signal processing and perception model. J. Acoust. Soc. Am.,
%     146(5), 2019.
%     
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am., 124(1),
%     2008.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/relanoiborra2019_preproc.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Stats M-Signal M-Control M-Communication Systems
%   #Author: Helia Relano Iborra (March 2019): v4.0 provided to the AMT team
%   #Author: Clara Hollomey (2021): adapted to the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   ------------ Change log: -------------------
%
%    v1.0 August 2017, Helia Relano Iborra:
%
%    - Added input N_org (number of time samples in the original sentence)
%    and added cutting of prepanned sentence after adaptation loops for
%    initialization.
%    - Changed the maximum modulation frequency for consistency with
%    Jepsen et al. 2008.
%
%   -------------------------------
%    v2.0 October 2018, Helia Relano Iborra:
%
%    - Added internal noise and removed minimum level check.
%    -------------------------------
%
%    v3.0 December 2018, Helia Relano Iborra:
%
%    - Changed drnl parameters for consistency with Jepsen (2008)
%    - Changed middle ear filter parameters from Lopez-Poveda & Meddis
%    (2001)
%    -------------------------------
%
%    v4.0 March 2019, Helia Relano Iborra:
%
%    - Adjust internal noise levels for NH free field thresholds
%    - 
%    -------------------------------

mfc = [];
fc = [];

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop','modfilterbank'};
definput.importdefaults={'drnl_relanoiborra2019','ihc_relanoiborra2019','adt_relanoiborra2019','mfb_jepsen2008'}; 
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

%% Outer- and middle-ear filtering
%C.H.==================================================================
if flags.do_outerear
    b_hp  = headphonefilter(fs); % calc headphone filtercoeffs
    insig = filter(b_hp,1,insig); % Outer-ear filterring
end
if flags.do_middleear
    b_me  = middleearfilter(fs);
    insig = filter(b_me,1,insig); % middle-ear-ear filterring
end

%%  -------  Case of DRNL

% % find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(keyvals.flow, keyvals.fhigh, keyvals.bwmul, keyvals.basef);

% -------------------------------------------------------------------------
nchannels = size(fc,2);

% Pre-allocate memory for auditory bands
outsig = zeros(length(insig),nchannels);

for n = 1:nchannels % Filter
    outsig(:,n) = relanoiborra2019_drnl(insig,fc(n),fs, keyvals.hearing_profile); % 'NH' or 'HIx'
end
% End of auditory filtering

if flags.do_internalnoise
    % Internal noise related to absolute threshold of hearing
    N_samples = size(outsig,1);
    int_noise = relanoiborra2019_internalnoise(N_samples, fc);    
    outsig    = outsig  +int_noise;
end

%% 'Haircell' envelope extraction
if flags.do_ihc
    outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
end

% For use with DRNL, this cannot be turned off:
outsig   = outsig * 10^(50/20); % Gain to compensate for the Outer/middle ear attenuation

%% Expansion (Auditory-nerve inspired) and Non-linear adaptation loops:
if flags.do_adt || flags.do_mfb
    % Expansion - Motivation, from Jepsen et al. 2008: '' However, a squaring 
    %   expansion was introduced in the model after hair-cell transduction, reflecting 
    %   the square-law behavior of rate-versus-level functions of the neural 
    %   response in the AN (Yates et al., 1990; Muller et al., 1991)''.
    outsig = outsig.^2;
    
    % non-linear adaptation loops
    outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);
end

if flags.do_mfb
    %% Modulation processing:
    % modfilterbank.m is equivalent to 'relanoiborra2019_mfbtd.m' and 'mfbtdpp.m' togheter:
    [outsig,mfc] = local_modfilterbank(outsig,fs,fc,'argimport',flags,keyvals);
    % [outsig,mfc] = modfilterbank(ousig,fs,fc,'argimport',flags,keyvals);
end

function [outsig,mfc] = local_modfilterbank(insig,fs,fc,varargin)
%MODFILTERBANK  Modulation filter bank
%   Usage: [outsig, mfc] = modfilterbank(insig,fs,fc);
%
%   Input parameters:
%      insig  : Input signal(s)
%      fs     : Sampling rate in Hz,
%      fc     : Center frequencies of the input signals
%
%   Output parameters:
%      outsig : Modulation filtered signals
%      mfc    : Center frequencies of the modulation filters.
%
%   `modfilterbank(insig,fs,fc)` applies a modulation filterbank to the input
%   signals *insig* which are sampled with a frequency of *fs* Hz. Each column in
%   *insig* is assumed to be bandpass filtered with a center frequency stored in *fc*.
%
%   By default, the modulation filters will have center frequencies
%   $0,5,10,16.6,27.77,...$ where each next center frequency is 5/3 times the
%   previous one. For modulation frequencies below (and including) 10 Hz,
%   the real value of the filters are returned, and for higher modulation
%   center frequencies, the absolute value (the envelope) is returned.
%  
%   References: fassel1993modulation dau1997mapI verhey1999 kohlrausch2000
%               jepsen2008cmh
%
%   See also: breebaart2001_preproc
  
% AUTHOR: Stephan Ewert
%
% Modifications by Morten L. Jepsen, Peter L. Søndergaard, and Alejandro Osses.
% Piotr Majdak (21.4.2021) adapted to be a local file

definput.keyvals.mfc=[];
[flags,kv]=ltfatarghelper({},definput,varargin);

nfreqchannels=length(fc);

Q = 2;
bw = 5;
ex=(1+1/(2*Q))/(1-1/(2*Q));

outsig=cell(nfreqchannels,1);

startmf = 5;

% second order modulation Butterworth lowpass filter with a cut-off frequency 
% of 2.5 Hz.
[b_lowpass,a_lowpass] = butter(2,2.5/(fs/2));

% first order modulation Butterworth lowpass filter with a cut-off
% frequency of 150 Hz. This is to remove all modulation frequencies
% above 150 Hz. The motivation behind this filter can be found in kohlrausch2000
if flags.do_LP_150_Hz 
    [b_highest,a_highest] = butter(1,150/(fs/2));
end

% Set the highest modulation frequency as proportion of the corresponding
% center frequency. (see Verhey1999)
if flags.do_mfc_upper_limit
    umf = min(fc.*0.25, 1000);  
end
if flags.do_no_mfc_upper_limit
    umf = 1000*ones(1,length(fc));
end
if flags.do_att_factor
    Factor = 1/sqrt(2); % This is according to Jepsen2008, page 426
end
if flags.do_no_att_factor
    Factor = 1;
end

for freqchannel=1:nfreqchannels
  
  if flags.do_no_LP_150_Hz % this should be the default for dau1997a,b
    outtmp = insig(:,freqchannel);
  end
  % Cut away highest modulation frequencies
  if flags.do_LP_150_Hz
    outtmp = filter(b_highest,a_highest,insig(:,freqchannel));
  end
  
  if umf(freqchannel)==0
    % ----------- only lowpass ---------------------
    outsigblock = filter(b_lowpass,a_lowpass,outtmp);
    mfc = 0;

  else                
    tmp = fix((min(umf(freqchannel),10) - startmf)/bw);
    tmp = 0:tmp;
    mfc = startmf + 5*tmp;
    tmp2 = (mfc(end)+bw/2)/(1-1/(2*Q));
    tmp = fix(log(umf(freqchannel)/tmp2)/log(ex));
    tmp = 0:tmp;
    tmp = ex.^tmp;
    mfc=[0 mfc tmp2*tmp];

    % --------- lowpass and modulation filter(s) ---
    outsigblock = zeros(length(insig),length(mfc));
    outsigblock(:,1) = filter(b_lowpass,a_lowpass,outtmp);

    for nmfc=2:length(mfc)
      w0 = 2*pi*mfc(nmfc)/fs;
      if mfc(nmfc) < 10
   	[b3,a3] = efilt(w0,2*pi*bw/fs);
      else
        [b3,a3] = efilt(w0,w0/Q);
      end
      
      outsigblock(:,nmfc) = 2*filter(b3,a3,outtmp);
    end
    
  end
  
  %% ------------ post-processing --------------------
  
  for nmfc=1:length(mfc) % v2 MJ 17. oct 2006
    if mfc(nmfc) <= 10
      outsigblock(:,nmfc) = 1*real(outsigblock(:,nmfc));
    else
      outsigblock(:,nmfc) = Factor*abs(outsigblock(:,nmfc));
    end
  end
  

  outsig{freqchannel}=outsigblock;
end;


% complex frequency shifted first order lowpass
function [b,a] = efilt(w0,bw);

e0 = exp(-bw/2);

b = 1 - e0;
a = [1, -e0*exp(1i*w0)];


