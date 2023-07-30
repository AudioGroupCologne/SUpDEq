function [outsig,mfc,params] = modfilterbank(insig,fs,fc,varargin)
%MODFILTERBANK  Modulation filter bank
%   Usage: [outsig, mfc,params] = modfilterbank(insig,fs,fc);
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
%   MODFILTERBANK(insig,fs,fc) applies a modulation filterbank to the input
%   signals insig which are sampled with a frequency of fs Hz. Each column in
%   insig is assumed to be bandpass filtered with a center frequency stored in fc.
%
%   By default, the modulation filters will have center frequencies
%   0,5,10,16.6,27.77,... where each next center frequency is 5/3 times the
%   previous one. For modulation frequencies below (and including) 10 Hz,
%   the real value of the filters are returned, and for higher modulation
%   center frequencies, the absolute value (the envelope) is returned.
%  
%   References:
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%     R. Fassel and D. Pueschel. Modulation detection and masking using
%     deterministic and random maskers. Contributions to Psychological
%     Acoustics, edited by A. Schick (Universitaetsgesellschaft Oldenburg,
%     Oldenburg), pages 419--429, 1993.
%     
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am.,
%     124(1):422--438, 2008.
%     
%     A. Kohlrausch, R. Fassel, and D. Torsten. The influence of carrier
%     level and frequency on modulation and beat-detection thresholds for
%     sinusoidal carriers. J. Acoust. Soc. Am., 108:723--734, 2000.
%     
%     J. Verhey, T. Dau, and B. Kollmeier. Within-channel cues in
%     comodulation masking release (cmr): experiments and model predictions
%     using a modulation-filterbank model. J. Acoust. Soc. Am.,
%     106:2733--2745, 1999.
%     
%
%   See also: breebaart2001
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/modfilterbank.php


%   #Author: Stephan Ewert (1999-2004) and Morten L. Jepsen: Original version
%   #Author: Peter L. Søndergaard (2009-2013): adapted to AMT
%   #Author: Piotr Majdak (2013-): adapted to AMT 1.0
%   #Author: Alejandro Osses (2020): extended with LP_150_Hz

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.keyvals.mfc=[];
definput.import = {'modfilterbank'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if nargout >= 3
    bStore_params = 1;
else
    bStore_params = 0;
end

nfreqchannels=length(fc);
if nfreqchannels == 0
    % This can be the case if the filter bank stage was by-passed, then:
    nfreqchannels = size(insig,2);
    fc = fs/2; % frequency that is high enough to ensure that all modfilters are assessed
end
Q = kv.Q_mfb; % Q = 2 is the default for this function
bw = 5;
ex=(1+1/(2*Q))/(1-1/(2*Q));

outsig=cell(nfreqchannels,1);

startmf = 5;

% second order modulation Butterworth lowpass filter with a cut-off frequency 
% of 2.5 Hz.
[b_lowpass,a_lowpass] = butter(2,2.5/(fs/2));

if bStore_params
    params.mfb_b(1,:) = b_lowpass;
    params.mfb_a(1,:) = a_lowpass;
    
    params.fs = fs;
    params.description = ['Coefficients for all modulation filters used in ',mfilename,' obtained at fs=',num2str(fs)];
    amt_disp(params.description);
end
        
% first order modulation Butterworth lowpass filter with a cut-off
% frequency of 150 Hz. This is to remove all modulation frequencies
% above 150 Hz. The motivation behind this filter can be found in kohlrausch2000
if flags.do_LP_150_Hz || flags.do_LP_150_Hz_att
    [b_lp_150_Hz,a_lp_150_Hz] = butter(1,150/(fs/2));
    
    if bStore_params
        params.b_lp_150_Hz = b_lp_150_Hz;
        params.a_lp_150_Hz = a_lp_150_Hz;
    end
end

% Set the highest modulation frequency as proportion of the corresponding
% center frequency. (see Verhey1999)
mfc_upper_limit_max = kv.mfc_upper_limit_max; % Hz
if flags.do_mfc_upper_limit
    umf = min(fc.*0.25, mfc_upper_limit_max);
end
if flags.do_no_mfc_upper_limit
    umf = mfc_upper_limit_max*ones(1,nfreqchannels);
end
if flags.do_att_factor
    Factor = 1/sqrt(2); % This is according to Jepsen et al. (2008), page 426
end
if flags.do_no_att_factor
    Factor = 1;
end

if flags.do_LP_150_Hz_att
    K =  8192; % arbitrary number
    f = (0:K-1)/K * fs/2;
    hlpf    = freqz(b_lp_150_Hz,a_lp_150_Hz, K); % frequency response of the 150 Hz filter
    hlpf_dB = 20*log10(abs(hlpf));
end

for freqchannel=1:nfreqchannels
  
  if flags.do_no_LP_150_Hz || flags.do_LP_150_Hz_att % this should be the default for dau1997a,b
    outtmp = insig(:,freqchannel);
  end
  % Cut away highest modulation frequencies
  if flags.do_LP_150_Hz
    outtmp = filter(b_lp_150_Hz,a_lp_150_Hz,insig(:,freqchannel));
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

    if flags.do_LP_150_Hz_att
        mfc_gains = interp1(f(:)',hlpf_dB(:)',mfc);
    end
    
    for nmfc=2:length(mfc)
        w0 = 2*pi*mfc(nmfc)/fs;
        if mfc(nmfc) < 10
            [b3,a3] = efilt(w0,2*pi*bw/fs);
        else
            [b3,a3] = efilt(w0,w0/Q);
        end
      
        if bStore_params
            params.mfb_b(nmfc,1:length(b3)) = b3;
            params.mfb_a(nmfc,1:length(a3)) = a3;
        end
        
        outsigblock(:,nmfc) = 2*filter(b3,a3,outtmp); % emphasis of 6 dB
        
        if flags.do_LP_150_Hz_att
            outsigblock(:,nmfc) = gaindb( outsigblock(:,nmfc), mfc_gains(nmfc));
        end
    end
    
  end
  
  %% ------------ post-processing --------------------
  % If enabled, the phase information for modulation filters with mfc below
  %   or equal to 10 Hz is kept, and above 10 Hz the envelope is assessed.
  if flags.do_phase_insens_hilbert
      for nmfc=1:length(mfc) % v2 MJ 17. oct 2006
        if mfc(nmfc) <= 10
          outsigblock(:,nmfc) = 1*real(outsigblock(:,nmfc));
        else
          outsigblock(:,nmfc) = Factor*abs(outsigblock(:,nmfc));
        end
      end
  end
  
  % If disabled, the bandpassed signals are returned:
  if flags.do_no_phase_insens
      outsigblock(:,2:end) = outsigblock(:,2:end)/2; % this factor of 2 is to
      % to comensate for the gain of 2 in Line L149: 'outsigblock(:,nmfc) = 2*filter(b3,a3,outtmp);'
      % which is applied to all modulation bandpass filters (excluding the first
      % LPF)
  end
      
  outsig{freqchannel}=outsigblock;
end;

%% ------------ subfunctions ------------------------


% complex frequency shifted first order lowpass
function [b,a] = efilt(w0,bw);

e0 = exp(-bw/2);

b = 1 - e0;
a = [1, -e0*exp(1i*w0)];


