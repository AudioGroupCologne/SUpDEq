function sig = sig_hofman1998(varargin)
%SIG_HOFMAN1998 - Stimuli from Hofman and Opstal (1998) 
%   Usage: sig = sig_hofman1998(nprime,f0,Obj,azi,dur,condition)
%   
%
%   Input parameters:
%     nprime : Highest altered harmonic. Default is 0.
%     f0     : Fundamental frequency. Default is 125 Hz.
%     Obj    : HRTFs as SOFA object.
%     azi    : Azimuth. Default is -37 deg.
%     dur    : Duration in seconds. Default is 0.1 s.
%     fs     : Sampling rate in Hz. Default is 48 kHz.
%
%
%   Output parameters:
%     sig    : signal wave form
%
%   SIG_HOFMAN1998 generates three types of stimuli having different 
%   spectro-temporal characteristics but identical long-term power spectra: 
%   noise bursts, frequency-modulated tones, and trains of short noise bursts.
%
%   The additional optional parameter condition flag may be one of:
%
%     'ILD'   ILDs up to nprime set to zero. This is the default.
%     'ISLD'  ISLDs (interaural spectral level differeces) maintained while
%             flattening right-ear HRTF up to nprime.
%
%   References:
%     P. M. Hofman and A. J. Van Opstal. Spectro-temporal factors in
%     two-dimensional human sound localization. J. Acoust. Soc. Am.,
%     103(2634), 1998.
%     
%
%   See also: sig_schroeder1970
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_hofman1998.php


%   #Author: Robert Baumgartner (2021), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.keyvals.nprime = 0;
definput.keyvals.f0 = 125; % Hz
definput.keyvals.Obj = [];
definput.keyvals.dur = 1;
definput.keyvals.fs = 48e3;
definput.keyvals.azi = -37; % deg (right-hand side)
definput.flags.expirement = {'ILD','ISLD'};

[flags,kv]=ltfatarghelper({'nprime','f0','Obj','azi','dur'},definput,varargin);

%% Fixed settings
N = 38; % # harmonics

%% Amplitude and phase data from Appendix A
Xabs = [...
 0.6047,0.3739,0.3056,0.2865,0.3104,0.6210,0.3676,0.2397,0.1573,0.6016,...
 0.3057,0.1411,0.0494,0.0351,0.0231,0.0554,0.0767,0.0688,0.0619,0.0645,...
 0.0791,0.1389,0.1545,0.1345,0.1296,0.1150,0.0754,0.0624,0.0733,0.0209,...
 0.0363,0.0476,0.0352,0.0159,0.0011,0.0028,0.0106,0.0197];
ph0 = [...
 0.000,6.202,0.098,2.290,0.033,4.097,2.487,2.583,5.533,0.604,...
 3.471,3.962,3.471,4.863,1.693,1.884,1.215,3.804,0.276,2.377,...
 2.888,5.156,4.635,1.311,1.809,2.188,5.119,5.482,3.679,6.115,...
 2.105,3.814,0.889,0.611,2.697,4.510,5.335,1.330];

%% HRTFs
if not(isempty(kv.Obj))
  kv.fs = kv.Obj.Data.SamplingRate;
  Nfft = round(kv.fs/kv.f0);
  ir = SOFAspat([1;zeros(Nfft-1,1)],kv.Obj,kv.azi,0);
  H = fftreal(ir,Nfft);
  Habs = abs(H(2:N+1,:));
  Hang = angle(H(2:N+1,:));
else
  Habs = zeros(N,2);
  Hang = zeros(N,2);
end

%% Alterations
if kv.nprime > 0
  if flags.do_ILD
    Habs(1:kv.nprime,:) = repmat(mean(Habs(1:kv.nprime,:),2),[1,2]);
  elseif flags.do_ISLD
    HabsR = Habs(:,2);
    HabsL = Habs(:,1);
    ISLD = HabsR./HabsL;
    HabsR(1:kv.nprime) = mean(HabsR(1:kv.nprime));
    HabsL = HabsR./ISLD;
    Habs = [HabsL,HabsR];
  end
end

%% Signal synthesis acc. to Eq. (1)
t = linspace(0,kv.dur,kv.fs*kv.dur)';
x = zeros(length(t),2);
for n = 1:N
  for ch = 1:2
    x(:,ch) = x(:,ch) + (Xabs(n)*Habs(n,ch)) * cos(n*2*pi*kv.f0*t + ph0(n)+Hang(n,ch));
  end
end

%% Fade in/out
lenRamp = 2*kv.fs*kv.dur/10;
ramp = hann(lenRamp);
win = [ramp(1:lenRamp/2);ones(length(x)-lenRamp,1);ramp(lenRamp/2+1:end)];
x = x.*repmat(win,[1,2]);

%% Normalization
sig = x/max(x(:));
end


