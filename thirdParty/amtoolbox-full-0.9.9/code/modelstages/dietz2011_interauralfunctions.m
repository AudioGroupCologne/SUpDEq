function [outp] = dietz2011_interauralfunctions(insig,fs,fc,varargin)
%dietz2011_interauralfunctions  Interaural stages of Dietz 2011
%
%   Input parameters
%     insig  : input signal
%     fs     : sampling frequencies
%     fc     : center frequencies
%
%   Output parameters:
%     outp   : Structure containing the output. See the description
%              below.
%
%   Calculate interaural parameters like interaural transfer function (ITF),
%   interaural phase difference (IPD), interaural time difference (ITD) and
%   interaural coherence (IC) for the given input signals. The ITD is calculated
%   by dividing the IPD by the instantaneous frequency of the signal. In
%   addition lowpassed filtered version of the interaural parameters are
%   calculated (ending with _lp) to simulate a finite time resolution of the
%   binaural system.
%
%   The output structure outp contains the following fields:
%     itf       : transfer function
%     ipd       : phase difference in rad
%     itd       : interaural time difference based on instantaneous frequency
%     itd_C     : interaural time difference based on center frequency
%     f_inst_1  : instantaneous frequencies of left ear signal
%     f_inst_2  : instantaneous frequencies of right ear canal signal
%     f_inst    : instantaneous frequencies (average of f_inst1 and 2)
%     ic        : interaural coherence
%     rms       : rms value of frequency channels for weighting
%     ild_lp    : based on low passed-filtered insig, level difference in dB
%     ipd_lp    : based on lowpass-filtered itf, phase difference in rad
%     itd_lp    : based on lowpass-filtered itf, interaural time difference
%     itd_C_lp  : based on lowpass-filtered itf, interaural time difference
%
%   The _lp values are not returned if the 'nolowpass' flag is set.
%
%   DIETZ2011_INTERAURALFUNCTIONS accepts the following optional parameters:
%
%     'compression_power',cpwr
%                    Applied compression of the signal on the cochlea with
%                    ^compression_power. The default value is 0.4.
%
%     'tau_cycles',tau_cycles
%                    Temporal resolution of binaural processor in terms of
%                    cycles per frequency channel. The default value is 5.
%
%     'signal_level_dB_SPL',signal_level
%                    Sound pressure level of left channel. Used for data
%                    display and analysis. Default value is 70.
%
%     'lowpass'      Calculate the interaural parameters of the lowpassed
%                    signal/ITF (_lp return values). This is the default.
%
%     'nolowpass'    Don't calculate the lowpass based interaural parameters.
%                    The _lp values are not returned.
%
%   See also: dietz2011, dietz2011_filterbank
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592-605, 2011. [1]http ]
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S016763931000097X
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/dietz2011_interauralfunctions.php

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

% AUTHOR: Mathias Dietz, Martin Klein-Hennig (for AMT), Hagen Wierstorf (for AMT)

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig)
  error('%s: insig has to be a numeric signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

if ~isnumeric(fc)
  error('%s: fc has to be a numeric signal!',upper(mfilename));
end

definput.import = {'dietz2011_interauralfunctions'};
[flags,kv]  = ltfatarghelper({},definput,varargin);

% if signal contains no frequency channels return empty
if size(insig,2)==0, outp = []; return; end

% ----- interaural parameters -----
% interaural transfer function (ITF), eq. 2 in Dietz (2011)
outp.itf = insig(:,:,2) .* conj(insig(:,:,1));
% interaural phase difference (IPD), eq. 3 in Dietz (2011) but without low pass
% filtering
outp.ipd = angle(outp.itf);
% instantaneous frequency (f_inst), eq. 4 in Dietz (2011)
for k = 1:length(fc)
  outp.f_inst_left(:,k)  = calc_f_inst(insig(:,k,1),fs);
  outp.f_inst_right(:,k) = calc_f_inst(insig(:,k,2),fs);
end
outp.f_inst = max(eps,0.5*(outp.f_inst_left + outp.f_inst_right)); % to avoid division by zero
% interaural time difference (ITD), based on instantaneous frequencies
outp.itd = 1/(2*pi)*outp.ipd./outp.f_inst;
% interaural time difference (ITD), based on center frequencies <== is this really needed?
for k = 1:length(fc)
    outp.itd_C(:,k) = 1/(2*pi)*outp.ipd(:,k)/fc(k);
end

% ----- low pass versions of interaural parameters -----
% The low pass is simulating a finite time resolution of the binaural system and
% is given by the coh_cycles parameter
tau = kv.tau_cycles./fc;
if flags.do_lowpass
  a = exp( -1./(fs.*tau) );
  % interaural phase difference (IPD) of lowpassed signals, eq. 3 in Dietz (2011)
  outp.ipd_lp = angle(lowpass(outp.itf,a));
  outp.f_inst_lp = lowpass(outp.f_inst,a);
  outp.itd_lp = 1/(2*pi)*outp.ipd_lp./outp.f_inst_lp;
  for k = 1:length(fc)
    outp.itd_C_lp(:,k) = 1/(2*pi)*outp.ipd_lp(:,k)/fc(k);
  end
  % interaural level difference
  inoutsig = lowpass(abs(insig),a);
  inoutsig(abs(inoutsig)<eps) = eps; % avoid division by zero and log(0)
  outp.ild_lp = 20./kv.compression_power.*log10(inoutsig(:,:,2)./inoutsig(:,:,1));
end

% interaural coherence (IC) estimated by interaural vector strength (IVS), see
% eq. 7 in Dietz (2011)
outp.ic = interaural_vector_strength(outp.itf, tau, fs);

% weighting of channels for cumulative ixd determination
% sqrt(2) is due to half-wave rectification (included 28th Sep 07)
outp.rms = kv.signal_level_dB_SPL*kv.compression_power + 20*log10(sqrt(2)*min(rms(squeeze(abs(insig(:,:,1)))),rms(squeeze(abs(insig(:,:,2))))));
outp.rms = max(outp.rms,eps); % avoid negative weights

end


%% ----- internal functions -----
% lowpass
function outsig = lowpass(sig, a)
  % outsig = lowpass(sig, a)
  % This is a simple low-pass filter y(n) = (1-a)*x(n) - a*y(n-1)
  % Meaning of parameter a:
  %   a - damping coefficient (0 - no filtering, 1 - flat output)
  %   tau = 1/(2*pi*f_c)      where f_c is the cutoff frequency of the filter
  %   a = exp(-1/(fs*tau))   where fs - sampling frequency
  %
  if ndims(sig)==3
    sig = permute(sig,[1 3 2]); % => [samples ears channels]
    channels = size(sig,3);
    outsig = zeros(size(sig));
    for ii=1:channels
      outsig(:,:,ii) = filter([1-a(ii)], [1, -a(ii)], sig(:,:,ii));
    end
    outsig = permute(outsig,[1 3 2]); % => [samples channels ears]
  else
    channels = size(sig,2);
    outsig = zeros(size(sig));
    for ii=1:channels
      outsig(:,ii) = filter([1-a(ii)], [1, -a(ii)], sig(:,ii));
    end
  end
end


% Hagen: I have removed the temporal smoothing from this function. The lowpass
% has to be applied afterwards. The positive effect is a little bit smoother ITD
% value. The negative effect is a longer time after the onset at the beginning
% to reach the actual value.
function f_inst = calc_f_inst(sig,fs)
  % function f_inst = calc_f_inst(sig,fs);
  %
  % Calculates instantaneous frequency from a complex (analytical) signal
  % using first order differences
  %
  % input parameters:
  %   sig  : complex (analytical) input signal
  %   fs   : sampling frequency of sig
  %
  % output values:
  %   f_inst:   vector of estimated inst. frequency values
  %
  % copyright: Universitaet Oldenburg
  % author   : volker hohmann
  % date     : 12/2004
  sig = sig./(abs(sig)+eps);
  f_inst = [0; sig(2:end).*conj(sig(1:end-1))];
  f_inst = angle(f_inst)/2/pi*fs;
end


% Calculate the interaural coherence (IC) with the help of the interaural vector
% strength (IVS) as defined by eq. 7 in Dietz (2011)
function ivs = interaural_vector_strength(itf,tau_coherence,fs)
  %tau_coherence = 15e-3; % good value for ipd_fine
  c_coh = exp(-1./(fs.*tau_coherence));
  if length(tau_coherence)==1
    ivs = abs(filter(1-c_coh,[1 -c_coh],itf))./abs(filter(1-c_coh,[1 -c_coh],abs(itf)));
  elseif length(tau_coherence)==size(itf,2)
    ivs = zeros(size(itf));
    for n=1:length(tau_coherence)
      ivs(:,n) = abs(filter(1-c_coh(n),[1 -c_coh(n)],itf(:,n)))./ ...
          abs(filter(1-c_coh(n),[1 -c_coh(n)],abs(itf(:,n))));
    end
  else
    error('wrong number of tau_coherence values')
  end
end

% vim: set sw=2 ts=2 expandtab textwidth=80:

