function [b,a]=gammachirp(fc,fs,varargin)
%GAMMACHIRP  Gammachirp filter coefficients
%   Usage: [b,a] = gammachirp(fc,fs,n,betamul);
%          [b,a] = gammachirp(fc,fs,n);
%          [b,a] = gammachirp(fc,fs);
%
%   Input parameters:
%      fc    : center frequency in Hz.
%      fs    : sampling rate in Hz.
%
%   Output parameters:
%      b     :  nominator coefficients.
%      a     :  denominator coefficients.
%
%   gammachirp takes the following key-value pairs:
%
%     'order',n       filter order (order of Gamma function t^(OrderG-1) )
%
%     'beta',b        bandwidth of the filter (exp(-2*pi*CoefERBw*ERB(f)))
%
%     'c',c           c-coefficient exp(j*2*pi*Frs + CoefC*ln(t))
%
%     'phase',phase   initial phase (0 ~ 2*pi)
%
%   gammachirp takes the following flags:
%
%     'carrier'    Carrier (cos,sin,complex,envelope: 3 letters)
%
%     'norm'       Normalization of peak spectrum level (no, peak)
% 
%
%
%   GAMMACHIRP(fc,fs,n,betamul) computes the filter coefficients of a
%   digital gammachirp filter with center frequency fc, order n, sampling
%   rate fs and bandwith determined by betamul. The bandwidth beta of
%   each filter is determined as betamul times audfiltbw of the center
%   frequency of corresponding filter.
%
%   By default, the returned filter coefficients comes from the all-pole
%   approximation described in Lyon (1997). The filters are normalized to
%   have a 0 dB attenuation at the center frequency (another way of
%   stating this is that their impulse responses will have unit area).
%
%   GAMMACHIRP(fc,fs) will do as above for a 4th order filter.
%
%   If fc is a vector, each entry of fc is considered as one center
%   frequency, and the corresponding coefficients are returned as row
%   vectors in the output.
%
%   The impulse response of the gammachirp filter is given by:
%
%        g(t) = a*t^(n-1)*cos(2*pi*fc*t)*exp(-2*pi*beta*t)
%
%  
%
%   To create the filter coefficients of a 1-erb spaced filter bank using
%   gammachirp filters use the following construction:
%
%     [b,a] = gammachirp(erbspacebw(flow,fhigh),fs,'complex');
%
%   To apply the (complex valued) filters to an input signal, use
%   FILTERBANKZ:
%
%     outsig = 2*real(ufilterbankz(b,a,insig));
%  
%   References:
%     T. Irino and R. D. Pattersion. A time-domain, level-dependent auditory
%     filter: The gammachirp. J. Acoust. Soc. Am., 101(412), 1997.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/gammachirp.php


%   #Author : Toshio Irino (2006)
%   #Author: Clara Hollomey (2021): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% ------ Checking of input parameters ---------
  
if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(fc) || ~isvector(fc) || any(fc<0) || any(fc>fs/2)
  error(['%s: fc must be a vector of positive values that are less than half ' ...
         'the sampling rate.'],upper(mfilename));
end;

definput.keyvals.n=4;
definput.keyvals.CoefERBw=1.019; % Default GammaTone value
ERB = ones(fc,1);
definput.keyvals.c=1;
definput.keyvals.phase=0;
definput.flags.carrier={'cos','sin','complex','envelope'};
definput.flags.norm={'no', 'peak'};

[flags,kv,n, c]  = ltfatarghelper({'n','c'},definput,varargin);

if ~isnumeric(n) || ~isscalar(n) || n<=0 || fix(n)~=n
  error('%s: n must be a positive, integer scalar.',upper(mfilename));
end

  b = exp(-2*pi*kv.CoefERBw*ERB(fc));

ERBrate = fc2erb(fc);
ERBw = f2erb(fc);

LenGC1kHz = (40*max(n)/max(kv.CoefERBw) + 200)*fs/16000;  % 2 Aug 96 
ERBw1kHz = f2erb(1000);	

if flags.do_sin kv.phase = kv.phase - pi/2*ones(1,length(fc)); end;
%%% Phase compensation
phase = kv.phase + c.*log(fc/1000); % relative phase to 1kHz

LenGC = fix(LenGC1kHz*ERBw1kHz./ERBw);

%%%%%  Production of GammaChirp  %%%%%
GC       = zeros(length(fc),max(LenGC));
if nargout > 2
 ERBwfc = f2erb(fc);
 fpeak = fc + c.*ERBwfc.*kv.CoefERBw./n; 
end

if nargout > 3, InstFreq = zeros(length(fc),max(LenGC));        end


for nch = 1:length(fc)
  
  t = (1:LenGC(nch)-1)/fs;

  GammaEnv = t.^(n(nch)-1).*exp(-2*pi*kv.CoefERBw(nch)*ERBw(nch)*t);
  GammaEnv = [ 0 GammaEnv/max(GammaEnv)];

  if flags.do_envelope
    carrier = ones(size(GammaEnv));
  elseif flags.do_complex
    carrier = [ 0 exp(1i * (2*pi*fc(nch)*t + c(nch)*log(t) +phase(nch)) )];
  else
    carrier = [ 0 cos(2*pi*fc(nch)*t + c(nch)*log(t) +phase(nch))];
  end;

  GC(nch,1:LenGC(nch)) = GammaEnv.*carrier;

  if nargout > 3, 
    InstFreq(nch,1:LenGC(nch)) = [0, [fc(nch) + c(nch)./(2*pi*t)]];
  end
  
  if flags.do_peak  % peak gain normalization
     [frsp, freq] = freqz(GC(nch,1:LenGC(nch)),1,LenGC(nch),SR);
     ERBwp = f2erb(fc(nch));
     fp = fc(nch) + c.*ERBwp.*kv.CoefERBw(nch)./n(nch);
     
     [~, np] = min(abs(freq-fp));
     GC(nch,:) = GC(nch,:)/abs(frsp(np));
  end;
  b = GC(nch, 1:LenGC(nch));
  a = 1;
end 


