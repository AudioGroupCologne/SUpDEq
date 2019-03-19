function y = sig_hassager2016(x,B,fs)
%sig_hassager2016 - spectral smoothing algorithm from Hassager et al. (2016)
%   Usage:    y = sig_hassager2016(x,B,fs)
%
%   Input parameters:
%     x   : binaural impulse response(s); SOFA object or matrix with time 
%           as first matrix dimension.
%     B   : bandwidth factor. B=1 represents 1 ERB
%     fs  : sampling rate; required only of x is not a SOFA object
%
%   Output parameters:
%     y   : smoothed version of x.
%
%   SIG_HASSAGER2016(...) is an algortihm for spectral
%   smoothing based on gammatone filters with a variable bandwidth B.
%
%   References:
%     H. G. Hassager, F. Gran, and T. Dau. The role of spectral detail in the
%     binaural transfer function on perceived externalization in a
%     reverberant environment. J. Acoust. Soc. Am., 139(5):2992-3000, 2016.
%     [1]arXiv | [2]www: ]
%     
%     References
%     
%     1. http://arxiv.org/abs/http://dx.doi.org/10.1121/1.4950847
%     2. http://dx.doi.org/10.1121/1.4950847
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_hassager2016.php

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

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

if isstruct(x) % input: SOFA object
  Obj = x;
  fs = Obj.Data.SamplingRate;
  x = shiftdim(Obj.Data.IR,2);
end

[~,excessPhaseTF] = RB_minphase(x,1,'freq');
Nfft = 2.^nextpow2(size(x,1));
X = fftreal(x,Nfft);
f = linspace(0,fs/2,Nfft/2+1);
b = B*1.149*24.7*(4.37*f/1000+1);
Yabs = zeros(size(X));
for idf = 1:length(f)
  Habs2 = abs((1./(1+1i*(f-f(idf))/b(idf))).^4).^2;
  Habs2 = repmat(Habs2(:),[1,size(X,2),size(X,3)]);
  Yabs(idf,:,:) = sqrt(sum(abs(X).^2.*Habs2)./sum(Habs2));
end
% YrandPhase = Yabs.*exp(1i*angle(X));
yRandPhase = ifftreal(Yabs,Nfft);
YminPhase = RB_minphase(yRandPhase,1,'freq');
y = ifft(YminPhase.*excessPhaseTF,Nfft);

if exist('Obj','var') % output: SOFA object
  Obj.Data.IR = shiftdim(y,1);
  y = Obj;
end

end

function [minPhase,excessPhase] = RB_minphase(IR,dim,TFdomainFlag)
% RB_minphase - create minimum-phase filter via causal cepstrum
%
% Usage: [minPhase,excessPhase] = RB_minphase(IR,dim,TFdomainFlag)

% RB, 2016/6/3

Nfft = 2.^nextpow2(size(IR,dim));

TF = fft(IR,Nfft,dim);
logTF = log(abs(TF)+eps);

cep = ifft(logTF,Nfft,dim);
Nshift = mod(dim-1,ndims(cep));
cep1 = shiftdim(cep,Nshift);
cep1(Nfft/2+2:Nfft,:,:,:,:) = 0;    % set non-causal part to zero and 
cep1(2:Nfft/2,:,:,:,:) = 2*cep1(2:Nfft/2,:,:,:,:);    % multiply causal part by 2 (due to symmetry)
cepMinPhase = shiftdim(cep1,ndims(cep)-Nshift);

logTFminPhase = fft(cepMinPhase,Nfft,dim);
TFminPhase = exp(logTFminPhase);

switch TFdomainFlag
  case 'freq'
    minPhase = TFminPhase;
  case 'time'
    minPhase = ifft(TFminPhase,Nfft,dim);
end


if nargout == 2
  switch TFdomainFlag
    case 'freq'
      excessPhase = TF./TFminPhase;
    case 'time'
      excessPhase = ifft(TF./TFminPhase,Nfft,dim);
  end
end

end
