function c = phaseunlockreal(c,a,M,varargin)
%PHASEUNLOCKREAL  Undo phase lock of Gabor coefficients
%   Usage:  c=phaseunlockreal(c,a,M);
%
%   PHASEUNLOCKREAL(c,a,M) removes phase locking from the Gabor coefficients c.
%   The coefficient must have been obtained from a DGTREAL with parameter a*
%   and using the 'timeinv' flag.
%
%   Phase locking the coefficients modifies them so as if they were obtained
%   from a frequency-invariant Gabor system. A filter bank produces phase locked
%   coefficients. 
%
%   See also: dgt, phaselockreal
%
%   References:
%     M. Puckette. Phase-locked vocoder. Applications of Signal Processing to
%     Audio and Acoustics, 1995., IEEE ASSP Workshop on, pages 222 --225,
%     1995.
%     
%
%   Url: http://ltfat.github.io/doc/gabor/phaseunlockreal.html

% Copyright (C) 2005-2016 Peter L. Soendergaard <peter@sonderport.dk>.
% This file is part of LTFAT version 2.2.0
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

%   AUTHOR:    Peter Balazs, Peter L. Søndergaard, Zdenek Prusa
%   TESTING:   OK
%   REFERENCE: OK

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if  ~isscalar(a) || ~isnumeric(a) || rem(a,1)~=0
  error('a must be integer');
end;

if  ~isscalar(M) || ~isnumeric(M) || rem(M,1)~=0
  error('M must be integer');
end;

definput.flags.accuracy={'normal', 'precise'};
flags=ltfatarghelper({},definput,varargin);

M2=size(c,1);
M2user = floor(M/2) + 1;

if M2~=M2user
    error('%s: Size of s does not comply with M.',upper(mfilename));
end

N=size(c,2);

if flags.do_normal
    TimeInd = (0:(N-1))*a;
    FreqInd = (0:(M2-1));

    phase = FreqInd'*TimeInd;
    phase = mod(phase,M);
    phase = exp(-2*1i*pi*phase/M);

    % Handle multisignals
    c=bsxfun(@times,c,phase);
elseif flags.do_precise
    
    frames = ifftreal(c,M);
    for n=0:N-1
        frames(:,n+1,:) = circshift(frames(:,n+1,:),n*a);
    end
    c = fftreal(frames);
    
end

