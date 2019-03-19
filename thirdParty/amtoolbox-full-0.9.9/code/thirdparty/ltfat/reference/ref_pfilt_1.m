function h=ref_pfilt_1(f,g,a)
%REF_PFILT_1  Reference PFILT implementation by FFT
%
%   This is the old reference pfilt from before the struct filters where
%   introduced.
%
%   Url: http://ltfat.github.io/doc/reference/ref_pfilt_1.html

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

[L W]=size(f);

g=fir2long(g,L);

% Force FFT along dimension 1, since we have permuted the dimensions
% manually
if isreal(f) && isreal(g)
  h=ifftreal(fftreal(f,L,1).*repmat(fftreal(g,L,1),1,W),L,1);
else
  h=ifft(fft(f,L,1).*repmat(fft(g,L,1),1,W),L,1);
end;

h=h(1:a:end,:);



