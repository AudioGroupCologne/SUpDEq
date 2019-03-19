function c=ref_fftreal(f)
%REF_FFTREAL  Reference FFTREAL
%
%  FFTREAL is computed by doing an FFT, and keeping only the positive
%  frequencies.
%
%   Url: http://ltfat.github.io/doc/reference/ref_fftreal.html

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
  

L=size(f,1);
c=fft(f);

c=c(1:floor(L/2)+1,:);



