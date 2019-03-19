function z = comp_fftanalytic(s)
%COMP_FFTANALYTIC Compute analytic representation
%
%   Usage: z = comp_fftanalytic(s);
%
%   COMP_FFTANALYTIC(s) computes the analytic representation of s.  
%   The analytic representation is computed through the FFT of f.
%
%
%   Url: http://ltfat.github.io/doc/comp/comp_fftanalytic.html

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

% AUTHOR: Jordy van Velthoven

Ls = size(s,1);
H = floor(Ls/2);

z = fft(s);
z(2:Ls-H,:) = 2*z(2:Ls-H,:);
z(H+2:Ls,:) = 0;
z = ifft(z);




