function [c,Ls]=wmdct2(f,g1,p3,p4,p5)
%WMDCT2  2D Discrete windowed MDCT transform
%   Usage: c=wmdct2(f,g,M); 
%          c=wmdct2(f,g1,g2,[M1,M2]);
%          c=wmdct2(f,g1,g2,[M1,M2],[L1,L2]);
%          [c,L]=wmdct2(f,g1,g2,[M1,M2],[L1,L2]);
%
%   Input parameters:
%         f        : Input data, matrix.
%         g,g1,g2  : Window functions.
%         M,M1,M2  : Number of bands.
%         L1,L2    : Length of transform to do.
%   Output parameters:
%         c        : array of coefficients.
%         Ls       : Original size of input matrix.
%
%   WMDCT2(f,g,M) calculates a two dimensional Modified Discrete Cosine
%   transform of the input signal f using the window g and parameter 
%   M along each dimension.
%
%   For each dimension, the length of the transform will be the smallest
%   possible that is larger than the length of the signal along that
%   dimension. f will be appropriately zero-extended.
%
%   All windows must be whole-point even.
% 
%   WMDCT2(f,g,M,L) computes a 2D windowed MDCT as above, but does a
%   transform of length L along each dimension. f will be cut or
%   zero-extended to length L before the transform is done.
%
%   [c,Ls]=wmdct(f,g,M) or [c,Ls]=wmdct(f,g,M,L) additionally return the
%   length of the input signal f. This is handy for reconstruction.
%
%   c=WMDCT2(f,g1,g2,M) makes it possible to use different windows along
%   the two dimensions.
%
%   The parameters L, M and Ls can also be vectors of length 2. In this case
%   the first element will be used for the first dimension and the second
%   element will be used for the second dimension. 
%
%   The output c has 4 or 5 dimensions. The dimensions index the following
%   properties:
%
%      1. Number of translation along 1st dimension of input.
%
%      2. Number of channel along 1st dimension  of input
%
%      3. Number of translation along 2nd dimension of input.
%
%      4. Number of channel along 2nd dimension  of input
%
%      5. Plane number, corresponds to 3rd dimension of input. 
% 
%   See also:  wmdct, iwmdct2, dgt2, wildual
%
%   Url: http://ltfat.github.io/doc/gabor/wmdct2.html

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

%   AUTHOR : Peter L. Søndergaard.

error(nargchk(3,5,nargin));

L=[];

if prod(size(p3))>2
  % Two windows was specified.
  g2=p3;
  M=p4;
  if nargin==5
    L=p5;
  end;
else
  g2=g1;
  M=p3;
  if nargin==4
    L=p4;
  end;
end;
  
if isempty(L)
  L1=[];
  L2=[];
else
  L1=L(1);
  L2=L(2);
end;

% Expand M if necessary to two elements
M=bsxfun(@times,M,[1 1]);

Ls=size(f);
Ls=Ls(1:2);

c=wmdct(f,g1,M(1),L1);
c=wmdct(c,g2,M(2),L2,'dim',3);

