function f=involute(f,dim);
%INVOLUTE  Involution 
%   Usage: finv=involute(f);
%          finv=involute(f,dim);
%
%   INVOLUTE(f) will return the involution of f.
%
%   INVOLUTE(f,dim) will return the involution of f along dimension dim.
%   This can for instance be used to calculate the 2D involution:
%
%     f=involute(f,1);
%     f=involute(f,2);
%
%   The involution finv of f is given by:
%
%     finv(l+1)=conj(f(mod(-l,L)+1));
%
%   for l=0,...,L-1.
%
%   The relation between conjugation, Fourier transformation and involution
%   is expressed by:
%
%     conj(dft(f)) == dft(involute(f))
%
%   for all signals f. The inverse discrete Fourier transform can be
%   expressed by:
%
%     idft(f) == conj(involute(dft(f)));
%
%   See also:  dft, pconv
%
%   Url: http://ltfat.github.io/doc/fourier/involute.html

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

% Assert correct input.

%   AUTHOR : Peter L. Søndergaard
%   TESTING: TEST_INVOLUTE
%   REFERENCE: OK

error(nargchk(1,2,nargin));

if nargin==1
  dim=[];
end;

L=[];
[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'INVOLUTE');

% This is where the calculation is performed.
% The reshape(...,size(f) ensures that f will keep its
% original shape if it is multidimensional.
f=reshape(conj([f(1,:); ...
	  flipud(f(2:L,:))]),size(f));

f=assert_sigreshape_post(f,dim,permutedsize,order);

