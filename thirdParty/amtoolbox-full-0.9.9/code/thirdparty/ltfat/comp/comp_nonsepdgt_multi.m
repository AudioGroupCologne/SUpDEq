function c=comp_nonsepdgt_multi(f,g,a,M,lt)
%COMP_NONSEPDGT_MULTI  Compute Non-separable Discrete Gabor transform
%   Usage:  c=comp_nonsepdgt_multi(f,g,a,M,lt);
%
%   This is a computational subroutine, do not call it directly.
%
%   Url: http://ltfat.github.io/doc/comp/comp_nonsepdgt_multi.html

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

%   AUTHOR : Nicki Holighaus and Peter L. Søndergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: REF_NONSEPDGT

% Assert correct input.

L=size(f,1);
W=size(f,2);
N=L/a;

% ----- algorithm starts here, split into sub-lattices ---------------

c=zeros(M,N,W,assert_classname(f,g));

mwin=comp_nonsepwin2multi(g,a,M,lt,L);

% simple algorithm: split into sublattices

for ii=0:lt(2)-1
    c(:,ii+1:lt(2):end,:)=comp_dgt(f,mwin(:,ii+1),lt(2)*a,M,[0 1],0,0,0);
end;

% Phase factor correction 
E = zeros(1,N,assert_classname(f,g));
for win=0:lt(2)-1
    for n=0:N/lt(2)-1
        E(win+n*lt(2)+1) = exp(-2*pi*i*a*n*rem(win*lt(1),lt(2))/M);
    end;
end;

c=bsxfun(@times,c,E);

