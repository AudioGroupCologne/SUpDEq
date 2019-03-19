function f=ref_igdgt(coef,g,a,M,c_t,c_f,c_w)
%REF_GDGT  Reference generalized DGT
%   Usage:  f=ref_igdgt(c,g,a,M,c_t,c_f,c_w);
%
%   Linear algebra version of the algorithm. Create big matrix
%   containing all the basis functions and multiply with the transpose.
%
%   Url: http://ltfat.github.io/doc/reference/ref_igdgt.html

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

N=size(coef,1)/M;
L=N*a;

F=zeros(L,M*N);

l=(0:L-1).';

if length(g)<L
  g=fir2long(g,L);
end;

for n=0:N-1	   
  for m=0:M-1
    F(:,M*n+m+1)=exp(2*pi*i*(m+c_f)*(l+c_t)/M).*circshift(g,n*a+c_w);
  end;
end;

f=F*coef;




