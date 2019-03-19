function c=ref_rdftii(f)
%REF_RDFTII  Reference Real DFT type II
%   Usage:  c=ref_rdftii(f);
%
%   Compute RDFTII by explicit formulas.
%
%   The transform is orthonormal
%
%   Url: http://ltfat.github.io/doc/reference/ref_rdftii.html

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
Lhalf=ceil(L/2);

F=zeros(L);

F(:,1)=ones(L,1);

l=(0:L-1).';
for m=1:Lhalf-1
  F(:,2*m)=sqrt(2)*cos(2*pi*m*(l+.5)/L);

  F(:,2*m+1)=sqrt(2)*sin(2*pi*m*(l+.5)/L);

end;

if mod(L,2)==0
  F(:,L)=cos(pi*l);
end;

F=F/sqrt(L);

% dot-transpose will work because F is real.
c=F.'*f;




