function c=ref_rdgt_1(f,g,a,M)
%REF_RDGT_1  Reference Real DGT by fac. and RDFT
%   Usage:  c=ref_rdgt_1(f,g,a,M);
%
%   Compute the factorization and use RDFT
%
%   Url: http://ltfat.github.io/doc/reference/ref_rdgt_1.html

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
W=size(f,2);
R=size(g,2);

N=L/a;

gf=comp_wfac(g,a,M);

% Compute the window application and the DFT modulation.
c=zeros(M*N*R,W);
c(:)=ref_rdft(comp_dgt_fw(f,gf,L,a,M));



