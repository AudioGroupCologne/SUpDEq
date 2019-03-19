function cout=ref_col2diag_1(cin);
%REF_COL2DIAG_1  Compute matrix represenation from spreading symbol
%
%  This function is its own inverse.
%
%   Url: http://ltfat.github.io/doc/reference/ref_col2diag_1.html

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
  
L=size(cin,1);
cout=zeros(L);

for jj=0:L-1
  for ii=0:jj-1
    cout(ii+1,jj+1)=cin(ii+1,ii-jj+L+1);
  end;
    for ii=jj:L-1
    cout(ii+1,jj+1)=cin(ii+1,ii-jj+1);
  end;
end;

% The second code also works.
if 0

  for ii=0:L-1
    for jj=0:ii
      cout(ii+1,jj+1)=cin(ii+1,ii-jj+1);
    end;
    for jj=ii+1:L-1
      cout(ii+1,jj+1)=cin(ii+1,ii-jj+L+1);
    end;  
  end;
end;



