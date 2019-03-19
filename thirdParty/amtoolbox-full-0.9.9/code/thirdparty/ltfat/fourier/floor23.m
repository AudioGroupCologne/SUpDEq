function [nfft,tableout]=floor23(n)
%FLOOR23  Previous number with only 2,3 factors
%   Usage: nceil=floor23(n);
%
%   FLOOR23(n) returns the first number less than or equal to n,
%   which can be written as a product of powers of 2 and 3.
%
%   The algorithm will look up the best size in a table, which is computed
%   the first time the function is run. If the input size is larger than the
%   largest value in the table, the input size will be reduced by factors of
%   2, until it is in range.
%
%   [nceil,table]=FLOOR23(n) additionally returns the table used for lookup.
%
%   Examples:
%   ---------
%
%   Return the first number smaller or equal to 26 that can be written
%   solely as products of powers of 2 and 3*:
% 
%     floor23(26)
%
%   This plot shows the behaviour of FLOOR23 and CEIL23 for numbers
%   up to 100:
%
%     x=1:100;
%     plot(x,floor23(x),x,ceil23(x));
%     legend('floor23','ceil23','Location','Northwest');
%
%   See also: ceil23, floor235, nextfastfft
%
%   Url: http://ltfat.github.io/doc/fourier/floor23.html

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
  
%   AUTHOR: Peter L. Søndergaard
  
  
persistent table;
  
maxval=2^20;

if isempty(table)
    % Compute the table for the first time, it is empty.
    l2=log(2);
    l3=log(3);
    l5=log(5);
    lmaxval=log(maxval);
    table=zeros(143,1);
    ii=1;
    prod2=1;
    for i2=0:floor(lmaxval/l2)
        prod3=prod2;
        for i3=0:floor((lmaxval-i2*l2)/l3)               
            table(ii)=prod3; 
            prod3=prod3*3;
            ii=ii+1;
        end;
        prod2=prod2*2;            
    end;
    table=sort(table);
end;

% Copy input to output. This allows us to efficiently work in-place.
nfft=n;

% Handle input of any shape by Fortran indexing.
for ii=1:numel(n)
  n2reduce=0;
  
  if n(ii)>maxval
    % Reduce by factors of 2 to get below maxval
    n2reduce=ceil(log2(nfft(ii)/maxval));
    nfft(ii)=nfft(ii)/2^n2reduce;
  end;
  
  % Use a simple bisection method to find the answer in the table.
  from=1;
  to=numel(table);
  while from<=to
    mid = round((from + to)/2);    
    diff = table(mid)-nfft(ii);
    if diff<0
      from=mid+1;
    else
      to=mid-1;                       
    end
  end
  if nfft(ii)~=table(from)
      nfft(ii)=table(from-1);
  end;
  
  % Add back the missing factors of 2 (if any)
  nfft(ii)=nfft(ii)*2^n2reduce;
  
end;

tableout=table;


