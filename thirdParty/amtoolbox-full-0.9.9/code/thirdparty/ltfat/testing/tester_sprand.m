function f=tester_sprand(varargin);
%TESTER_SPRAND   Sparse random numbers for testing.
%   Usage: f=tester_sprand(p1,p2);
%
%   Url: http://ltfat.github.io/doc/testing/tester_sprand.html

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

global LTFAT_TEST_TYPE

if isempty(LTFAT_TEST_TYPE)
    LTFAT_TEST_TYPE='double'
end;

f=sprand(varargin{:});

if strcmp(LTFAT_TEST_TYPE,'single')
    f=single(f);
end;
    


