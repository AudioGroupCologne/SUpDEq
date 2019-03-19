function f=iufilterbank(varargin);  
%IUFILTERBANK  Filter bank inversion, DEPRECATED
%   Usage:  f=iufilterbank(c,g,a);
%
%   IUFILTERBANK has been deprecated by IFILTERBANK. Call IFILTERBANK
%   with the exact same parameters as the old call to IUFILTERBANK.
%
%   See also: ifilterbank
%
%   Url: http://ltfat.github.io/doc/deprecated/iufilterbank.html

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

warning(['LTFAT: IUFILTERBANK has been deprecated, please use IFILTERBANK ' ...
         'instead.']);
  

f=ifilterbank(varargin{:});


