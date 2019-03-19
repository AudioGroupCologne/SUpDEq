function c = vect2cell(x,idx)
%VECT2CELL Vector to cell
%
%   Works exactly like mat2cell(x,idx,size(x,2))
%   but it is faster.
%
%   Url: http://ltfat.github.io/doc/comp/vect2cell.html

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

if sum(idx)~=size(x,1)
    error('%s: Sizes do not comply.',upper(mfilename));
end

idxEnd = cumsum(idx(:));
idxStart = [1;1+idxEnd(1:end-1)];
c = arrayfun(@(idS,idE) x(idS:idE,:),idxStart,idxEnd,'UniformOutput',0);

