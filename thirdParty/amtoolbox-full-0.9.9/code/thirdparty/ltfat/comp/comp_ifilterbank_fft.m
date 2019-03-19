function F = comp_ifilterbank_fft(c,G,a)
%COMP_IFILTERBANK_FFT  Compute filtering in FD
%
%
%
%   Url: http://ltfat.github.io/doc/comp/comp_ifilterbank_fft.html

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

W = size(c{1},2);
M = numel(G);
L = numel(G{1});
F = zeros(L,W,assert_classname(c{1},G{1}));

for m=1:M
   for w=1:W
      % This repmat cannot be replaced by bsxfun
      F(:,w)=F(:,w)+repmat(fft(c{m}(:,w)),a(m),1).*conj(G{m});
   end;
end

