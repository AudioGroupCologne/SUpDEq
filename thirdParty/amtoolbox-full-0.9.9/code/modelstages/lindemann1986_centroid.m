function d = lindemann1986_centroid(cc)
%lindemann1986_centroid Calculates the centroid for a cross-correlation
%   Usage: d = lindemann1986_centroid(cc)
%
%   Input parameters:
%       cc  : Lindemann cross-correlation. Dim: 1 x delay line length
%
%   Output parameters:
%       d   : lindemann1986_centroid in the range -1..1~ms
%
%   LINDEMANN1986_CENTROID(cc) calculates the centroid for a given
%   cross-correlation from the Lindemann model.
%
%   The centroid is computed by (see Lindemann (1986a), page 1613, eq. 22):
%
%              M                  M
%       d = ( sum m*Psi(m) ) / ( sum Psi(m) )
%             m=-M               m=-M
%
%   where M is half the length of the delay line -M,...,M.
%
%   See also: lindemann1986
%
%   References:
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608-1622, 1986.
%     
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/lindemann1986_centroid.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
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

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input  parameters -----------------------------------

error(nargchk(1,1,nargin));

if ~isnumeric(cc) || ~isvector(cc)
    error('%s: cc has to be a numeric vector signal!',upper(mfilename));
end
% Ensure size(cc) = delay line length x 1
if size(cc,1)==1
    cc = cc';
end


% ------ Computation -----------------------------------------------------
% Calculate the length of the delay line as -M:M
m = linspace(-1,1,length(cc))';
% Calculate the centroid using the -M:M delay line
d = sum(m.*cc) / sum(cc);


