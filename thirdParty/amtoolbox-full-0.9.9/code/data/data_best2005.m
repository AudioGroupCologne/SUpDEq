function data = data_best2005
%DATA_BEST2005  Listener averages of absolute polar angle error and SCC
%   Usage: data = data_best2005
%
%   Output parameters:
%     data.qe     : quadrant error rates
%     data.ape    : absolute polar angle error
%     data.seape  : standard error of absolute polar angle error
%     data.meta   : condition labels corresponding to data entries
%
%   DATA_BEST2005 returns listener averages of absolute polar angle error
%   and SCCs from Best et al. (2005): Fig.5(b), Fig.10(b), and Tab.2
%
%
%   References best2005highfrequency
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_best2005.php

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

% AUTHOR: Robert Baumgartner

% Mean quadrant error rates from Tab.2 (Exp. I) and Tab. 4 (Exp. II)
data.qe =     [3.0,mean([8.4,6.7]),12.0,13.8,16.4];

% Mean absolute polar angle errors from Fig.5b (Exp. I) and Fig.10b (Exp. II)
data.ape =    [21,30,38,42.5,46];
data.seape =  ones(1,5);

data.meta =   {'BB noise','BB speech','-20dB','-40dB','-60dB'};
end
