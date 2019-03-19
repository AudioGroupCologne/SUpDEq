function definput=arg_baumgartner2014_pmv2ppp(definput)

definput.keyvals.p=ones(72,44);
definput.keyvals.rang=-90:5:269;
definput.keyvals.tang=[-30:5:70,80,100,110:5:210];
definput.keyvals.exptang=[];

definput.flags.print = {'noprint','print'};
definput.flags.chance = {'','chance'};
definput.flags.ppp = {'','QE_PE_EB','QE','PE','EB','absPE'};

%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/defaults/arg_baumgartner2014_pmv2ppp.php

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

