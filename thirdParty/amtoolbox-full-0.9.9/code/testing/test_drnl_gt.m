% This script compares the implentation of gammatone filter coefficiens
% supplied for the DRNL with the standard gammatone implementation.
%
% The DRNL filters seems to have must higher order, and they are not
% all-pole.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_drnl_gt.php

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


fc=1000;
fs=16000;
n=4;
bw=100;

[b_ref,a_ref]=coefGtDRNL(1000,100,4,16000);

[b,a]=gammatone(fc,fs,n,bw,'complex');
b=real(b);
a=real(a);




