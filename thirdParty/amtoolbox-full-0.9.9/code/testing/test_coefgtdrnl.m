function test_failed=test_coefgtdrnl
%TEST_COEFGTDRNL
%
%  Test how to replace the coefGtDRNL function by gammatone.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_coefgtdrnl.php

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

amt_disp(' ===============  TEST_COEFGTDRNL ================');

test_failed=0;

fc=2000;
fs=10000;
n=4;
bw=200;

[b1,a1]=coefGtDRNL(fc,bw,n,fs);
[b2,a2]=gammatone(fc,fs,n,bw/audfiltbw(fc),'classic');

res=norm(b1-b2);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GTDRNL B %f %s\n', res,fail);

res=norm(a1-a2);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GTDRNL A %f %s\n', res,fail);



