function test_failed=test_outermiddle
%TEST_DAU96  Test outermiddle
%
%  Test outermiddel by comparison with reference implementation.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_outermiddle.php

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

amt_disp(' ===============  TEST_OUTERMIDDLE ================');
  
test_failed=0;
  
fs=22050;

L=10*fs;

f=randn(L,1);

ref_HeadphoneFilter(fs);
ref_MiddleEarFilter(fs);

fr=OuterMiddleFilter(f);

outer_ear_fir_coeff=headphonefilter(fs);
mid_ear_fir_coeff=middleearfilter(fs,'jepsenmiddleear');

ff = fftfilt(outer_ear_fir_coeff,f);
ff = fftfilt(mid_ear_fir_coeff,ff);

res=norm(ff(:)-fr(:));

[test_failed,fail]=ltfatdiditfail(res,test_failed);

amt_disp(sprintf('OUTER+MIDDLE %0.5g %s',res,fail));


