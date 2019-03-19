function test_failed=test_dbsplsafety
%TEST_DBSPLSAFETY  Test dbsplsafety
%
%  Test if the non-linear stages are safe with respect to changes in the
%  reference level.
%
%  Any function using DBSPL or SETDBSPL should be tested by this tester.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_dbsplsafety.php

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

amt_disp(' ===============  TEST_DBSPLSAFETY ================');

test_failed=0;

ref_level=77;

refl=[-20*log10(20e-6),100,0];

save_dboffset = dbspl(1);

ltfatsetdefaults('dbspl','dboffset',ref_level);

% ----------------- test basic functionality 

for ii=1:length(refl)
  
  dboffset=refl(ii);

  ltfatsetdefaults('dbspl','dboffset',dboffset);

  x=randn(200,1);
  x_orig=dbspl(x);
  y=setdbspl(x,70);
  y=setdbspl(y,x_orig);
  
  res=norm(y-x);
  [test_failed,fail]=ltfatdiditfail(res,test_failed);  
  s=sprintf('BASIC     %0.5g %0.5g %s',dboffset,res,fail);
  amt_disp(s);   
  
end;

% ----------------- test adaptloop ----------------------

ltfatsetdefaults('dbspl','dboffset',ref_level);

f=randn(8000,1);

ref_adapt = adaptloop(setdbspl(f,80),8000);

for ii=1:length(refl)
  
  dboffset=refl(ii);
  
  ltfatsetdefaults('dbspl','dboffset',dboffset);

  f_adapt = adaptloop(setdbspl(f,80),8000);
  
  res=f_adapt-ref_adapt;
  res=norm(res(:));
  [test_failed,fail]=ltfatdiditfail(res,test_failed);
  s=sprintf('ADAPTLOOP %0.5g %0.5g %s',dboffset,res,fail);
  amt_disp(s);

end;

% ----------------- test DRNL ----------------------

ltfatsetdefaults('dbspl','dboffset',ref_level);

f=randn(10000,1);
fs=32000;

ref_drnl = drnl(setdbspl(f,80),fs);

for ii=1:length(refl)
  
  dboffset=refl(ii);
  
  ltfatsetdefaults('dbspl','dboffset',dboffset);

  f_drnl = drnl(setdbspl(f,80),fs);
  
  res=f_drnl-ref_drnl;
  res=norm(res(:));
  [test_failed,fail]=ltfatdiditfail(res,test_failed,1e-6);
  s=sprintf('DRNL      %0.5g %0.5g %s',dboffset,res,fail);
  amt_disp(s);

end;

ltfatsetdefaults('dbspl','dboffset',save_dboffset);



