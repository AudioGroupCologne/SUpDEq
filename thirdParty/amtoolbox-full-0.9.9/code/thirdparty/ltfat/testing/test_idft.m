function test_failed=test_idft
Lr=[1, 19, 20];


test_failed=0;

disp(' ===============  TEST_IDFT ==============');

for jj=1:length(Lr)
  L=Lr(jj);
    for n = 1:2
    
    if (n==1)
       type = 'complex';
       c=tester_crand(L,1);
    elseif (n==2)
       type = 'real';
       c=tester_rand(L,1);      
    end
    
    f1=idft(c);
    f2=ref_idft(c);
    
    res=norm(f1-f2);
    [test_failed,fail]=ltfatdiditfail(res,test_failed);        
    s=sprintf('IDFT %6s  L:%3i %0.5g %s',type,L,res,fail);
    disp(s);
    end
  end;
end;



%
%   Url: http://ltfat.github.io/doc/testing/test_idft.html

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

