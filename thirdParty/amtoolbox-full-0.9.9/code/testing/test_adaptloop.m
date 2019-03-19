function test_failed=test_adaptloop
%TEST_ADAPTLOOP  Test adaptloop
%
%  Test the adaptation loop implementation by comparing it to a reference
%  implementation.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_adaptloop.php

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

% The reference implementations are not dbspl-safe, so switch to the
% dboffset=100 setting.
  
save_dboffset = dbspl(1);

ltfatsetdefaults('dbspl','dboffset',100);

  
if 0
  f  = greasy;
  L=length(f);
  f=[randn(L,1), f];
  fs = 16000;
else
  fs=100;
  L=500;
  f=randn(L,2);
end;

limitr = [0.5, 10];

minlvl=0;

amt_disp(' ===============  TEST_ADAPTLOOP ================');

amt_disp('--- Used subroutines ---');

which comp_adaptloop

test_failed = 0;
for w=1:2
  for tauchoose=1:3
    switch(tauchoose)
     case 1
      tau=[0.005 0.050 0.129 0.253 0.500];
      tauname='DAU';
     case 2
      tau=linspace(0.005,0.5,5);
      tauname='BRE';
     case 3
      tau=[0.005 0.050 0.500];
      tauname='3RD';      
    end;
    
    for ii=1:length(limitr)
      
      limit=limitr(ii);
      
      slimit='     ';
      if limit>1
        slimit='LIMIT';
      end;
      
      

            
      s1=adaptloop(f(:,1:w),fs,limit,minlvl,tau);
            
      s2=ref_adaptloop_1(f(:,1:w),fs,limit,setdbspl(minlvl),tau);


      res=norm(s1(:)-s2(:))/norm(s1);
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf('REF_ADAPT %s %s L: %3i W:%2i %0.5g %s',slimit,tauname,L,w,res,fail);
      amt_disp(s);

      
      if tauchoose==1
        % Compute the reference solution, this is nlal_lim 
        % nlal_lim cannot handle matrices, so we do it manually.
        % nlal_lim can only handle the Dau loops.
        s3=zeros(L,w);
        for kk=1:w
          s3(:,kk)=nlal_lim(f(:,kk),fs,limit);    
        end;        
      
        res=norm(s1(:)-s3(:))/norm(s3);
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        s=sprintf('NLAL_LIM  %s %s L: %3i W:%2i %0.5g %s',slimit,tauname,L,w,res,fail);
        amt_disp(s);
      end;
                                          
      
      
      state = adaptloop_init(w,fs,limit,minlvl,tau);
      s4    = zeros(L,w);
      part  = floor(L/2)+10;
      
      [s4(1:part,:),  state] = adaptloop_run(state,f(1:part,1:w));
      [s4(part+1:L,:),state] = adaptloop_run(state,f(part+1:L,1:w));
      
      res=norm(s1(:)-s4(:))/norm(s1);
      
      %[9999,norm(s2(:)-s4(:))]
      
      [test_failed,fail]=ltfatdiditfail(res,test_failed);
      s=sprintf('BLOCK     %s %s L: %3i W:%2i %0.5g %s',slimit,tauname,L,w,res,fail);
      amt_disp(s);
      
    end;
  end;
end;

ltfatsetdefaults('dbspl','dboffset',save_dboffset);



