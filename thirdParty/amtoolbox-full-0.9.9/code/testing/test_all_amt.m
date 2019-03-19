
tests_todo={ 'adaptloop', 'outermiddle', 'coefgtdrnl','dbsplsafety'};


total_tests_failed=0;
list_of_failed_tests={};

for ii=1:length(tests_todo)
  test_failed=feval(['test_',tests_todo{ii}]);
  total_tests_failed=total_tests_failed+test_failed;
  if test_failed>0
    list_of_failed_tests{end+1}=['test_',tests_todo{ii}];
  end;
end;

amt_disp(' ');
if total_tests_failed==0
  amt_disp('ALL TESTS PASSED');
else
  s=sprintf('%i TESTS FAILED',total_tests_failed);
  amt_disp(s);
  amt_disp('The following test scripts contained failed tests');
  for ii=1:length(list_of_failed_tests)
    amt_disp(['   ',list_of_failed_tests{ii}]);
  end;
end;



%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_all_amt.php

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

