%% test installation mode
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/test_amt_start.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: All messages displayed, install toolboxes');
amt_start('install');
disp('TEST AMT_START: DEMO_BAUMGARTNER2013 should run and display results now'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display a line'); 
amt_stop;

%% test default verbose mode
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: All messages displayed');
amt_start;
disp('TEST AMT_START: DEMO_BAUMGARTNER2013 should display results now'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display a line'); 
amt_stop;

%% test silent mode
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: SILENT MODE: no messages displayed');
amt_start('silent');
disp('DTEST AMT_START: DEMO_BAUMGARTNER2013 should display nothing'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display nothing'); 
amt_stop;

%% test documentation mode
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: DOCUMENTATION MODE should display as normal');
amt_start('documentation');
disp('TEST AMT_START: DEMO_BAUMGARTNER2013 should display nothing'); 
demo_baumgartner2013;
disp('TEST AMT_START: AMT_STOP should display nothing'); 
amt_stop;

%% test cache
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: NOW TESTING CACHE SETTINGS');
disp(' ');disp(' ');disp(' ');
amt_start('silent', 'redo');
[flags, keyvalues] = amt_configuration;
assert(strcmp(flags.cachemode, 'redo'))
assert(strcmp(flags.disp, 'silent'))
amt_stop

amt_start('redo');
[flags, keyvalues] = amt_configuration;
assert(strcmp(flags.cachemode, 'redo'))
assert(strcmp(flags.disp, 'verbose'))
amt_stop

amt_start('localonly');
[flags, keyvalues] = amt_configuration;
assert(strcmp(flags.cachemode, 'localonly'))
assert(strcmp(flags.disp, 'verbose'))
amt_stop

amt_start;
[flags, keyvalues] = amt_configuration;
assert(strcmp(flags.cachemode, 'normal'))
assert(strcmp(flags.disp, 'verbose'))
amt_stop

amt_start('cached');
[flags, keyvalues] = amt_configuration;
assert(strcmp(flags.cachemode, 'cached'))
assert(strcmp(flags.disp, 'verbose'))
amt_stop

disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: CACHE SETTINGS FINISHED');
disp(' ');disp(' ');disp(' ');

%% test done
disp(' ');disp(' ');disp(' ');
disp('TEST AMT_START: FINISHED');


