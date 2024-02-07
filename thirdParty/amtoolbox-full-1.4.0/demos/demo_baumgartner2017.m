%DEMO_BAUMGARTNER2017
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_baumgartner2017.php


%   #Author: Clara Hollomey (2020): for testing purposes

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

flags = {'boyd2012', 'hartmann1996', 'hassager2016', 'baumgartner2017'};

flag = flags{1};

data = exp_baumgartner2017(flag);

