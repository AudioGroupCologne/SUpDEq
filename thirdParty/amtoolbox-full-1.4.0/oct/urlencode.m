function utf=urlencode(str)
% This function emulates the MATLAB function urlencode by encoding 
% the string STR to UTF-8. All special characters are encoded to 
% %HH (ASCII in hex), except: . _  and -. Further,   is encoded to +.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/oct/urlencode.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% #Author: Piotr Majdak (2017)
  utf = '';
  
if isoctave
  for ii = 1:length(str),
    if isalnum(str(ii)) || str(ii)=='.' || str(ii)=='_' || str(ii)=='*' || str(ii)=='-' 
      utf(end+1) = str(ii);
    elseif str(ii)==' '
      utf(end+1) = '+';
    else
      utf=[utf,'%',dec2hex(str(ii)+0)];
    end 	
  end  
else
  for ii = 1:length(str),
    if isletter(str(ii)) || str(ii)=='.' || str(ii)=='_' || str(ii)=='*' || str(ii)=='-' 
      utf(end+1) = str(ii);
    elseif str(ii)==' '
      utf(end+1) = '+';
    elseif str2num(str(ii)) < 10
      utf(end+1) = str(ii);  
    else
      utf=[utf,'%',dec2hex(str(ii)+0)];
    end 	
  end  
end


