function [f,stat]=urlwrite(webfn,fn)

%   #Author: The AMT Team
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/legacy/urlwrite.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
if exist('OCTAVE_VERSION','builtin') ~= 0,
    f=dbstack('-completenames');
    fx=f(2).file;
    if ~strcmp('urlwrite.m', fx(end-9:end))
      [f,stat]=builtin('urlwrite' ,webfn,fn);
    end
else
  f=websave(fn,webfn);
  stat=1;
end

