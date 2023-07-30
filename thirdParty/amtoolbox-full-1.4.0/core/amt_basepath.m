function bp = amt_basepath;
%AMT_BASEPATH  The base path of the AMT installation
%   Usage: bp = amt_basepath;
%
%   AMT_BASEPATH returns the top level directory in which the AMT
%   files are installed.
%
%   See also: amt_start amt_auxdatapath
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_basepath.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
%   #Author: Piotr Majdak (2015)

basepath = which('amt_start.m');
bp = basepath(1:end-numel('amt_start.m'));


