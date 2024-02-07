function varargout = comp_may2011_gammatone(varargin)
%comp_gammatone   Private method, which is normally shadowed by the
%   corresponding MEX routine. It will ONLY be executed, if the underlying
%   MEX routine is not available.  
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/mex/comp_may2011_gammatone.php


%   #Author: Tobias May (2009)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% Get M-file name
mfile = mfilename;

% Show help
try
    % Display help file with html link
    showHelp(mfile(1:end-3));
catch exeption %#ok
    % Conventional help message
    help(mfile(1:end-3));
end
    
% Report error message
error(['MEX function "',mfilename,'" is not available for your ',...
       'operating system (',computer,'). Run "amt_mex" to built ',...
       'required MEX binaries.'])

%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************

