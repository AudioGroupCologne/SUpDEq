function [flags,kv]=amt_flags(varargin)
%amt_flags Returns the start-up flags of the AMT
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/amt_flags.php

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
  
%   AUTHOR : Piotr Majdak 

persistent AMT;
  
if ~isempty(varargin),
    % create persistent variable with flags
  if iscell(varargin{1}), varargin=varargin{1}; end
  definput.import={'amt_cache','amt_disp'};
  [flags,kv]=ltfatarghelper({},definput,varargin);
  AMT.flags=flags;
  AMT.kv=kv;
elseif isempty(AMT)
  definput.import={'amt_cache','amt_disp'};
  [flags,kv]=ltfatarghelper({},definput,[]);
  AMT.flags=flags;
  AMT.kv=kv;  
end
flags=AMT.flags;
kv=AMT.kv;

