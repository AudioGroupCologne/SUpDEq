function amt_disp(msg,flag)
%amt_disp AMT-specific overload of the function 'disp'
%   Usage: amt_disp(X);
%     amt_disp(X,'progress');
%     amt_disp(X,'volatile');
%
%   AMT_DISP(X); can be used to show message X in the command
%   window. The output of AMT_DISP depends on the start-up
%   configuration of the AMT. 
%     
%   When the AMT is started in the verbose mode (default mode), AMT_DISP 
%   will always display. When the AMT is started in the documentation mode, 
%   AMT_DISP will display only if no flag is provided. When the AMT is 
%   started in the silent mode, AMT_DISP will never display. See
%   amt_start for further explanation on the start-up configurations. 
%
%   AMT_DISP(X,'progress'); can be used as progress indicator. It will be
%   shown during the normal operation but supressed when used to create the
%   documentation.
%
%   AMT_DISP(X,'volatile'); can be used as volatile progress indicator.
%   Any subsequent call of the AMT_DISP will delete the previous volatile
%   message. This way a changing progress can be clearly shown even in loops. 
% 
%   See also: amt_start
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/amt_disp.php

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

  
%   Author: Piotr Majdak, 2016
%   last change: 17.5.2017, volatile added

persistent CachedMsg;

if exist('flag','var')
  if ~strcmp(flag,'progress') && ~strcmp(flag,'volatile'),
    error(['Unsupported flag ' flag]);
  end
else
  flag='';
end

flags=amt_flags;

if flags.do_verbose
  if strcmp(flag,'volatile')
    if ~isempty(CachedMsg),
      reversemsg = repmat(sprintf('\b'), 1, length(CachedMsg));
      fprintf(reversemsg);
    end
    fprintf(strrep(msg,'\','\\'));
    CachedMsg=msg;
  else
    if ~isempty(CachedMsg),
      reversemsg = repmat(sprintf('\b'), 1, length(CachedMsg));
      fprintf(reversemsg);
      CachedMsg=[];
    end
    disp(msg);
  end
end

if flags.do_documentation
  if ~strcmp(flag,'progress'), 
    disp(msg); 
  end
end

if flags.do_silent
  % do nothing
end
