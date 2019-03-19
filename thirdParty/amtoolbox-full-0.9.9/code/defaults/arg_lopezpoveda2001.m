function definput=arg_lopezpoveda2001(definput)
 
  
  definput.keyvals.flow=80;
  definput.keyvals.fhigh=8000;
  definput.keyvals.basef=[];
  definput.keyvals.bwmul=1;

  % parameters according to Lopez-Poveda and Meddis 2001
  definput.keyvals.lin_ngt = 2; 
  definput.keyvals.lin_nlp = 4; 
  definput.keyvals.lin_fc = [-0.06762 1.01679];
  definput.keyvals.lin_bw = [  .03728 .78563];
  definput.keyvals.lin_gain = [4.20405 -.47909];
  definput.keyvals.lin_lp_cutoff = [-0.06762 1.01679 ];
  
  definput.keyvals.nlin_ngt_before = 3;
  definput.keyvals.nlin_ngt_after = [];
  definput.keyvals.nlin_nlp = 3;
  definput.keyvals.nlin_fc_before = [-0.05252 1.01650];
  definput.keyvals.nlin_fc_after  = [];
  definput.keyvals.nlin_bw_before = [-0.03193 .77426 ];
  definput.keyvals.nlin_bw_after  = [];
  definput.keyvals.nlin_lp_cutoff = [-0.05252 1.01650];
  
  definput.keyvals.nlin_a = [1.40298 .81916 ];
  definput.keyvals.nlin_b = [1.61912 -.81867 ];
  definput.keyvals.nlin_c = [log10(.25) 0];
  definput.keyvals.nlin_d = 1;

  definput.keyvals.compresslimit = [];
  
  definput.flags.middleear={'middleear','nomiddleear','jepsenmiddleear'};
  
  definput.flags.path = {'bothparts','linonly','nlinonly'};

  % This parameter set is not supported anymore, as there is no evince as
  % to whether or not this is actually the dataset that whas used in the paper.  
  definput.groups.jepsen2008={...
      'lin_bw',          [ .03728   .75 ],...
      'lin_lp_cutoff',   [-0.06762 1.01 ],...  
      'nlin_bw_before' , [-0.03193  .77 ],...
      'compresslimit', 1500, ...
      'jepsenmiddleear', ...
      'nlin_ngt_before', 2, ...
      };

%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/defaults/arg_lopezpoveda2001.php

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

