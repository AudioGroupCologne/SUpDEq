function data = data_brimijoin2013
%DATA_BRIMIJOIN2013 Data from Brimijoin et al. (PLoS One, 2013)
%   Usage: data = data_brimijoin2013 
%
%   Output parameters:
%     data    : structure containing the data
%
%
%   data_brimjoin2013 provides the mean externalization likelihood of NH listeners extracted from Fig. 6
%
%   The data struct comprises:
%
%     'response'    mean externalization likelihood
%     'legend'      conditions (columns of response matrix)
%     'mix'         head-present re head-absent mix
%     'xlabel'      abscissa label of original figure
%     'ylabel'      ordinate label of original figure
%
%   References:
%     W. O. Brimijoin, A. W. Boyd, and M. Akeroyd. The contribution of head
%     movement to the externalization and internalization of sounds. PLOS
%     ONE, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_brimijoin2013.php


%   #Author: Robert Baumgartner

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


data.mix = 0:.2:1;
data.xlabel = 'Head-present / Head-absent Mix';
data.ylabel = 'Proportion externalized responses';
data.legend = {'Head: fixed; Signal: fixed';'Head: fixed; Signal: move';...
  'Head: move; Signal: fixed';'Head: move; Signal: move'};
data.response = ...
  [.08,.09,.18,.55,.70,.78;...
   .07,.08,.29,.51,.79,.73;...
   .13,.07,.08,.17,.32,.38;...
   .27,.51,.54,.81,.89,.89]';


end


