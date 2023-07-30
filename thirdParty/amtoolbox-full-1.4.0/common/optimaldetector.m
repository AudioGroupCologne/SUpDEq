function mue = optimaldetector(ir_stim,template)
%OPTIMALDETECTOR  Generic optimal detector for the CASP and Breebaart models
%
%   This is a correlation-based optimal detector for a signal known exactly.
%   See Green & Swets (1966) for more details.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/optimaldetector.php


%   #Author: Peter L. Soendergaard (2011)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

corrmue = ir_stim.*template;
optfactor = sqrt(numel(corrmue));

% Take mean over all dimensions of internal representation and correct for
% optimalityfactor.
mue = mean(corrmue(:))*optfactor;


%OLDFORMAT


