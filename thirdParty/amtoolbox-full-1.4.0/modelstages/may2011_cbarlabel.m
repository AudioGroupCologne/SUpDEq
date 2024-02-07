function may2011_cbarlabel(szString,fig)
%MAY2011_CBARLABEL sets the labels of a cbar plot
%
%   Usage: 
%     may2011_cbarlabel(szString,fig)
%
%   Input parameters:
%     szString : ordinate label
%     fig      : figure number
%
%   Output parameters:
%      OUT : results [nFrames x 1]
%
%
%   MAY2011_CBARLABEL is a helper function to modify cbar plots.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/may2011_cbarlabel.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Tobias May (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
if nargin < 1 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Set default values
if nargin < 2 || isempty(fig); fig = gcf; end

% Look for color handel
allH = get(fig,'children');

% Number of handels
nHandles = length(allH);

% Allocate logical operator
isCBar = false(nHandles,1);

% Loop over objects
for ii = 1 : length(allH)
    if strmatch('Colorbar',get(allH(ii),'tag'));
        isCBar(ii) = true;
    end       
end

% Set colormap label
switch isequal(1,sum(isCBar))
    case 0
        warning('No colorbar was detected ...')
    case 1
        set(get(allH(isCBar),'ylabel'),'String',szString);
    otherwise
        error('Confusion due to multiple colorbars ...')
end

h = allH(isCBar);


