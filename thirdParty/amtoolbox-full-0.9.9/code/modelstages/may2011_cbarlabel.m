function cbarlabel(szString,fig)
%name   Short description...
%
%USAGE 
%   OUT = NAME(IN1,IN2)
%
%INPUT ARGUMENTS
%   IN1 : input parameter
%   IN2 : second parameter, (default, IN2 = 5)
%
%OUTPUT ARGUMENTS
%   OUT : results [nFrames x 1]
%
%NOTE
%   This is a template.
% 
%REFERENCES
%   [1] P. Boersma, "Accurate Short-Term Analysis of the Fundamental
%       Frequency and the Harmonics-to-Noise Ratio of Sampled Sounds", IFA
%       Proceedings 17, 1993.
%
%   See also ...
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/may2011_cbarlabel.php

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

%   Developed with Matlab 7.8.0.347 (R2009a). Please send bug reports to:
%   
%   Author  :  Tobias May, 2009 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :
%   v.1.0   2009/08/6
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
errMessage = nargchk(1,2,nargin);

% Display error message
if ~isempty(errMessage); 
    % Display help file with html link ... or conventional help message
    if   exist('showHelp.m','file'); showHelp(mfilename); 
    else help(mfilename); end
    % Display error message
    error(errMessage); 
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
        error('No colorbar was detected ...')
    case 1
        set(get(allH(isCBar),'ylabel'),'String',szString);
    otherwise
        error('Confusion due to multiple colorbars ...')
end

h = allH(isCBar);

