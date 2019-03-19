% ci = AKci(data, ciLevel, dim)
% calculated the upper and lower bound of the confidence interval for
% values in data
%
% e.g.
% AK(data, .05)
% calculated th 95% confidence interval
%
% I N P U T:
% data    - vector or matrix of size [N M]. Confidence intervals are
%           caluclated for each row of data
% ciLevel - confidence interval level. 
%           The default of 95 computes the 95% condidence interval
% dim     - specifies the dimension to operate along
%           (default first non-singleton dimension)
%
% O U T P U T:
% ci   - convidence intervals. Matrix of size [N 2] where each row
%        specifies the ci for each row or column of data (depending on dim)
%
% 11/2016 - fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function ci = AKci(data, ciLevel, dim)

% check matlab dependencies
if ~license('test', 'Statistics_Toolbox') || ~exist('norminv.m', 'file')
    error('AKtools:MatlabToolboxes', 'AKci needs ''norminv.m'' from the Statistics and Machine Learning Toolbox')
end

% set default values
if ~exist('alpha', 'var')
    ciLevel = 95;
end
if ~exist('dim', 'var')
    dim = size(data);
    dim = find(dim>1, 1, 'first');
end

% mean and standard deviation
m = mean(data, dim);
s = std(data, 0, dim);

% correct the shape
m = reshape(m, [numel(m), 1]);
s = reshape(s, [numel(s), 1]);

% z-value
alpha  = 1-ciLevel/100;
z      = norminv([alpha/2 1-alpha/2],0,1);

% get the confidence intervall: m +/- z*s
ci = AKm(m, AKm(z, s, '*'), '+');

