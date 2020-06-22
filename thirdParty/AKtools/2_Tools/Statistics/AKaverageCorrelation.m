% [r, rAll, R] = AKaverageCorrelation(data, do_plot)
%
% calculates and averages intersubject-correlations for 1-factorial
% repeated measures ANOVA ratings. This is needed for caluclation of
% optimal sample size. Correlations are transformed to fisher's z values
% bevore averging and inverse transformed afterwards.
% If you have multiple factorial ANOVA, calculate  correlations seperately
% for each factor.
%
% INPUT
% data      - [subjects x ratings], this means a column contains ratings of
%             all subjects for one factor level; a row contains ratings
%             of one subject for all factor levels
% do_plot   - plot correlations (true, false=default)
%
% OUTPUT
% r         - mean intersubject correlation across factor levels
% rAll      - development of mean intersubject correlation across subjects
% R         - matrix of correlation between subjects
%
% 01/2014 - fabian.brinkmann@tu-berlin.de

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
function [r, rAll, R] = AKaverageCorrelation(data, do_plot)

if ~exist('do_plot', 'var')
    do_plot = false;
end

r        = zeros(size(data, 1)-1,1);
highCorr = false;

for n = 2:size(data, 1)
    % number of correlations to average
    N = sum(sum( triu(ones(n), 1) ));
    % get correlations
    tmp_r = corrcoef(data(1:n,:)');
    if n == size(data, 1)
        R = tmp_r;
    end
    % check for 1 or -1
    id = abs( triu(tmp_r, 1) ) == 1;
    if any(id)
        tmp_r(id) = tmp_r(id) * .999;
        highCorr = true;
    end
    % get upper triangular part
    tmp_r = triu(tmp_r, 1);
    % fisher's z-transformation
    tmp_r = 1/2 * log((1+tmp_r) ./ (1-tmp_r));
    % arithmetic mean
    tmp_r = sum(sum(tmp_r)) / N;
    % inverse fisher's z-transformation
    tmp_r = (exp(2*tmp_r)-1) ./ (exp(2*tmp_r)+1);
    % output
    r(n-1) = tmp_r;   
end

if highCorr
    warning('AKaverageCorrelation:Corr', 'Correlations of 1 or -1 appeared.')
end

if do_plot
    AKf(20,10);
    plot(2:size(data,1), r, '.-k')
    xlabel('subjects')
    set(gca, 'xTick', 2:size(data,1))
    ylabel('correlation')
    axis tight
    ylim([0 1])
end

% set output
rAll = r;
r    = r(end);

% clear output if not desired
if nargout == 0
    clear
end