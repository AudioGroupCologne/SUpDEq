% [N_crit, alpha, beta, power] = AKnAFCstatistics(p_guess, p_pop, N_trials, alpha)
% this script calculates the number of correct answers that is needed in an
% N-AFC-test to assume the H1, that differences between two stimuli are
% audible [1].
% It further gives the actual alpha- and beta-error-levels as well as the
% power (1-beta). Calling this function without output arguments will print
% results to the command window.
%
% I N P U T:
% p_guess  - guessing probability between 0 and 1. This is 1/N
% p_pop    - assumed detection rate in the population between 0 and 1. This
%            affects the beta error (type II), and the power
% N_trials - number of trials
% alpha    - desired alpha error level (type I error level)
%
% O U T P U T:
% N_crit   - number of correct answers that is at least needed to assume
%            audible differences
% alpha    - actual alpha level (in binomial test statistics the alpha
%            level is quantized)
% beta     - beta error level
% power    - test power (1-beta)
%
%
% Caluclations are implemented and validet using
% [1] Les Leventhal: "Type I and type 2 errors in the statistical analysis
%     of listening tests." J. Audio Eng. Soc. 34(6):437-453 (1986).
%
% 07/2013 - fabian.brinkmann@tu-berlin.de

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
function [N_crit, alpha, beta, power] = AKnAFCstatistics(p_guess, p_pop, N_trials, alpha)
% ---------------------------------------------------------------- fill out
if nargin == 0
    p_guess  = .5;  % guessing propability
    p_pop    = .66; % assumed propbility for detection in population
    N_trials = 20;  % number of trials
    alpha    = .05; % alpha-error level
end
% -------------------------------------------------------------------------

% calculation using binocdf and binopdf (robust with large N_trials)
p_N_guessed_cumsum = binocdf(N_trials:-1:0, N_trials, 1-p_guess)';

N_crit = find(p_N_guessed_cumsum < alpha, 1, 'first') - 1;
alpha  = p_N_guessed_cumsum(N_crit+1);

beta  = binocdf(N_crit-1, N_trials, p_pop);
power = 1-beta;

if nargout == 0
    disp(['Assume H1, if ' num2str(N_crit) '/' num2str(N_trials) ' or more answers are correct (' num2str(round(N_crit/N_trials*100000)/1000) ' percent).'])
    disp(['(alpha=' num2str(round(alpha*10000)/10000)    ...
         ', beta=' num2str(round(beta*10000)/10000)      ...
         ', power=' num2str(round((1-beta)*10000)/10000) ...
         ', p_pop=' num2str(p_pop) ')'])
    clear N_crit alpha beta power
end