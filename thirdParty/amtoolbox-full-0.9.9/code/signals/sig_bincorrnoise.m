function outsig = sig_bincorrnoise(siglen,coher,varargin)
% sig_bincorrnoise  Binaurally correlated noise
%   Usage: outsig = sig_bincorrnoise(siglen,coher)
%
%   Input parameters:
%       siglen    : Number of samples of outsig
%       coher     : Interaural coherence of the produced signal.
%
%   Output parameters:
%       outsig    : nsig x2 correlated noise signal
%
%   SIG_BINCORRNOISE(siglen,coher) will generate a interaurally correlated noise signal 
%   with coherence coher. The output is a 2 column matrix of length siglen.
%
%   SIG_BINCORRNOISE(siglen,coher,...) will pass all additional parameters
%   onto the noise function to select between different types of stochastic
%   noise.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_bincorrnoise.php

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

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameters ------------------------------------

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ( ~isnumeric(siglen) || ~isscalar(siglen) || siglen<0 )
    error('%s: siglen has to be a positive scalar.',upper(mfilename));
end

if ( ~isnumeric(coher) || ~isscalar(coher) || coher<0)
    error('%s: coher has to be a positive scalar.',upper(mfilename));
end


% ------ Computation -----------------------------------------------------

% Generate correlation matrix
R = [1 coher; coher 1];
% Eigen decomposition
[V,D] = eig(R);

% Form correlating filter
W = V*sqrt(D);

% Generate uncorrelated noise
n = noise(siglen,2,varargin{:});

% Correlate the noise
outsig = n * W';

