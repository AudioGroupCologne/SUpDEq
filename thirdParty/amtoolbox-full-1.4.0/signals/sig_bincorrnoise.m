function outsig = sig_bincorrnoise(siglen,coher,varargin)
%SIG_BINCORRNOISE  Binaurally correlated noise
%   Usage: outsig = sig_bincorrnoise(siglen,coher)
%
%   Input parameters:
%       siglen    : Number of samples of outsig
%       coher     : Interaural coherence of the produced signal.
%
%   Output parameters:
%       outsig    : nsig x 2 correlated noise signal
%
%   SIG_BINCORRNOISE(siglen,coher) will generate a interaurally correlated noise signal 
%   with coherence coher. The output is a 2 column matrix of length siglen.
%
%   SIG_BINCORRNOISE(siglen,coher,...) will pass all additional parameters
%   onto the noise function to select between different types of stochastic
%   noise.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_bincorrnoise.php


%  #Author: Hagen Wierstorf (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


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


