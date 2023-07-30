function [varargout] = barumerli2023_metrics(varargin)
%BARUMERLI2023_METRICS extract localization metrics
%   Usage: metrics = barumerli2023_metrics(m, 'middle_metrics') 
%           metrics = barumerli2023_metrics(m, metric) 
%
%   Input parameters:
%     m       : (matrix) table organized as in localizationerror with actual
%               and predicted directions. The model barumerli2023 can
%               output directly such matrix.
%
%     metric  : (string) string indicating which metric has to be computed.
%
%
%   Output parameters (optional):
%      varagout : output as in localizationerror-m. If 'middle_metrics' is
%                 provided then a struct will be provided.
%
%   BARUMERLI2023_METRICS(...) returns psychoacoustic performance 
%   parameters for experimental response patterns. 
%   This script is a wrapper for localizationerror. This function is also
%   used internally in barumerli2023 which behavior is not here
%   described. 
%   For a complete list of supported metrics, please consider localizationerror. Moreover,
%   if 'middle_metrics' is provided the function returns a struct
%   containing the four metrics used in the paper Middlebrooks (1999).
%   There are: accuracy and root mean squared error for both
%   the lateral and polar dimensions and the quadrant error. 
%
%   See also: demo_barumerli2023 barumerli2023 localizationerror
%
%   References:
%     P. Majdak, M. J. Goupell, and B. Laback. 3-D localization of virtual
%     sound sources: Effects of visual environment, pointing method and
%     training. Atten Percept Psycho, 72:454--469, 2010.
%     
%     J. C. Middlebrooks. Virtual localization improved by scaling
%     nonindividualized external-ear transfer functions in frequency. J.
%     Acoust. Soc. Am., 106:1493--1510, 1999.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/barumerli2023_metrics.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: MATLAB SOFA M-STATISTICS M-Control M-Signal
%   #Author: Roberto Barumerli (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% parameters 
% m, 'middle_metrics'       return middlebrooks metrics as a struct
% doa, doa_real, 'm'        return m matrix
% doa, doa_real, <error>    compute error from `localizationerror`
% m, <error>                as above

if strcmp(varargin{2}, 'middle_metrics')
    assert(size(varargin{1}, 2) == 8, 'Please provide m matrix')
    m = varargin{1};
    % lateral_bias 
    exp.accL = localizationerror(m, 'accL'); 
    % lateral_rms_error
    exp.rmsL = localizationerror(m, 'rmsL'); 
    % elevation_bias
    exp.accP = localizationerror(m, 'accP'); 
    % local_rms_polar
    exp.rmsP = localizationerror(m, 'rmsPmedianlocal'); 
    % quadrant_err
    exp.querr = localizationerror(m, 'querrMiddlebrooks'); 
    varargout{1} = exp;
elseif strcmp(varargin{3}, 'm')
    assert(isfield(varargin{1}, 'estimations') & isa(varargin{2}, 'barumerli2023_coordinates'), 'If looking for m matrix please give doa as a struct and doa_real as barumerli2023_coordinates object(see barumerli2023)')
    varargout{1} = local_returnmatrixlocalizationerror(varargin{1}, varargin{2});
else
    if isfield(varargin{1}, 'estimations') && isa(varargin{2}, 'barumerli2023_coordinates')
        m = local_returnmatrixlocalizationerror(varargin{1}, varargin{2});
        errorflag = varargin{3};
    elseif size(varargin{1}, 2) == 8
        m = varargin{1};
        errorflag = varargin{2};
    else
        error('something went wrong!')
    end

    [varargout{1}, meta, par] = localizationerror(m, errorflag);
    
    if length(varargout) > 1
        varargout{2}=meta;
    end
    if length(varargout) > 2
        varargout{3}=par; 
    end

end 

function m = local_returnmatrixlocalizationerror(doa, doa_real)
    assert(size(doa.estimations, 3) == 3)
   
    doa_est_cart = barumerli2023_coordinates(reshape(doa.estimations, [], 3), 'cartesian');
    
    %% compute the metric relying on `localizationerror`
    doa_real_sph = doa_real.return_positions('spherical');
    doa_est_sph = doa_est_cart.return_positions('spherical');
    doa_real_hor = doa_real.return_positions('horizontal-polar');
    doa_est_hor = doa_est_cart.return_positions('horizontal-polar');
    
    num_rep = size(doa_est_cart.pos, 1)/size(doa_real.pos, 1);
    
    m = zeros(size(doa_real.pos, 1)*num_rep, 8);
    m(:, 1:2) = repmat(doa_real_sph(:, [1 2]), num_rep, 1);
    m(:, 3:4) = doa_est_sph(:, [1 2]);
    m(:, 5:6) = repmat(doa_real_hor(:,[1 2]), num_rep, 1);
    m(:, 7:8) = doa_est_hor(:, [1 2]);


