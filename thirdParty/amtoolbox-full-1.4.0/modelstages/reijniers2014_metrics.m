function [varargout] = reijniers2014_metrics(doa, varargin)
%REIJNIERS2014_METRICS - extract localization metrics
%   Usage: [mean_error, bias] = reijniers2014_metrics(doa, parameters) 
%
%   Input parameters:
%     doa: Struct in returned from the reijiniers2014's model with
%          estimated and real directions of arrival
%
%   REIJNIERS2014_METRICS(...) returns psychoacoustic performance 
%   parameters for experimental response patterns. 
%   doa is a struct where actual and estimated directions of arrival must
%   be provided. If no input params are provided the returned metrics
%   resemble the ones provided in the original paper, see Reijiners et al. (2014).
%   This script is a wrapper for localizationerror.
%
%   If parameter is provided, reijniers2014_metric is a wrapper for localizationerror
%   with the parameter as the localization error. 
%
%   See also: demo_reijniers2014 reijniers2014 localizationerror
%
%   References:
%     R. Barumerli, P. Majdak, R. Baumgartner, J. Reijniers, M. Geronazzo,
%     and F. Avanzini. Predicting directional sound-localization of human
%     listeners in both horizontal and vertical dimensions. In Audio
%     Engineering Society Convention 148. Audio Engineering Society, 2020.
%     
%     R. Barumerli, P. Majdak, R. Baumgartner, M. Geronazzo, and F. Avanzini.
%     Evaluation of a human sound localization model based on bayesian
%     inference. In Forum Acusticum, 2020.
%     
%     J. Reijniers, D. Vanderleist, C. Jin, C. S., and H. Peremans. An
%     ideal-observer model of human sound localization. Biological
%     Cybernetics, 108:169--181, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/reijniers2014_metrics.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB SOFA M-Stats M-Image
%   #Author: Michael Sattler 
%   #Author: Roberto Barumerli (2020)
%   #Author: Clara Hollomey (2021)
% (adapted from code provided by Jonas Reijniers)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if isempty(varargin)
    nargout=2;
    num_src = size(doa.est, 1);
    num_exp = size(doa.est, 2);
    
    error = zeros(num_src, num_exp);

    for e = 1 : num_exp
        for s = 1 : num_src
            % calculate great circle error 
            % (x*y)=|x|*|y|cos(theta)
            x = doa.real(s,:);
            y = squeeze(doa.est(s,e,:));
            % NOTE: not accounting for distance! 
            % error(i,j) =  acos(x*y'/(norm(x)*norm(y))); 
            error(s, e) =  rad2deg(acos(x*y));
        end
    end

    % calculate ensemble values 
    mean_error = sum(error,2)/num_exp;
    
    % calculate summed estimated direction
    bias = squeeze(sum(doa.est,2));
    if size(bias,2) == 1
        bias = bias';
    end

    for s=1:num_src
        bias(s,:) = bias(s,:)/norm(bias(s,:)) - doa.real(s,:);
    end

    varargout{1}=mean_error;
    varargout{2}=bias;
elseif strcmp(varargin{1}, 'middle_metrics')
    % lateral_bias 
    exp.accL = reijniers2014_metrics(doa, 'accL'); 
    % lateral_rms_error
    exp.rmsL = reijniers2014_metrics(doa, 'rmsL'); 
    % elevation_bias
    exp.accE = reijniers2014_metrics(doa, 'accE'); 
    % local_rms_polar
    exp.rmsP = ... 
        reijniers2014_metrics(doa, 'rmsPmedianlocal'); 
    % quadrant_err
    exp.querr = ...
        reijniers2014_metrics(doa, 'querrMiddlebrooks'); 
    varargout{1} = exp;
else
    %% compute the metric relying on localizationerror.m
    definput.import={'localizationerror'};
    [flags,kv]=ltfatarghelper({'f','r'},definput,varargin);

    doa_real = SOFAconvertCoordinates(doa.real, 'cartesian', 'spherical');
    [lat_real,pol_real]=sph2horpolar(doa_real(:,1),doa_real(:,2));
    
    doa.est = squeeze(doa.est);
    
    if size(doa.est, 3) ~= 1 % remove third dimension relative to repetition of experiments
        num_exp = size(doa.est, 2);
        doa_est = [];
        
        for i = 1:size(doa.est, 2)
            est = squeeze(doa.est(:,i,:));
            doa_est = [doa_est; est];
        end
    else
        if(size(doa.est, 1) == size(doa.real, 1))
            num_exp = 1;
        else
            num_exp = size(doa.est, 1);
        end
        doa_est = doa.est;
    end
    
    doa_est = SOFAconvertCoordinates(doa_est, 'cartesian', 'spherical');
    [lat_est,pol_est] = sph2horpolar(doa_est(:,1),doa_est(:,2));
    
    m = zeros(size(doa_real, 1)*num_exp, 8);
    m(:, 1:2) = repmat(doa_real(:, [1 2]), num_exp, 1);
    m(:, 3:4) = doa_est(:, [1 2]);
    m(:, 5:6) = repmat([lat_real,pol_real], num_exp, 1);
    m(:, 7:8) = [lat_est,pol_est];

    % workaround since sph2horpolar can return complex numbers due to
    % numerical approximations
    m = real(m);
    
    if (strcmp(flags.errorflag, 'perMacpherson2003'))
        [varargout{1}, meta, par] = localizationerror(m, kv.f, kv.r, flags.errorflag);
    else
        [varargout{1}, meta, par] = localizationerror(m, flags.errorflag);
    end
    
    if (strcmp(flags.errorflag, 'sirpMacpherson2000'))
        varargout{2}=meta;
    end
    if length(varargout) > 1
        varargout{2}=meta;
    end
    if length(varargout) > 2
        varargout{3}=par; 
    end
end


