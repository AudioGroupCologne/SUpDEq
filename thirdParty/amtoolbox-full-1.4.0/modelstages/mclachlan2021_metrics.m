function [metric, m] = mclachlan2021_metrics(doa,varargin)
%MCLACHLAN2021_METRICS - extract localization metrics
%   Usage: metric = mclachlan2021_metrics(doa) 
%
%   Input parameters:
%     doa               : Struct returned from the mclachlan2021 model with
%                         estimated and real directions of arrival
%
%   Output parameters:
%     metric            : metrics for evaluation of localisation
%                         performance
%
%   MCLACHLAN2021_METRICS(...) returns psychoacoustic performance 
%   parameters for experimental response patterns. 
%   doa is a struct where actual and estimated directions of arrival must
%   be provided.
%
%   The metrics struct contains the following fields:
%   
%     .reversal_fb      Percentage of front-back reversals, omitting
%                       small errors occurring within 10 degrees of the
%                       plane dividing the hemifields
%
%     .reversal_ud      Percentage of up-down reversals, omitting
%                       small errors occurring within 10 degrees of the
%                       plane dividing the hemifields
%
%     .rmsL             Lateral root mean squared error
%
%     .rmsP             Polar root mean squared error, omitting reversals
%                       and restricted to directions within 35 degrees of
%                       the vertical midline
%
%   
%   See also: demo_mclachlan2021 mclachlan2021
%
%   References:
%     D. Schonstein, L. Ferre, and B. F. G. Katz. Comparison of headphones
%     and equalization for virtual auditory source localization. Proc.
%     Euronoise, 2008.
%     
%     J. Reijniers, D. Vanderleist, C. Jin, C. S., and H. Peremans. An
%     ideal-observer model of human sound localization. Biological
%     Cybernetics, 108:169--181, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/mclachlan2021_metrics.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal M-Image
%   #Author: Michael Sattler
%   #Author: Roberto Barumerli
%   #Author: Glen McLachlan (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% real doa, in spherical and lat-pol coordinates
sph_real = SOFAconvertCoordinates(doa.real, 'cartesian', 'spherical');
[lat_real,pol_real]=sph2horpolar(sph_real(:,1),sph_real(:,2));
num_dirs = size(doa.real,1);
doa.est = squeeze(doa.est);

if size(doa.est, 3) ~= 1 % remove third dimension relative to repetition of experiments
    num_exp = size(doa.est, 2);
    doa_est = [];

    for i = 1:size(doa.est, 2)
        est = squeeze(doa.est(:,i,:));
        doa_est = [doa_est; est];
    end
else
    if(size(doa.est, 1) == num_dirs)
        num_exp = 1;
    else
        num_exp = size(doa.est, 1);
    end
    doa_est = doa.est;
end

% mean estimated doa over all experiments, in spherical and lat-pol coordinates
sph_est = SOFAconvertCoordinates(doa_est, 'cartesian', 'spherical');
[lat_est,pol_est] = sph2horpolar(sph_est(:,1),sph_est(:,2));

% matrix of all doas
m = zeros(size(sph_real, 1)*num_exp, 14);
m(:, 1:2) = repmat(sph_real(:, [1 2]), num_exp, 1);     % spherical real
m(:, 3:4) = sph_est(:, [1 2]);                          % spherical estimate
m(:, 5:6) = repmat([lat_real,pol_real], num_exp, 1);    % lat-pol real
m(:, 7:8) = [lat_est,pol_est];                          % lat-pol estimate
m(:, 9:11) = repmat(doa.real(:,1:3),num_exp, 1);        % cartesian real
m(:, 12:14) = doa_est(:,1:3);                           % cartesian estimate

% ensure real values, as sph2horpolar can return complex numbers 
% due to numerical approximations
m = real(m);

%% metrics

% lateral RMS
metric.rmsL=sqrt(sum((mynpi2pi(m(:,7)-m(:,5))).^2)/num_dirs);

% percentage of quadrant errors according to Schonstein et al. (2008)
idx_fb = find(acos(sum(m(:,9:11).*m(:,12:14),2))>... % fb reversals
    acos(sum([-m(:,9) m(:,10) m(:,11)].*m(:,12:14),2))...
    & (abs(m(:,6))<80 | (m(:,6)>100 & m(:,6)<260))); % drop dirs near midline

idx_ud = find(acos(sum(m(:,9:11).*m(:,12:14),2))>... % ud reversals
    acos(sum([m(:,9) m(:,10) -m(:,11)].*m(:,12:14),2))...
    & ((abs(m(:,6))>10 & m(:,6)<170)| m(:,6)>190)); % drop dirs near midline

idx = find(abs(m(:,5))<=35); % ignore 10 degrees around interaural axis
idx_fbm = intersect(idx,idx_fb);
idx_udm = intersect(idx,idx_ud); 

metric.mean_fbc = size(idx_fbm,1)*100/size(m,1); %percent fb reversals
metric.mean_udc = size(idx_udm,1)*100/size(m,1); %percent ud reversals

% polar RMS excluding quadrant errors and 10 degrees around interaural axis
m_rmsP = m;
m_rmsP([idx_fb; idx_ud],:) = [];
idx = abs(m_rmsP(:,5))<=35; % ignore 10 degrees around interaural axis
m_rmsP = m_rmsP(idx,:);

metric.rmsP=sqrt(sum((mynpi2pi(m_rmsP(:,8)-m_rmsP(:,6))).^2)/size(m_rmsP,1));

    % plottable data, omitting dirs near midline
    pol_err = abs(mynpi2pi(m(:,8)-m(:,6))); % plottable polar error
    lat_err = abs(mynpi2pi(m(:,7)-m(:,5))); % plottable lateral error
    p_err = sum(reshape(pol_err,num_dirs,num_exp),2)/num_exp;
    l_err = sum(reshape(lat_err,num_dirs,num_exp),2)/num_exp;
    n_fb = zeros(size(m,1),1);
    n_ud = zeros(size(m,1),1);
    n_fb(idx_fb) =  1;
    n_ud(idx_ud) =  1;
    metric.n_fbc = sum(reshape(n_fb,num_dirs,num_exp),2)*100/num_exp; % front-back
    metric.n_udc = sum(reshape(n_ud,num_dirs,num_exp),2)*100/num_exp; % up-down

    if isempty(varargin)
        error = 'none';
    else
        error = varargin{1};
    end

    switch error
        case 'polar'
            figure;plot_reijniers2014(doa.real,p_err); % plot polar
            title('Polar error in [deg]');
            caxis([0, 50]);
        case 'lateral'
            figure; plot_reijniers2014(doa.real,l_err); % plot lateral error
            title('Lateral error in [deg]');
        case 'fbc'
            figure;plot_reijniers2014(doa.real,metric.n_fbc); % plot fb errors
            title('Percentage of front-back reversals');
            caxis([0 50]);
        case 'udc'
            figure;plot_reijniers2014(doa.real,metric.n_udc); % plot ud errors
            title('Percentage of up-down reversals');
            caxis([0 50]);
        otherwise
    end
end

function out_deg=mynpi2pi(ang_deg)
ang=ang_deg/180*pi;
out_rad=sign(ang).*((abs(ang)/pi)-2*ceil(((abs(ang)/pi)-1)/2))*pi;
out_deg=out_rad*180/pi;
end

