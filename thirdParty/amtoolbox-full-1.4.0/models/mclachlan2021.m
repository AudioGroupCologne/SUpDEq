function [doa, params] = mclachlan2021(template,target,varargin)
%MCLACHLAN2021 A dynamic ideal-observer model of human sound localization
%   Usage: [results,template,target] = mclachlan2014(template,target,'num_exp',20,'sig_S',4.2);
%
%   Input parameters:
%
%     template         : struct with the fields
%
%                        - fs : sampling rate (Hz)
%
%                        - fc : ERB frequency channels (Hz)
%
%                        - itd0 : itd computed for each hrir (samples)
%
%                        - H : Matrix containing absolute values of HRTFS for all grid points
%
%                        - coords : Matrix containing cartesian coordinates of all grid points, normed to radius 1m
%
%                        - T : angular template for each coordinate
%
%     target           : struct with the fields
%
%                        - fs : sampling rate
%
%                        - fc : ERB frequency channels
%
%                        - itd0 : itd corresponding to source position
%
%                        - S : sound source spectrum
%
%                        - H : Matrix containing absolute values of HRTFS for all source directions
%
%                        - coords : Matrix containing cartesian coordinates of all source positions to be estimated, normed to radius 1m
%
%                        - T : angular template for each coordinate
%
%   Output parameters:
%
%     doa                 : directions of arrival in spherical coordinates
%
%                           - est : estimated [num_sources, num_repetitions, 3]
%
%                           - real : actual [num_sources, 3]
%
%     params              : additional model's data computed for estimations
%
%                           - est_idx : Indices corresponding to template direction where the maximum probability density for each source position is found     	                  
%                           - est_loglik : Log-likelihood of each estimated direction
%                           
%                           - post_prob : Maximum posterior probability density for each target source
%
%                           - freq_channels : number of auditory channels
%
%                           - T_template : Struct with template data elaborated by the model
%
%                           - T_target : Struct with target data elaborated by the model
%
%
%   MCLACHLAN2021(...) is a dynamic ideal-observer model of human sound 
%   localization, by which we mean a model that performs optimal 
%   information processing within a Bayesian context. The model considers
%   all available spatial information contained within the acoustic
%   signals encoded by each ear over a specified hear rotation. Parameters 
%   for the optimal Bayesian model are determined based on psychoacoustic 
%   discrimination experiments on interaural time difference and sound 
%   intensity.
%
%
%   MCLACHLAN2021 accepts the following optional parameters:
%
%     'num_exp',num_exp   Set the number of localization trials. Default is num_exp = 500.
%
%     'SNR',SNR           Set the signal to noise ratio corresponding to
%                         different sound source intensities.
%                         Default value is SNR = 75 [dB]
%
%     'dt',dt             Time between each acoustic measurement in seconds.
%                         Default value is dt = 0.005.
%
%     'sig_itd0',sig      Set standard deviation for the noise on the initial
%                         itd. Default value is sig_itd0 = 0.569.
%
%     'sig_itdi',sig      Set standard deviation for the noise on the itd
%                         change per time step. Default value is sig_itdi = 1.
%
%     'sig_I',sig         Set standard deviation for the internal noise.
%                         Default value is sig_I = 3.5.     
%
%     'sig_S',sig         Set standard deviation for the variation on the 
%                         source spectrum. Default value is sig_S = 3.5.
%
%     'rot_type',type     Set rotation type. Options are 'yaw', 'pitch' and
%                         'roll'. Default value is 'yaw'.
%
%     'rot_size',size     Set rotation amount in degrees. Default value is 
%                         rot_size = 0.
%
%     'stim_dur',dur      Set stimulus duration in seconds. Default value is
%                         stim_dur = 0.1.
%
%   Further, cache flags (see amt_cache) can be specified.
%
%
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v1.1  or higher from 
%      http://sourceforge.net/projects/sofacoustics for Matlab (e.g. in 
%      thirdparty/SOFA)
%
%   See also: exp_reijniers2014 plot_reijniers2014
%   reijniers2014_metrics
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/mclachlan2021.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal M-Image
%   #Author: Glen McLachlan (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



%% Check input options
definput.import={'amt_cache'};
definput.import={'mclachlan2021'}; % looks for arg_mclachlan2021
% definput.importdefaults={'afb_osses2021','ihc_breebaart2001', 'adt_osses2021','mfb_jepsen2008'}; 

[flags,kv]  = ltfatarghelper({}, definput, varargin);

% prevent infinite matrix if rot_size==0                        
if kv.rot_size == 0
    kv.rot_size = 0.1;
end


%% derived parameters
rot_speed = kv.rot_size/kv.stim_dur;    % rotation speed
steps = floor(kv.stim_dur/kv.dt);       % amount of measurements steps during stimulus
alpha = linspace(0,kv.rot_size,steps);  % vector with angles at each measurement step                                               

%% define templates 
% T_template has size [Q x (2xfc+1)], where Q is number of sampled points
% and fc = number of frequency channels
T_template=[template.itd; template.ditd;...
    squeeze(template.H(:,:,1)-template.H(:,:,2)); ...
    0.5.*squeeze(template.H(:,:,1)+template.H(:,:,2))]'; 

T_target=[target.itd; target.ditd; ...
    squeeze(target.H(:,:,1)-target.H(:,:,2)); ...
    0.5.*squeeze(target.H(:,:,1)+target.H(:,:,2))]'; 


%% define error covariance matrix
sig_itd0 = kv.sig_itd0; %0.569;
sig_itdi = kv.sig_itdi; %still unknown, currently set to 1
sig_I = kv.sig_I; % 3.5; Intensity discrimination for broadband signal
sig_S = kv.sig_S; %3.5; Source's template error
sig = kv.sig; % Expected variance on the source strength - interchannel noise communication
sig_i = [sig_itd0, repmat(sig_itdi,1,steps-1)];   % var on ITD at each time step

% create M_beta covariance matrix
W = diag(1./(sig_i.^2));    % weight matrix
X = ones(steps,2);          % each column a slope of a parameter beta
X(:,2) = alpha;

M_beta = inv(X.'*W*X); % covariance matrix of ITD0 and ITDd

Sig = blkdiag(M_beta, 2*sig_I^2*eye(length(template.fc)), ...
    (sig_I^2/2 + sig_S^2)*eye(length(template.fc)) + sig^2); % full covariance matrix

%% simulate num_exp experimental trials 
num_exp = kv.num_exp;
invSig = inv(Sig);
num_targ = size(target.coords,1);
num_temp = size(template.coords,1);
log_lik = zeros(num_targ, num_exp);
doa_idx = zeros(num_targ, num_exp);
post_prob = zeros(num_targ, num_exp, num_temp);
doa_estimations = zeros(num_targ, num_exp, 3);
entropy = zeros(num_targ,1); % entropy in bits
Entropy = zeros(num_targ,1);
if nargout > 1
    X_all = zeros(num_targ, num_exp, size(T_target, 2));
end

for e = 1:num_exp
    amt_disp(sprintf('experiment %i', e))
    X = mvnrnd(T_target,Sig);
    if nargout > 1
        X_all(:,e,:) = X;
    end
    
    for s = 1:num_targ
        for d = 1:num_temp
            % Formula R
            u_diff = (X(s,:)-T_template(d,:));
            post_prob(s,e,d) = abs(exp(-0.5* u_diff*invSig*transpose(u_diff)));
        end
        % normalize
        post_prob(s,e,:) = post_prob(s,e,:)/sum(post_prob(s,e,:) + eps); 
        % maximum a posteriori
        [log_lik(s,e), doa_idx(s,e)] = max(post_prob(s,e,:));
        entropy(s)= - squeeze(post_prob(s,e,:))'*log2(squeeze(post_prob(s,e,:)) + eps);
        doa_estimations(s,e,:) = template.coords(doa_idx(s,e), :);
    end
    Entropy = Entropy + entropy; %cumulative entropy over experiments
end

Entropy = Entropy/num_exp;
Information = log2(num_temp) - Entropy; 

%% results
if (size(doa_estimations,1)==1)
    doa.est = squeeze(doa_estimations)';
else
    doa.est = doa_estimations;
end
doa.real = target.coords;

% user required more than the estimations
if nargout > 1
    params.template_coords = template.coords;
    params.post_prob = post_prob;
    params.entropy = Entropy;
    params.information = Information;
    params.est_idx = doa_idx;
    params.est_loglik = log_lik;
    params.X = X_all;
    params.T_template = T_template;
    params.T_target = T_target;
    params.freq_channels = template.fc;
    params.Tidx.itd = 1;
    assert(length(target.fc)==length(template.fc))
    params.Tidx.Hp = params.Tidx.itd + (1:length(target.fc));
    params.Tidx.Hm = params.Tidx.Hp(end) + (1:length(target.fc));
else
    clear X_all post_prob doa_idx log_lik
end

end

function Y_N = SH(N, dirs)
% calculate spherical harmonics up to order N for directions dirs [azi ele;...] (in radiant)
% 
    N_dirs = size(dirs, 1);
    N_SH = (N+1)^2;
	dirs(:,2) = pi/2 - dirs(:,2); % convert to inclinations

    Y_N = zeros(N_SH, N_dirs);

	  % n = 0
	Lnm = legendre(0, cos(dirs(:,2)'));
	Nnm = sqrt(1./(4*pi)) * ones(1,N_dirs);
	CosSin = zeros(1,N_dirs);
	CosSin(1,:) = ones(1,size(dirs,1));
	Y_N(1, :) = Nnm .* Lnm .* CosSin;
	
	  % n > 0
	idx = 1;
    for n=1:N
        
        m = (0:n)';            

		Lnm = legendre(n, cos(dirs(:,2)'));
		condon = (-1).^[m(end:-1:2);m] * ones(1,N_dirs);
		Lnm = condon .* [Lnm(end:-1:2, :); Lnm];
		
		mag = sqrt( (2*n+1)*factorial(n-m) ./ (4*pi*factorial(n+m)) );
		Nnm = mag * ones(1,N_dirs);
		Nnm = [Nnm(end:-1:2, :); Nnm];
		
		CosSin = zeros(2*n+1,N_dirs);
			% m=0
		CosSin(n+1,:) = ones(1,size(dirs,1));
			% m>0
		CosSin(m(2:end)+n+1,:) = sqrt(2)*cos(m(2:end)*dirs(:,1)');
			% m<0
		CosSin(-m(end:-1:2)+n+1,:) = sqrt(2)*sin(m(end:-1:2)*dirs(:,1)');

		Ynm = Nnm .* Lnm .* CosSin;
        Y_N(idx+1:idx+(2*n+1), :) = Ynm;
        idx = idx + 2*n+1;
    end
    
    Y_N = Y_N.';
    
end


