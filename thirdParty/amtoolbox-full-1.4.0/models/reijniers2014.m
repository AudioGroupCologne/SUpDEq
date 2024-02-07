function [doa, params] = reijniers2014(template,target,varargin)
%REIJNIERS2014   Bayesian spherical sound localization model (basic)
%
%   Usage: [results,template,target] = reijniers2014(template,target,'num_exp',20,'sig_S',4.2);
%
%   Input parameters:
%     template       : struct, with the fields
%                      .fs      sampling rate (Hz)
%                      .fc      ERB frequency channels (Hz)
%                      .itd     itd computed for each hrir (samples)
%                      .H       Matrix containing absolute values of HRTFS for all grid points
%                      .coords  Matrix containing cartesian coordinates of all grid points, normed to radius 1m
%                      .T       angular template for each coordinate
%     target        : struct, with the fields
%                     .fs       sampling rate
%                     .fc       ERB frequency channels
%                     .itd      itd corresponding to source position
%                     .S        sound source spectrum
%                     .H        Matrix containing absolute values of HRTFS for all source directions
%                     .coords   Matrix containing cartesian coordinates of all source positions to be estimated, normed to radius 1m
%                     .T        angular template for each coordinate
%
%
%   Output parameters:
%     doa               : directions of arrival in spherical coordinates
%     params            : additional model's data computed for estimations
%
%
%   'doa' contains the following fields:
%
%     .est            estimated [num_sources, num_repetitions, 3]
%     .real           actual    [num_sources, 3]
%
%   'params' contains the following fields:
%
%     .est_idx        : Indices corresponding to template direction where the maximum probability density for each source position is found
%     .est_loglik     : Log-likelihood of each estimated direction
%     .post_prob      : Maximum posterior probability density for each target source
%     .freq_channels  : number of auditory channels
%     .T_template     : Struct with template data elaborated by the model
%     .T_target       : Struct with target data elaborated by the model
%     .Tidx           : Helper with indexes to parse the features from T and X
%
%
%   REIJNIERS2014 accepts the following optional parameters:
%
%     'num_exp',num_exp   Set the number of localization trials. Default is num_exp = 500.
%
%     'SNR',SNR           Set the signal to noise ratio corresponding to
%                         different sound source intensities.
%                         Default value is SNR = 75 [dB]
%
%     'sig_itd',sig       Set standard deviation for the noise on the itd.
%                         Default value is sig_itd = 0.569.
%
%     'sig_I',sig         Set standard deviation for the internal noise.
%                         Default value is sig_I = 3.5.     
%
%     'sig_S',sig         Set standard deviation for the variation on the 
%                         source spectrum. Default value is sig_I = 3.5.
%
%
%   Further, cache flags (see amt_cache) can be specified.
%
%
%   Description: 
%   ------------
%
%   REIJNIERS2014(...) is an ideal-observer model of human sound 
%   localization, by which we mean a model that performs optimal 
%   information processing within a Bayesian context. The model considers
%   all available spatial information contained within the acoustic
%   signals encoded by each ear. Parameters for the optimal Bayesian model
%   are determined based on psychoacoustic discrimination experiments on 
%   interaural time difference and sound intensity.
%
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v1.1  or higher from 
%      http://sourceforge.net/projects/sofacoustics for Matlab (e.g. in 
%      thirdparty/SOFA)
%
%   See also: exp_reijniers2014 plot_reijniers2014 reijniers2014_featureextraction
%   reijniers2014_metrics
%   exp_engel2021
%   baumgartner2014
%   mclachlan2021
%   baumgartner2013
%   demo_reijniers2014
%   demo_mclachlan2021
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/reijniers2014.php


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




%% Check input options
definput.import={'amt_cache'};
definput.flags.disp = {'no_debug','debug'};

definput.flags.type = {'fig2'};

definput.keyvals.num_exp = 500;
definput.keyvals.SNR = 75;

% parameters of the model computed in the supplementary material
definput.keyvals.sig_itd = 0.569;
definput.keyvals.sig_I = 3.5;
definput.keyvals.sig_S = 3.5; 
definput.keyvals.sig = 5; 

[flags,kv]  = ltfatarghelper({'num_exp','SNR','sig_itd','sig_I','sig_S', 'sig'}, ...
                             definput, varargin);
        
%% sample uniformly over sphere with N is number of directions
% NOTE: amt_load('reijniers2014', 'dirs.mat') contains the sampled point on a unitary
% sphere
dirs=amt_load('reijniers2014','dirs.mat');
dirs=dirs.dirs;
if(isempty(dirs))
    error('New directions grid not available. Please check your internet connection!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   if there is the need to generate a different number of directions
%   this toolbox is required and following commented code need to be
%   executed
%   S2-Sampling-Toolbox-master V. 79cc337 
%   from 13 June 2019 or higher from 
%   https://github.com/AntonSemechko/S2-Sampling-Toolbox
% if isempty(dirs)
%     num_dirs = 2000; 
%     [dirs,~,~,~] = ParticleSampleSphere('N',num_dirs); 
%     save('AUX DIRECTORY/reijniers2014/dirs.mat','dirs.mat',dirs);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% remove the points from the unitary sphere below HRTF lowest elevation
idx = find(dirs(:,3) > min(template.coords(:, 3))); 
dirs = dirs(idx,:);
num_dirs = length(idx); 


%% interpolate at uniformly distributed directions and update H and itd

% calculate spherical harmonic coefficients of H and itd, using tikonov regularization
SHorder = 15; % spherical harmonic order 

[AZ,EL] = cart2sph(template.coords(:,1),template.coords(:,2),template.coords(:,3));
Y_N = SH(SHorder, [AZ EL]); 

% tikonov
lambda = 4.;
SIG = eye((SHorder+1)^2);
SIG(1:(2+1)^2,1:(2+1)^2) = 0;

cH(:,:,1) = transpose((Y_N'*Y_N+lambda*SIG)\Y_N'*squeeze(template.H(:,:,1))');
cH(:,:,2) = transpose((Y_N'*Y_N+lambda*SIG)\Y_N'*squeeze(template.H(:,:,2))');
citd = (Y_N'*Y_N+lambda*SIG)\Y_N'*template.itd(:);

% interpolate at uniformly distributed directions and update
[AZ,EL] = cart2sph(dirs(:,1),dirs(:,2),dirs(:,3));
Y_N = SH(SHorder, [AZ EL]); 
template.H = [];
template.itd = [];
template.H(:,:,1) = transpose(Y_N*squeeze(cH(:,:,1))');
template.H(:,:,2) = transpose(Y_N*squeeze(cH(:,:,2))');  
template.itd = Y_N*citd;
template.coords = dirs;


%% transform HRTF and itd to perceptually relevant units 
% itd trasformation through jnd - see supplementary material for parameters
a = 32.5e-6;
b = 0.095;
template.itd = sign(template.itd) .* ((log(a + b * abs(template.itd)) - log(a)) / b); 
target.itd = sign(target.itd) .* ((log(a + b * abs(target.itd)) - log(a)) / b); 

% account for SNR and frequency-dependent hearing sensitivity (see section 2.1 in SI)
% add source spectrum to target and to template 
% see last formula in the supplementary materials
temp_H = template.H + repmat(target.S(:), 1, size(template.H, 2), 2);
targ_H = target.H + repmat(target.S(:), 1, size(target.H, 2), 2);

SNR = kv.SNR; % defined as maximal SNR (in interval 2kHz-7kHz)

temp_H = max(temp_H ,-SNR);
temp_H(template.fc<=2000,:,:) = max(temp_H(template.fc<=2000,:,:),-SNR + 10);
temp_H(template.fc>=7000,:,:) = max(temp_H(template.fc>=7000,:,:),-SNR + 20);

targ_H = max(targ_H,-SNR);
targ_H(target.fc<=2000,:,:) = max(targ_H(target.fc<=2000,:,:),-SNR + 10);
targ_H(target.fc>=7000,:,:) = max(targ_H(target.fc>=7000,:,:),-SNR + 20);

%% define templates 
T_template=[template.itd, ...
    squeeze(temp_H(:,:,1)-temp_H(:,:,2))', ...
    0.5.*squeeze(temp_H(:,:,1)+temp_H(:,:,2))']; 

T_target=[target.itd, ...
    squeeze(targ_H(:,:,1)-targ_H(:,:,2))', ...
    0.5.*squeeze(targ_H(:,:,1)+targ_H(:,:,2))']; 


%% define covariance matrix
sig_itd = kv.sig_itd; %0.569;
sig_I = kv.sig_I; % 3.5; Intensity discrimination for broadband signal
sig_S = kv.sig_S; %3.5; Source's template error
sig = kv.sig; % Expected variance on the source strength - interchannel noise communication
Sig = blkdiag(sig_itd^2, 2*sig_I^2*eye(length(template.fc)), ((sig_I^2)/2 + sig_S^2)*eye(length(template.fc)) + sig^2);

%% simulate num_exp experimental trials 
num_exp = kv.num_exp;
invSig = inv(Sig);
num_src = size(target.coords,1);
log_lik = zeros(num_src, num_exp);
doa_idx = zeros(num_src, num_exp);
post_prob = zeros(num_src, num_exp, num_dirs);
doa_estimations = zeros(num_src, num_exp, size(template.coords, 2));
if nargout > 1
    X_all = zeros(num_src, num_exp, size(T_target, 2));
end

for e = 1:num_exp
    X = mvnrnd(T_target,Sig);
    if nargout > 1
        X_all(:,e,:) = X;
    end
    
    for s = 1:num_src
        for d = 1:num_dirs
            % Formula R
            u_diff = (X(s,:)-T_template(d,:));
            post_prob(s,e,d) = abs(exp(-0.5* u_diff*invSig*transpose(u_diff)));
        end
        % normalize
        post_prob(s,e,:) = post_prob(s,e,:)/sum(post_prob(s,e,:) + eps); 
        % maximum a posteriori
        [log_lik(s,e), doa_idx(s,e)] = max(post_prob(s,e,:));    
        doa_estimations(s,e,:) = template.coords(doa_idx(s,e), :);
    end
end

%% results
doa.est = doa_estimations;
doa.real = target.coords;

% user required more than the estimations
if nargout > 1
    params.template_coords = template.coords;
    params.post_prob = post_prob;
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


