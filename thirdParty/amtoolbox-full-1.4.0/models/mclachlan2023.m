function [doa, params] = mclachlan2023(template, varargin)
%MCLACHLAN2023 - A dynamic ideal-observer model of human sound localisation
%   Usage: [results,template,target] = mclachlan2023(template,target,'num_exp',20,'sig_S',4.2);
%
%   Input parameters:
%
%     template.fs      : sampling rate (Hz)
%     template.fc      : ERB frequency channels (Hz)
%     template.itd0    : itd computed for each hrir (samples)
%     template.H       : Matrix containing absolute values of HRTFS for all grid points
%     template.coords  : Matrix containing cartesian coordinates of all grid points, normed to radius 1m
%     template.T       : angular template for each coordinate
%
%     'theta',theta  : Set pitch rotation amount in rad. Default value is 
%                      theta = 0.
%     'phi',phi      : Set yaw rotation amount in rad. Default value is 
%                      phi = 0.01/180*pi.
%
%   Output parameters:
%
%     doa               : directions of arrival in spherical coordinates
%
%       .est            : estimated [num_sources, num_repetitions, 3]
%       .real           : actual    [num_sources, 3]
%
%     params            : additional model's data computerd for estimations
%
%       .est_idx        : Indices corresponding to template direction where
%                         the maximum probability density for each source
%                         position is found
%
%     	.est_loglik     : Log-likelihood of each estimated direction
%
%       .post_prob      : Maximum posterior probability
%
%                         density for each target source
%
%       .freq_channels  : number of auditory channels
%
%       .T_template     : Struct with template data elaborated by the model
%
%       .T_target       : Struct with target data elaborated by the model
%
%     	.Tidx           : Helper with indexes to parse
%                         the features from T and X
%
%   MCLACHLAN2023 accepts the following optional parameters:
%
%     'num_exp',num_exp Set the number of localization trials.
%                       Default is num_exp = 500.
%
%     'dt',dt           Time between each acoustic measurement in seconds.
%                       I.e. time step. Default value is dt = 0.005.
%
%     'sig_itd0',sig    Set standard deviation for the noise on the initial
%                       itd. Default value is sig_itd0 = 0.569.
%
%     'sig_itdi',sig    Set standard deviation for the noise on the itd
%                       change per time step. Default value is sig_itdi = 1.
%
%     'sig_I',sig       Set standard deviation for the internal noise.
%                       Default value is sig_I = 3.5.     
%
%     'sig_S',sig       Set standard deviation for the variation on the 
%                       source spectrum. Default value is sig_S = 3.5.
%
%     'stim_dur',dur    Set stimulus duration in seconds. Default value is
%                       stim_dur = 0.1.
%
%   Further, cache flags (see amt_cache) can be specified.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.10.0/doc/models/mclachlan2021.php
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/mclachlan2023.php


% Copyright (C) 2009-2021 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.10.0
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

%
%   Description: 
%   ------------
%
%   `mclachlan2021(...)` is a dynamic ideal-observer model of human sound 
%   localization, by which we mean a model that simulates head rotations to 
%   perform optimal information processing within a Bayesian context. The 
%   model considers all available spatial information contained within the 
%   acoustic signals encoded by each ear over a specified hear rotation. 
%   Parameters for the optimal Bayesian model are determined based on 
%   psychoacoustic discrimination experiments on interaural time difference 
%   and sound intensity.
%
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v1.1  or higher from 
%      http://sourceforge.net/projects/sofacoustics for Matlab (e.g. in 
%      thirdparty/SOFA)
%
%   See also: exp_reijniers2014 plot_reijniers2014 reijniers2014_preproc
%   reijniers2014_metrics
%
%   References: reijniers2014 mclachlan2021

%   AUTHOR: Glen McLachlan and  Herbert Peremans



%% Check input options
%c = parcluster; parpool(c)
% sc = parallel.pool.Constant(RandStream('Threefry'))
definput.import={'amt_cache'};
definput.import={'mclachlan2023'};
definput.flags.type = {'fig2'};

[flags,kv]  = ltfatarghelper({}, definput, varargin);

SHorder = kv.SHorder;
template_coords = template.coords;
tempx=template_coords(1,:); tempy=template_coords(2,:); tempz=template_coords(3,:);
target_coords=load('dirs_1527.mat');
target_coords=target_coords.dirs';
c_itd = template.coef.itd;
c_H1 = template.coef.H(:,:,1);
c_H2 = template.coef.H(:,:,2);

%% misc parameters
M = floor(kv.stim_dur/kv.dt)+1;         % amount of measurements steps during stimulus
tmeas=(0:kv.dt:(M-1)*kv.dt);            % vector of time samples of measurement sequence
num_targets = size(target_coords,2);    % amount of test target directions
num_template = size(template_coords,2); % amount of template directions
num_fc=length(template.fc);

if strcmp(kv.rot_type,'yaw')
    u = [kv.rot_size*tmeas(2)/max(tmeas),0]; %motor command at each time step
    u_dir = [1,1];
elseif strcmp(kv.rot_type,'pitch')
    u = [0,kv.rot_size*tmeas(2)/max(tmeas)];
    u_dir = [1,1];
end


%% head orientation list for marginalisation
% range of considered initial head orientations
h0range=12;
azmin_h=-h0range; azmax_h=h0range;
elmin_h=-h0range; elmax_h=h0range;

% uniform distribution of intitial head orientations within the allowed range
[AZ,EL] = meshgrid(azmin_h:4:azmax_h,elmin_h:4:elmax_h);
tmp=cat(2,AZ',EL');
tmp=reshape(tmp,[],2);
list_dir0_h=tmp/180*pi;
[list_dir0_hc(:,1),list_dir0_hc(:,2),list_dir0_hc(:,3)] = sph2cart(list_dir0_h(:,1),list_dir0_h(:,2),ones(size(list_dir0_h,1),1));

% noise on head rotation angle sequence
sig_head2 = zeros(1,M);
sig_head2(1) = kv.sig_H^2;

for t=1:M-1
K(t) = (sig_head2(t) + kv.sig_u^2) / (sig_head2(t) + kv.sig_u^2 + kv.sig_H^2); %Kalman filter
sig_head2(t+1)=(1-K(t))*(sig_head2(t)+kv.sig_u^2);
end

sig_itd2=[(kv.sig_itd0)^2 ; kv.sig_itdi^2*ones(M-1,1)]; %noise on ITD sequence 

%define covar matrix for noise on spectral content
covar_spec(:,:,1) = blkdiag(2*kv.sig_I0^2*eye(num_fc), ...
     (kv.sig_I0^2/2 + kv.sig_S^2)*eye(num_fc) + kv.sig^2);
for t=1:M
covar_spec(:,:,t) = blkdiag(2*kv.sig_Ii^2*eye(num_fc), ...
     (kv.sig_Ii^2/2 + kv.sig_S^2)*eye(num_fc) + kv.sig^2);
end

az0=kv.az0;
el0=kv.el0;

%% posterior computation
num_dir0_h = size(list_dir0_h,1); %amount of initial head orientations

tic
for e=1:kv.num_exp
%     stream = sc.Value;        % Extract the stream from the Constant
%     stream.Substream = e;
%     % compute prior
    [AZ, EL] = cart2sph(tempx,tempy,tempz);
    prior_targdir = mclachlan2022_prior([AZ',EL'],'full'); %'horizon'
    prior_targdir = repmat(prior_targdir,1,num_targets);
%     % check prior
%     figure; scatter(EL,prior_targdir)
%     figure; scatter3(template_coords(:,1),template_coords(:,2),template_coords(:,3),prior_targdir)

    Theta_H=[az0,el0]; %initial head position
    y_H=zeros(M+1,2);
    y_H(1,:)=Theta_H+chol(kv.sig_H^2)*randn([1,2]).*u_dir; %first head position observation
    mu_H=y_H; %initiate head position estimate

    for t=1:M % per time step
        jointdistr_marg=zeros([num_template, num_targets]); %stored posterior probability per target
        
        tmp_TM = mclachlan2023_rotatedirs(Theta_H(t,1)/180*pi,Theta_H(t,2)/180*pi);% rotated source directions
        tmp_TM_inv=transpose(tmp_TM);   %world to head coords
        tmp_T = tmp_TM_inv*target_coords;

        [AZ,EL]= cart2sph(tmp_T(1,:),tmp_T(2,:),tmp_T(3,:)); %true source direction sequence
        dirs_seqT=[AZ;EL];

        % noise on acoustic measurements
        delta_itd=chol(sig_itd2(t))*randn([1,num_targets]); % noise on itd measurement
        delta_spec=chol(covar_spec(:,:,t))*randn([2*num_fc,num_targets]); % noise on spectral measurement

        %noise on sensorimotor information
        delta_H=chol(kv.sig_H^2)*randn([1,2]); %noise on measurement
        delta_u=chol(kv.sig_u^2)*randn([1,2]); %noise on motor control around axis of rotation
	if t<M
        Theta_H(t+1,:)=Theta_H(t,:)+u+delta_u.*u_dir;
        y_H(t+1,:)=Theta_H(t+1,:)+delta_H.*u_dir; %[M,2] noisy measurement of head orientations
	    mu_H(t+1,:)=(1-K(t))*(mu_H(t,:)+u)+K(t)*y_H(t+1,:); % estimate of true head orientation using Kalman filter
	end

        %noisy spectral measurements
        Y_N = SH(SHorder, dirs_seqT'); 
        X_H1 = Y_N*c_H1';
        X_H2 = Y_N*c_H2';
        X_spec = [(X_H1-X_H2)'; 0.5*(X_H1+X_H2)']; %transform to X-,X+
        y_spec=X_spec+delta_spec; %[2*nchan,ntarg] noisy measurement

        %noisy itd measurements
        X_itd = (Y_N*c_itd)';
        y_itd=(X_itd+delta_itd)'; %[M,ntarg] noisy measurement

        % rotate head orientations around time step
        rotM = mclachlan2023_rotatedirs(mu_H(t,1)/180*pi,mu_H(t,2)/180*pi);
        list_dir_hc=rotM*list_dir0_hc';
        [az,el]=cart2sph(list_dir_hc(1,:),list_dir_hc(2,:),list_dir_hc(3,:));
        list_dir_h=[az;el];

        for dirH_indx=1:num_dir0_h % marginalisation over all head orientations    

            dir0_h=list_dir_h(:,dirH_indx);    %head orientation [rad]
    
            %calculate rotation matrix corresponding with initial head orientation
            rotM = mclachlan2023_rotatedirs(dir0_h(1),dir0_h(2)); %head to world coords
            rotM_inv=transpose(rotM);   %world to head coords
        
            %transform template directions into head reference system:
            dirsT_tmp = (rotM_inv*template_coords)';
            [dirsT_az,dirsT_el] = cart2sph(dirsT_tmp(:,1),dirsT_tmp(:,2),dirsT_tmp(:,3));
    
            %interpolate acoustic cues to initial head orientation
            Y_N = SH(SHorder, [dirsT_az dirsT_el]); 
            itd_tmp = Y_N*c_itd;     %ITD in head reference system
            likelihood_itd=exp(-0.5*pdist2(itd_tmp,y_itd,'mahalanobis',sig_itd2(t)).^2); %ITD likelihood of current timestep
            
            spec_tmp = [transpose(Y_N*c_H1');...
                        transpose(Y_N*c_H2')];
            spec_tmp = [squeeze(spec_tmp(1:end/2,:)-spec_tmp(end/2+1:end,:))', ...
                0.5.*squeeze(spec_tmp(1:end/2,:)+spec_tmp(end/2+1:end,:))']';
            likelihood_spectral=exp(-0.5*pdist2(spec_tmp',y_spec','mahalanobis',covar_spec(:,:,t)).^2); %spectral likelihood
            likelihood_acoustic=likelihood_spectral.*likelihood_itd;  

            likelihood_acoustic = likelihood_acoustic ./ (sum(likelihood_acoustic, 1)); % normalise
%             makefig4(template.coords',likelihood_acoustic);
    
            % head postition likelihood with independent azimuth and elevation angle
            tmp=exp(-0.5*pdist2(mu_H(t,1),dir0_h(1)*180/pi,'mahalanobis',sig_head2(t)).^2);
            likelihood_motor=tmp*exp(-0.5*pdist2(mu_H(t,2),dir0_h(2)*180/pi,'mahalanobis',sig_head2(t)).^2);
    
            posterior_motor = likelihood_motor; %Bayes' rule, multiply with prior if it isn't uniform!!
    
            jointdistr = posterior_motor.*likelihood_acoustic; %joint PDF of acoustic & sensorimotor PDF
    
            %marginal joint PDF for each target direction: [n_temp, n_targ*n_exp]
            jointdistr_marg=jointdistr_marg+jointdistr;
            %figure; plot_reijniers2014(template_coords',jointdistr_marg(:,4));
        end %initial head direction list
   
        %update posterior
        jointdistr_marg=jointdistr_marg./(sum(jointdistr_marg,1)); 
        posterior_targdir=prior_targdir.*jointdistr_marg; 
        posterior_targdir = posterior_targdir./(sum(posterior_targdir,1)); %normalise
        prior_targdir = posterior_targdir; %update prior to posterior from t-1
        prior_targdir=prior_targdir./(sum(prior_targdir,1)); 
%         figure; plot_reijniers2014(template_coords',jointdistr_marg(:,500),[],target_coords(:,500)');
%         caxis([0,3e-3])
%         figure; plot_reijniers2014(template_coords',prior_targdir(:,500),[],target_coords(:,500)');
%         caxis([0,0.03])
    end %time step

    post_prob(:,:,e) = posterior_targdir;

    % maximum a posteriori
    [log_lik(:,e), doa_idx(:,e)] = max(posterior_targdir,[],1); 
end %experiment list
toc
post_prob=squeeze(mean(post_prob,3));
%estimated directions of arrival in cart coords per target per trial
doa_estimations = reshape(template_coords(:,doa_idx),[3,num_targets,kv.num_exp]);

% entropy & information
% entropy = - sum(sum(post_prob.*log2(post_prob+eps),3),2);
% 
% entropy = entropy/kv.num_exp;
% information = log2(num_template) - entropy;

%% results
doa.est = permute(doa_estimations,[2,3,1]);  % estimated doa [ntargets x nexp x 3]
doa.real = target_coords';   % reasl doa [ntargets x 3]

% user required more than the estimations
if nargout > 1
    params.template_coords = template_coords;   % template coordinates
    params.post_prob = post_prob;               % posterior probabilities [ncoords x nexp]
    %params.information = Information;
    params.est_idx = doa_idx;                   % index of estimated doa [ntargets x nexp]
    params.est_loglik = log_lik;                % likelihood of estimated doa [ntargets x nexp]
    params.freq_channels = template.fc;
    params.Tidx.itd = 1;
    params.Tidx.Hp = params.Tidx.itd + (1:num_fc);
    params.Tidx.Hm = params.Tidx.Hp(end) + (1:num_fc);
    %params.entropy = entropy;
    %params.information = information;
    params.sig.itd0 = kv.sig_itd0; params.sig.itdi = kv.sig_itdi; params.sig.I0 = kv.sig_I0; params.sig.Ii = kv.sig_Ii; params.sig.S = kv.sig_S; params.sig.b = kv.sig; params.sig.H = kv.sig_H; params.sig.u = kv.sig_u;
else
    clear post_prob doa_idx log_lik
end
end
 
function rotM = mclachlan2023_rotatedirs(az,el)
    % return rotation matrix for specified azimuth and elevation rotation
    
    rV_az=-az*[0,0,1];        % positive rotation is to the right
    rM_az = rotationVectorToMatrix(rV_az);
    rV_el=el*[0,1,0];        % positive rotation is upwards
    rM_el = rotationVectorToMatrix(rV_el);
    rotM=rM_az*rM_el;           %head to world coords
end
