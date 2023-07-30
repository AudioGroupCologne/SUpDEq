function [data,itd] = data_pausch2022(varargin)
%DATA_PAUSCH2022 - Results from Pausch et al. (2022)
%
%   Usage: [data,itd] = data_pausch2022(type_stf,type_mod,fig,plot)
%
%   Input parameters:
%     type_stf        : flag that may be used to choose between one of the following:
%
%                       - hrtf load individual HRTF datasets
%                       - hartf_front load invidual front HARTF datasets
%                       - hartf_rear load invidual rear HARTF datasets
%
%   Output parameters:
%     data            : data struct
%
%                       - (if type_stf=={'hrtf','hartf_front','hartf_rear'}) 
%                         Individual spatial transfer functions, with fields id, the participant ID
%                         and sofa, the sofa file 
%                       - (if weights=='weights') 
%                         Polynomial regression weights (polynomial degree of P=4), 
%                         to be applied on a subset of M=5 individual features
%                       - if type_mod=={'kuhn','woodworth','woodworth_ext'}: 
%                         vector of (M*P+1) x 1 weights,
%                         if type_mod=='pausch': matrix of (M*P+1) x 4 weights
%                       - (if features=='features') 
%                         Individual features as published in Pausch et al. (2022)
%
%     itd            : itd struct (if type_stf=={'hrtf','hartf_front','hartf_rear'})
%                      Estimated ITDs for directions in the horizontal plane
%
%                      - itd_hor: ITDs in s [double]
%                      - itd_max: maximum ITD in s [s]
%                      - itd_arg_max_idx: index of the argument of the 
%                        the maximum ITD [double]
%                      - itd_arg_max_phi: argument of the the maximum ITD 
%                        in deg [double]
%                      - itd_hor_mean: direction-dependent mean ITDs
%                        across participants in s [double] 
%                      - itd_hor_std: direction-dependent standard
%                        deviation of ITDs across participants in s [double]
%                      - bp_fc_low: lower cut-off frequency (Hz) 
%                        of the bandpass filter applied 
%                        before estimating the ITDs [double]
%                      - bp_fc_high: upper cut-off frequency (Hz) 
%                        of the bandpass filter applied 
%                        before estimating the ITDs [double]
% 
%   The type_mod flag may be used to select the ITD model for the ITD 
%   predictions (only required if weights=='weights'):
%
%     'pausch'         hybrid ITD model by Pausch et al. (2022) (default)
%     'kuhn'           analytic ITD model by Kuhn (1977)
%     'woodworth'      analytic ITD model by Woodworth and Schlosberg (1954)
%     'woodworth_ext'  analytic ITD model by Woodworth and Schlosberg (1954), extended 
%                      by Aaronson and Hartmann (2014)
%
%   The weights flag may be one of the following:
%
%     'no_weights'     do not load polynomial regression weights (default)
%     'weights'        load polynomial regression weights
%
%   The features flag may be one of the following:
%
%     'no_features'    do not the individual features (default)
%     'features'       load the individual features
%
%   Additional key/value pairs include:
%
%     'bp_fc_low'      lower cut-off frequency (Hz) of the bandpass filter 
%                      applied before estimating the ITDs (default: []) [double]
%     'bp_fc_high'     upper cut-off frequency (Hz) of the bandpass filter 
%                      applied before estimating the ITDs (default: []) [double]
%
%
%   Requirements:
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics 
%      for Matlab (in e.g. thirdparty/SOFA)
%
%   2) Data in auxdata/pausch2022 (downloaded on the fly)
%
%
%   References:
%     F. Pausch, S. Doma, and J. Fels. Hybrid multi-harmonic model for the
%     prediction of interaural time differences in individual behind-the-ear
%     hearing-aid-related transfer functions. Acta Acust., 6:34, 2022.
%     [1]http ]
%     
%     G. F. Kuhn. Model for the interaural time differences in the azimuthal
%     plane. The Journal of the Acoustical Society of America,
%     62(1):157--167, 1977. [2]arXiv | [3]http ]
%     
%     R. S. Woodworth and H. Schlosberg. Experimental psychology, Rev. ed.
%     Holt, Oxford, England, 1954.
%     
%     N. L. Aaronson and W. M. Hartmann. Testing, correcting, and extending
%     the Woodworth model for interaural time difference. The Journal of the
%     Acoustical Society of America, 135(2):817--823, 2014. [4]arXiv |
%     [5]http ]
%     
%     References
%     
%     1. https://doi.org/10.1051/aacus/2022020
%     2. http://arxiv.org/abs/https://doi.org/10.1121/1.381498
%     3. https://doi.org/10.1121/1.381498
%     4. http://arxiv.org/abs/https://doi.org/10.1121/1.4861243
%     5. https://doi.org/10.1121/1.4861243
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_pausch2022.php


%   #Author: Florian Pausch (2022): Integration in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Parse flags and keyvals

definput.import={'pausch2022','amt_cache'};
definput.flags.weights = {'no_weights','weights'};
definput.flags.features = {'no_features','features'};

[flags,kv] = ltfatarghelper({},definput,varargin);

%% Load the polynomials regression weights

if flags.do_weights

    switch flags.type_mod
        case {'kuhn','woodworth','woodworth_ext'}
            if strcmp(flags.type_stf,'hrtf')
                [data,flag_mod] = amt_cache('get', 'cache_data_pausch2022_weights_hrtf_mod1_2plus');
            else % {hartf_front, hartf_rear}
                [data,flag_mod] = amt_cache('get', 'cache_data_pausch2022_weights_hartf_mod1_2plus');
            end
        otherwise
            if strcmp(flags.type_stf,'hrtf')
                [data,flag_mod] = amt_cache('get', 'cache_data_pausch2022_weights_hrtf_mod3');
            elseif strcmp(flags.type_stf,'hartf_front')
                [data,flag_mod] = amt_cache('get', 'cache_data_pausch2022_weights_hartf_front_mod3');
            else % hartf_rear
                [data,flag_mod] = amt_cache('get', 'cache_data_pausch2022_weights_hartf_rear_mod3');
            end
    end

    if isempty(flag_mod)

        if strcmp(flags.type_mod,'kuhn') || strcmp(flags.type_mod,'woodworth') || strcmp(flags.type_mod,'woodworth_ext')

            switch flags.type_stf
                case 'hrtf'
                        data = amt_load('pausch2022','weights.mat','weights_hrtf_mod1_2plus');

                case {'hartf_front','hartf_rear'}
                        data = amt_load('pausch2022','weights.mat','weights_hartf_mod1_2plus');
            end

        else % strcmp(flags.type_mod,'pausch')

            switch flags.type_stf
                case 'hrtf'
                        data = amt_load('pausch2022','weights.mat','weights_hrtf_mod3');

                case 'hartf_front'
                        data = amt_load('pausch2022','weights.mat','weights_fhartf_mod3');

                case 'hartf_rear'
                        data = amt_load('pausch2022','weights.mat','weights_rhartf_mod3');
            end

        end

        amt_disp([mfilename,': Loaded polynomial regression weights for model ',flags.type_mod,...
            ' (type_mod) and ',flags.type_stf,' (type_stf).'])

    end

    switch flags.type_mod
        case {'kuhn','woodworth','woodworth_ext'}
            if strcmp(flags.type_stf,'hrtf')
                amt_cache('set', 'cache_data_pausch2022_weights_hrtf_mod1_2plus',data,flags.type_mod);
            else % {hartf_front, hartf_rear}
                amt_cache('set', 'cache_data_pausch2022_weights_hartf_mod1_2plus',data,flags.type_mod);
            end
        otherwise
            if strcmp(flags.type_stf,'hrtf')
                amt_cache('set', 'cache_data_pausch2022_weights_hrtf_mod3',data,flags.type_mod);
            elseif strcmp(flags.type_stf,'hartf_front')
                amt_cache('set', 'cache_data_pausch2022_weights_hartf_front_mod3',data,flags.type_mod);
            else % hartf_rear
                amt_cache('set', 'cache_data_pausch2022_weights_hartf_rear_mod3',data,flags.type_mod);
            end
    end

end


%% Neither weights nor features: Load the individual directional datasets as per flag type_stf

if ~flags.do_weights && ~flags.do_features

    switch flags.type_stf
        case 'hrtf'
            [itd,flag_stf,bp_fc_low,bp_fc_high] = amt_cache('get', ...
                ['cache_data_pausch2022_hrtf_bp_',...
                num2str(kv.bp_fc_low),'_',num2str(kv.bp_fc_high),'_Hz'], flags.cachemode);
        case 'hartf_front'
            [itd,flag_stf,bp_fc_low,bp_fc_high] = amt_cache('get', ...
                ['cache_data_pausch2022_hartf_front_bp_',...
                num2str(kv.bp_fc_low),'_',num2str(kv.bp_fc_high),'_Hz'], flags.cachemode);
        otherwise
            [itd,flag_stf,bp_fc_low,bp_fc_high] = amt_cache('get', ...
                ['cache_data_pausch2022_hartf_rear_bp_',...
                num2str(kv.bp_fc_low),'_',num2str(kv.bp_fc_high),'_Hz'], flags.cachemode);
    end

    fileNames = struct();
    num_part = 30;
    phi_vec = 0:180;

    if flags.do_hrtf
        for id=1:num_part
            fileNames(id).name = ['id',num2str(id,'%02d'),'_HRTF.sofa'];
        end
    elseif flags.do_hartf_front
        for id=1:num_part
            fileNames(id).name = ['id',num2str(id,'%02d'),'_HARTF_front.sofa'];
        end
    elseif flags.do_hartf_rear
        for id=1:num_part
           fileNames(id).name = ['id',num2str(id,'%02d'),'_HARTF_rear.sofa'];
        end
    end

    amt_disp([mfilename,': Loading ',upper(flags.type_stf),' datasets...'])
    data = struct('id',{},'sofa',{});
    for idx = 1:length(fileNames)
        data(idx).id = idx;
        data(idx).sofa = amt_load('pausch2022',fileNames(idx).name);
    end
    amt_disp([mfilename,': Finished loading of ',upper(flags.type_stf),' datasets.'])

    if isempty(flag_stf) || ~isequal(bp_fc_low,kv.bp_fc_low) || ~isequal(bp_fc_high,kv.bp_fc_high)

        % Estimate the ITDs for directions in the horizontal plane, and calculate mu+/-std
        idx_dir_hor = data(1).sofa.SourcePosition(:,2)==0;
        num_dir_hor = sum(idx_dir_hor);

        itd = struct('itd_hor',{},'itd_max',{},'itd_arg_max_idx',{},...
            'itd_arg_max_phi',{},'itd_hor_mean',{},'itd_hor_std',{});
        itd(idx).itd_hor = zeros(num_dir_hor,1);

        amt_disp([mfilename,': Estimating ITDs between ',num2str(kv.bp_fc_low),...
            '-',num2str(kv.bp_fc_high),' Hz in ',upper(flags.type_stf),' datasets...'])
        itd_hor = zeros(num_dir_hor,num_part);
        fs = data(1).sofa.Data.SamplingRate;
        for idx=1:num_part
            data_hor = data(idx).sofa.Data.IR(idx_dir_hor,:,:);

            if idx==8 && strcmp(flags.type_stf,'hrtf')
                itd(idx).itd_hor = itd_MaxIACCe(data_hor,fs,...
                    kv.bp_fc_low,kv.bp_fc_high,1);
            else
                itd(idx).itd_hor = itd_MaxIACCe(data_hor,fs,...
                    kv.bp_fc_low,kv.bp_fc_high);
            end

            % store max(itd)
            [itd(idx).itd_max,itd(idx).itd_arg_max_idx] = max(itd(idx).itd_hor);

            % store arg_max_phi(itd)
            itd(idx).itd_arg_max_phi = phi_vec(itd(idx).itd_arg_max_idx);

            % temporarily store ITD data for the calculation of mean/std
            itd_hor(:,idx) = itd(idx).itd_hor;

        end

        [itd.itd_hor_mean] = deal(mean(itd_hor,2));
        [itd.itd_hor_std] = deal(std(itd_hor,0,2));

        amt_disp([mfilename,': Finished ITD estimations between ',num2str(kv.bp_fc_low),...
            '-',num2str(kv.bp_fc_high),' Hz in ',upper(flags.type_stf),' datasets.'])

        switch flags.type_stf
            case 'hrtf'
                amt_cache('set', ['cache_data_pausch2022_hrtf_bp_',...
                    num2str(kv.bp_fc_low),'_',num2str(kv.bp_fc_high),'_Hz'], ...
                    itd, flags.type_stf, kv.bp_fc_low, kv.bp_fc_high);
            case 'hartf_front'
                amt_cache('set', ['cache_data_pausch2022_hartf_front_bp_',...
                    num2str(kv.bp_fc_low),'_',num2str(kv.bp_fc_high),'_Hz'], ...
                    itd, flags.type_stf, kv.bp_fc_low, kv.bp_fc_high);
            otherwise
                amt_cache('set', ['cache_data_pausch2022_hartf_rear_bp_',...
                    num2str(kv.bp_fc_low),'_',num2str(kv.bp_fc_high),'_Hz'], ...
                    itd, flags.type_stf, kv.bp_fc_low, kv.bp_fc_high);
        end

    end

end

%% Load the individual features

if flags.do_features
    data = amt_cache('get', 'cache_data_pausch2022_features');
    if isempty(data)
          amt_disp('The features are not available when offline.')
    end
%         amt_disp([mfilename,': Downloading individual features...'])
%         if ~exist(fullfile(amt_auxdatapath,'pausch2022'),'dir'); mkdir(fullfile(amt_auxdatapath,'pausch2022')); end
%         websave(fullfile(amt_auxdatapath,'pausch2022','anthropometrics_hearing_aid_features.xlsx'),...
%             'https://publications.rwth-aachen.de/record/844866/files/anthropometrics_hearing_aid_features.xlsx');
%         amt_disp([mfilename,': Finished downloading individual features.'])
% 
%         [~,~,data] = xlsread(fullfile(amt_auxdatapath,'pausch2022','anthropometrics_hearing_aid_features.xlsx'),3);
% 
%         amt_cache('set', 'cache_data_pausch2022_features', data);
%     end
end

end

%% ------------------------------------------------------------------------
%  ---- INTERNAL FUNCTIONS ------------------------------------------------
%  ------------------------------------------------------------------------

function itd = itd_MaxIACCe(data,fs,bp_fc_low,bp_fc_high,rm_outliers)
%%ITD_MAXIACCE - Function to estimate interaural time differences (ITDs) 
%                based on the maximum of interaural cross correlation of 
%                energy envelopes in head-related impulse responses [5].
%                Parameters settings as used in [1].
%
%   Usage: itd = itd_MaxIACCe(data, phi_vec)
%
%   Input parameters (required):
%     data       : matrix containing the directional impulse responses (IRs),
%                  size(data) = num_dir_azimuth x 2 x filter length, that is 
%                  (144 x 2 x 256) [double]
%     fs         : sampling frequency (Hz) [double]
%     bp_fc_low  : lower cut-off frequency (Hz) of the bandpass filter 
%                  applied before estimating the ITDs (default: 500) [double]
%     bp_fc_high : upper cut-off frequency (Hz) of the bandpass filter 
%                  applied before estimating the ITDs (default: 1.5e3) [double]
%     rm_outl    : flag to remove ITD outliers in HRTF datasets measured
%                  from participant 8
%
%   Note: Parts of the code used in this function were adapted from the
%         ITA-Toolbox [6].
%
% [5] Areti Andreopoulou and Brian F. G. Katz, "Identification of perceptually 
%     relevant methods of inter-aural time difference estimation", The Journal 
%     of the Acoustical Society of America 142, 588-598 (2017) 
%     https://doi.org/10.1121/1.4996457 
% [6] Berzborn, M., Bomhardt, R., Klein, J., Richter, J. G., & Vorländer, M. 
%     (2017, March). The ITA-Toolbox: An open source MATLAB toolbox for acoustic 
%     measurements and signal processing. In 43th Annual German Congress on 
%     Acoustics, Kiel (Germany) (Vol. 6, pp. 222-225).

%   #Author: Florian Pausch, Institute for Hearing Technology and
%            Acoustics, RWTH Aachen University

if nargin<5
    rm_outliers=false;
end

%% apply low-pass filter on data
filter_order = 10;
f_cutoff = 1.6e3;

h_lp = fdesign.lowpass('n,f3db',filter_order,f_cutoff,fs);
Filter.Hd = design(h_lp,'butter');

data_filt = filter(Filter.Hd,data,3);

%% apply band-pass filter on pre-filtered data
if ~isempty(bp_fc_low) || ~isempty(bp_fc_high)
    if isempty(bp_fc_low)
        h_bp  = fdesign.lowpass('n,f3dB',filter_order,bp_fc_high,fs);
    elseif isempty(bp_fc_high)
        h_bp  = fdesign.bandpass('n,f3dB1,f3dB2',filter_order,bp_fc_low,fs/2,fs);
    else
        h_bp  = fdesign.bandpass('n,f3dB1,f3dB2',filter_order,bp_fc_low,bp_fc_high,fs);
    end
    Filter.Hd = design(h_bp,'butter');

    data_filt_itd = filter(Filter.Hd,data_filt,3);
else
    data_filt_itd = data_filt;
    amt_disp([mfilename,': No additional band-pass filtering applied before estimating the ITDs.'])
end

%% upsample and interpolate data
data_time_len = size(data,3)/fs;
data_time_vec = (1:size(data,3))./fs;

fac_upsample = 5;
fs_upsampled = fac_upsample*fs;
time_interp = 0:1/fs_upsampled:data_time_len;

data_time_interp = interp1(data_time_vec,permute(data_filt_itd,[3,1,2]),time_interp,'spline');
data_time_interp_e = data_time_interp.^2;

%% estimate ITDs in the horizontal plane
corr_ir = zeros(2*size(data_time_interp_e,1)-1,size(data_time_interp_e,2));
for idx_dir = 1:size(data_time_interp_e,2)
    corr_ir(:,idx_dir) = xcorr(data_time_interp_e(:,idx_dir,1),data_time_interp_e(:,idx_dir,2));
end

[~, idx_max] = max(corr_ir);
itd = ((numel(time_interp) - idx_max)/fs_upsampled).';

% remove ITD outliers in participant 8
if rm_outliers
    itd(18:21) = NaN;
    temp = itd;
    x = 1:numel(temp);
    itd(18:21) = pchip(x(~isnan(temp)), temp(~isnan(temp)), x(isnan(temp)));
end

end

% function data = local_dataload(instring)
% 
% for ii = 1:30
%   if numel(num2str(ii)) == 1
%     jj = ([num2str(0),num2str(ii)]);
%   else
%     jj=ii;
%   end
%   data = amt_load('pausch2022',['id',num2str(jj),'_', instring, '.sofa']);
% end
%   data = amt_load('pausch2022',['id',num2str(31),'a_', instring, '.sofa']);
%   data = amt_load('pausch2022',['id',num2str(31),'b_', instring, '.sofa']);
%   data = amt_load('pausch2022',['id',num2str(31),'c_', instring, '.sofa']);
%   data = amt_load('pausch2022',['id',num2str(31),'d_', instring, '.sofa']);
% end

