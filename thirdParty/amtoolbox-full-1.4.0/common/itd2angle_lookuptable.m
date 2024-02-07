function lookup = itd2angle_lookuptable(hrtf,varargin)
%ITD2ANGLE_LOOKUPTABLE generates an ITD-azimuth lookup table for the given HRTF set
%   Usage: lookup = itd2angle_lookuptable(hrtf,fs,model);
%          lookup = itd2angle_lookuptable(hrtf,fs);
%          lookup = itd2angle_lookuptable(hrtf);
%
%   Input parameters:
%       hrtf   : HRTF data set (as SOFA file or struct)
%       fs     : sampling rate, (default: 44100) / Hz
%       model  : binaural model to use:
%
%                - dietz2011 uses the Dietz binaural model (default)
%                - lindemann1986 uses the Lindemann binaural model
%
%   Output parameters:
%       lookup : struct containing the polinomial fitting data for the
%                ITD -> azimuth transformation, p,MU,S, see help polyfit
%
%   ITD2ANGLE_LOOKUPTABLE(hrtf) creates a lookup table from the given HRTF data
%   set. This lookup table can be used by the dietz2011 or lindemann1986 binaural
%   models to predict the perceived direction of arrival of an auditory event.
%   The azimuth angle is stored in degree in the lookup table.
%
%   For the handling of the HRTF SOFA file format see
%   http://www.sofaconventions.org/
%
%   See also: dietz2011, lindemann1986, wierstorf2013
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592--605, 2011. [1]http ]
%     
%     H. Wierstorf, M. Geier, A. Raake, and S. Spors. A free database of
%     head-related impulse response measurements in the horizontal plane with
%     multiple distances. In Proceedings of the 130th Convention of the Audio
%     Engineering Society, 2011.
%     
%     H. Wierstorf, A. Raake, and S. Spors. Binaural assessment of
%     multi-channel reproduction. In J. Blauert, editor, The technology of
%     binaural listening, chapter 10. Springer, Berlin--Heidelberg--New York
%     NY, 2013.
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S016763931000097X
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/itd2angle_lookuptable.php


%   #Author: Hagen Wierstorf 
%   #Author: Robert Baumgartner (2017): SOFA compatibility

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% ===== Checking of input parameters ===================================
narginchk(1,3);

definput.flags.model = {'dietz2011','lindemann1986'};
definput.keyvals.fs = 44100;
[flags,kv]=ltfatarghelper({'fs'},definput,varargin);

% HRTF format
if isfield(hrtf,'GLOBAL_Conventions')
  flags.do_SOFA = true;
  flags.do_TUB = not(flags.do_SOFA);
  Obj = hrtf;
else
  flags.do_TUB = true;
  flags.do_SOFA = not(flags.do_TUB);
end

%% ===== Configuration ==================================================
% Samplingrate
fs = kv.fs;
% Time of noise used for the calculation (samples)
nsamples = fs;
% Noise type to use
noise_type = 'white';
% SFS Toolbox settings
conf.ir.useinterpolation = true;
conf.ir.useoriglength = true;
conf.fs = fs;
conf.c = 343;
conf.usefracdelay = false;


%% ===== Calculation ====================================================
% Generate noise signal
sig_noise = noise(nsamples,1,noise_type);

% get only the -90 to 90 degree part of the irs set
ele = 0;
idFrontHor = Obj.SourcePosition(:,2) == ele & ... % horizontal plane
    (((Obj.SourcePosition(:,1) >= -90 & Obj.SourcePosition(:,1) < 0) ...
    | Obj.SourcePosition(:,1) >= 270) | ... % front right
    (Obj.SourcePosition(:,1) >= 0 & Obj.SourcePosition(:,1) <= 90)); % front left
azi = Obj.SourcePosition(idFrontHor,1);
azi(azi>180) = azi(azi>180)-360;
% iterate over azimuth angles
nangles = length(azi);
% create an empty mod_itd, because the lindemann model didn't use it
mod_itd = [];

if flags.do_dietz2011

    itd = zeros(nangles,12);
    mod_itd = zeros(nangles,23);
    ild = zeros(nangles,23);
    for ii = 1:nangles
        % generate noise coming from the given direction
        sig = SOFAspat(sig_noise,hrtf,azi(ii),ele);
        % calculate binaural parameters
        [fine, cfreqs, ild_tmp, env] = dietz2011(sig,fs);
        % unwrap ITD
        itd_tmp = dietz2011_unwrapitd(fine.itd,ild_tmp(:,1:12),fine.f_inst,2.5);
        env_itd_tmp = dietz2011_unwrapitd(env.itd,ild_tmp(:,13:23),env.f_inst,2.5);
        % calculate the mean about time of the binaural parameters and store
        % them
        itd(ii,1:12) = median(itd_tmp,1);
        itd(ii,13:23) = median(env_itd_tmp,1);
        ild(ii,:) = median(ild_tmp,1);
    end

elseif flags.do_lindemann1986

    itd = zeros(nangles,36);
    ild = zeros(nangles,36);
    for ii = 1:nangles
        % generate noise coming from the given direction
        sig = SOFAspat(sig_noise,irs,azi(ii),ele);
        % Ten fold upsampling to have a smoother output
        %sig = resample(sig,10*fs,fs);
        % Calculate binaural parameters
        c_s = 0.3; % stationary inhibition
        w_f = 0; % monaural sensitivity
        M_f = 6; % decrease of monaural sensitivity
        T_int = inf; % integration time
        N_1 = 17640; % sample at which first cross-correlation is calculated
        [cc_tmp,dummy,ild(ii,:),cfreqs] = lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1);
        clear dummy;
        cc_tmp = squeeze(cc_tmp);
        % Calculate tau (delay line time) axes
        tau = linspace(-1,1,size(cc_tmp,1));
        % find max in cc
        for jj = 1:size(cc_tmp,2)
            [v,idx] = max(cc_tmp(:,jj));
            itd(ii,jj) = tau(idx)/1000;
        end
    end

end

% Fit the lookup data
for n = 1:size(itd,2)
    [p(:,n),S{n},MU(:,n)] = polyfit(itd(:,n),azi,12);
    [p_ild(:,n),S_ild{n},MU_ild(:,n)] = polyfit(ild(:,n),azi,12);
end
% Create lookup struct
lookup.p = p;
lookup.MU = MU;
lookup.S = S;
lookup.p_ild = p_ild;
lookup.MU_ild = MU_ild;
lookup.S_ild = S_ild;


