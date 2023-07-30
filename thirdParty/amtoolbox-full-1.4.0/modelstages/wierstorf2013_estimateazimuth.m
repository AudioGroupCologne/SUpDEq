function [phi,phi_std,itd,ild,cfreqs] = wierstorf2013_estimateazimuth(insig,lookup,varargin)
%WIERSTORF2013_ESTIMATEAZIMUTH Estimate the perceived azimuth using a binaural model
%   Usage: [phi,phi_std,itd,ild,cfreqs] = wierstorf2013_estimateazimuth(sig,lookup,'fs',44100,'dietz2011')
%          [phi,phi_std,itd,ild,cfreqs] = wierstorf2013_estimateazimuth(sig,lookup)
%
%   Input parameters:
%       sig                   : binaural singal
%       lookup                : lookup table to map ITDs to angles (struct)
%
%   Output parameters:
%       phi     : estimated azimuth / deg
%       phi_std : standard deviation of the estimated azimuth / deg
%       itd     : calculated ITD (s)
%       ild     : calculated ILD (dB)
%       cfreqs  : center frequencies of used auditory filters (Hz)
%
%   WIERSTORF2013_ESTIMATEAZIMUTH(sig,lookup) uses a binaural model to
%   estimate the perceived direction for a given binaural signal.  Therefore,
%   it needs the struct lookup, which maps ITD values to the corresponding
%   angles. This can be created with the ITD2ANGLE_LOOKUPTABLE function.
%   The azimuth values are first calculated for every frequency channel and
%   after that their median is calculated. In this process the different
%   frequency channels could be weighted and outlier could be removed, see the
%   options below. The default setting does not apply any weighting of the
%   frequency channels and removes outlier that deviate more than 30 degree from
%   the median.
%
%   WIERSTORF2013_ESTIMATEAZIMUTH accepts the following optional parameters:
%
%     'fs',fs                  Sampling rate
%
%     'dietz2011'              Use the dietz2011 binaural model to estimate the
%                              azimuth value. This is the default.
%
%     'lindemann1986'          Use the lindemann1986 binaural model to estimate
%                              the azimuth value.
%
%     'no_spectral_weighting'  Apply equal weighting of all frequency channels.
%                              This is the default behavior.
%
%     'rms_weighting'          Weight the frequency channels according their rms
%                              value of the signal.
%
%     'raatgever_weighting'    Weight the frequency channels after the empirical
%                              curve from Raatgever that has a maximum around
%                              600 Hz. Note, that this works well only for
%                              special stimuli.
%
%     'remove_outlier'         Remove frequency channels from azimuth
%                              calculation that deviate more than 30 degree from
%                              the median azimuth.
%
%     'include_outlier'        Use all azimuth values to calculate the median.
%                              Note, this can lead to NaN if one of the
%                              frequency channels has a NaN as direction.
%
%   See also: wierstorf2013, dietz2011, lindemann1986, itd2angle_lookuptable
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592--605, 2011. [1]http ]
%     
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608--1622, 1986.
%     
%     J. Raatgever. On the binaural processing of stimuli with different
%     interaural phase relations. PhD thesis, TU Delft, 1980.
%     
%     R. Stern, A. Zeiberg, and C. Trahiotis. Lateralization of complex
%     binaural stimuli: A weighted-image model. J. Acoust. Soc. Am.,
%     84:156--165, 1988.
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/wierstorf2013_estimateazimuth.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA SFS
%   #Author: Hagen Wierstorf

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end

if ~isstruct(lookup)
    error('%s: lookup has to be a struct!',upper(mfilename));
end

definput.keyvals.fs = 44100;
definput.flags.binaural_model = {'dietz2011','lindemann1986'};
definput.flags.spectral_weighting = {'no_spectral_weighting','rms_weighting','raatgever_weighting'};
definput.flags.outlier = {'remove_outlier','include_outlier'};

[flags,kv]  = ltfatarghelper({},definput,varargin);


%% ===== Computation ====================================================
%
% === Calculate azimuth values for every frequency channel ===
if flags.do_dietz2011
    ic_threshold=0.98;
    % Run the Dietz model on signal
    [fine,cfreqs,ild] = dietz2011(insig,kv.fs,'nolowpass','fhigh',1400);
    % Unwrap ITDs and get the azimuth values
    itd = dietz2011_unwrapitd(fine.itd,ild,fine.f_inst,2.5);
    phi = itd2angle(itd,lookup);
    % Calculate the median over time for every frequency channel of the azimuth
    for n = 1:size(phi,2)
        idx = fine.ic(:,n)>ic_threshold&[diff(fine.ic(:,n))>0; 0]; % compare eq. 9 in Dietz (2011)
        angle = phi(idx,n);
        idx = ~isnan(angle);
        if size(angle(idx),1)==0
            azimuth(n) = NaN;
            azimuth_std(n) = NaN;
        else
            azimuth(n) = median(angle(idx));
            azimuth_std(n) = std(angle(idx));
        end
    end
    % Calculate ITD and ILD values
    itd = median(itd,1);
    % weights for rms-weighting
    rms_weights = fine.rms;
elseif flags.do_lindemann1986
    % run Lindemann model on signal
    c_s = 0.3; % stationary inhibition
    w_f = 0; % monaural sensitivity
    M_f = 6; % decrease of monaural sensitivity
    T_int = 6; % integration time
    N_1 = 1764; % sample at which first cross-correlation is calculated
    [cc_tmp,~,ild,cfreqs] = lindemann1986(insig,kv.fs,c_s,w_f,M_f,T_int,N_1);
    cc_tmp = squeeze(cc_tmp);
    % Calculate tau (delay line time) axes
    tau = linspace(-1,1,size(cc_tmp,2));
    % find max in cc
    itd = zeros(size(cc_tmp,1),size(cc_tmp,3));
    for ii=1:size(cc_tmp,1)
        for jj=1:size(cc_tmp,3)
            [~,idx] = max(cc_tmp(ii,:,jj));
            itd(ii,jj) = tau(idx)/1000;
        end
    end
    azimuth = itd2angle(itd,lookup);
    azimuth_std = std(azimuth);
    azimuth = median(azimuth);
    % TODO: rms weights
end

% === Weights for cross-frequency integration ===
if flags.do_no_spectral_weighting
    w = ones(size(azimuth));
end
if flags.do_rms_weighting
    w = rms_weights;
end
if flags.do_raatgever_weighting
    % Calculate a spectral weighting after Stern1988, after the data of
    % Raatgever1980
    b1 = -9.383e-2;
    b2 =  1.126e-4;
    b3 = -3.992e-8;
    w = 10.^( -(b1*f+b2*(f).^2+b3*(f).^3)/10 );
end

if flags.do_remove_outlier
    % Remove outliers
    [azimuth,azimuth_std,itd,cfreqs,w] = remove_outlier(azimuth,azimuth_std,itd,cfreqs,w);
end

% Calculate mean about frequency channels
if length(azimuth)==0
    phi = NaN;
    phi_std = NaN;
else
    if flags.do_no_spectral_weighting
        phi = median(azimuth);
        phi_std = median(azimuth_std);
    else
        phi = sum(azimuth.*w)/sum(w);
        phi_std = sum(azimuth_std.*w)/sum(w);
    end
end

end % of main function

%% ===== Subfunctions ====================================================
function [azimuth,azimuth_std,itd,cfreqs,w] = remove_outlier(azimuth,azimuth_std,itd,cfreqs,w)
    % NOTE: the following was enabled for the original paper
    % remove unvalid ITDs
    %azimuth = azimuth(abs(itd)<0.001);
    %azimuth_std = azimuth_std(abs(itd)<0.001);
    %cfreqs = cfreqs(abs(itd)<0.001);
    %w = w(abs(itd)<0.001);
    %itd = itd(abs(itd)<0.001);
    % remove NaN
    w = w(~isnan(azimuth));
    itd = itd(~isnan(azimuth));
    cfreqs = cfreqs(~isnan(azimuth));
    azimuth = azimuth(~isnan(azimuth));
    azimuth_std = azimuth_std(~isnan(azimuth_std));
    % remove outliers more than 30deg away from median
    if length(azimuth)>0
        itd = itd(abs(azimuth-median(azimuth))<30);
        cfreqs = cfreqs(abs(azimuth-median(azimuth))<30);
        w = w(abs(azimuth-median(azimuth))<30);
        azimuth_std = azimuth_std(abs(azimuth-median(azimuth))<30);
        azimuth = azimuth(abs(azimuth-median(azimuth))<30);
    end
end


