function [out, clean, noisy] = relanoiborra2019(insig_clean, insig_noisy, fs, varargin)
%RELANOIBORRA2019 Modulation filterbank (based on DRNL)
%
%   Usage: out = relanoiborra2019(insig_clean, insig_noisy, fs, varargin)
%          out = relanoiborra2019(insig_clean, insig_noisy, fs, flow, fhigh, varargin)
%          [out, clean, noisy] = relanoiborra2019([..])
%
%   Input parameters:
%     insig_clean  :  clean speech template signal
%     insig_noisy :  noisy speech target signal
%     fs         :  Sampling frequency
%     flow       :  lowest center frequency of auditory filterbank
%     fhigh      :  highest center frequency of auditory filterbank
%     N_org      :  length of original sentence. Will be double the length 
%                   of 'insig_clean' if not provided. 
%     sbj        :  subject profile for drnl definition. default: 'NH' 
%
%   Output parameters:
%     out           : correlation metric structure
%                     The out structure has the following fields:
%
%                     - .dint : correlation values for each modulation band
%
%                     - .dsegments : correlation values from each time window and mod. band
%
%                     - .dfinal : final (averaged) correlation
%
%   This script builds the internal representations of the template and target signals
%   according to the CASP model (see references).
%   The code is based on previous versions of authors: Torsten Dau, Morten
%   Leve Jepsen, Boris Kowalesky and Peter L. Soendergaard.
%   The model has been optimized to work with speech signals, and the
%   preprocesing and variable names follow this principle. The model is
%   also designed to work with broadband signals. In order to avoid undesired
%   onset enhancements in the adaptation loops, the model expects to recive a
%   prepaned signal to initialize them.
%
%
%   References:
%     H. Relaño-Iborra, J. Zaar, and T. Dau. A speech-based computational
%     auditory signal processing and perception model. J. Acoust. Soc. Am.,
%     146(5), 2019.
%     
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am., 124(1),
%     2008.
%     
%
%   See also: ihcenvelope relanoiborra2019_drnl
%             relanoiborra2019_mfbtd joergensen2013_sim
%             exp_osses2022 dau1997
%
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/relanoiborra2019.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Stats M-Signal M-Control
%   #Author: Helia Relano Iborra (March 2019): v4.0 provided to the AMT team
%   #Author: Clara Hollomey (2021): adapted to the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Auditory filtering:
if isoctave
   warning(['Currently this model is only fully functional under MATLAB.']); 
end

definput.import={'relanoiborra2019'}; % load defaults from arg_relanoiborra2019
[flags,kv]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

if isempty(kv.N_org), N_org=length(insig_clean); else N_org=kv.N_org; end

[clean_mfb, fc_mod, clean_afb, fc] = relanoiborra2019_featureextraction(insig_clean, fs, 'argimport',flags,kv);
[noisy_mfb, ~, noisy_afb] = relanoiborra2019_featureextraction(insig_noisy, fs, 'argimport',flags,kv);

out = relanoiborra2019_decision(clean_mfb((N_org+1):end, :, :), noisy_mfb((N_org+1):end, :, :), fs, fc, fc_mod,'argimport',flags,kv);
clean.afb = clean_afb;
clean.fc = fc;
clean.mfb = clean_mfb;
clean.fmod = fc_mod;
noisy.afb = noisy_afb;
noisy.fc = fc;
noisy.mfb = noisy_mfb;
noisy.fmod = fc_mod;

%   [out,fc,mfc] = relanoiborra2019_preproc(insig, fs, varargin);
%     varargout{1} = fc;
%     varargout{2} = mfc;


