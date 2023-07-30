function definput=arg_auditorybrainstem(definput)
% ARG_AUDITORYBRAINSTEM
%
% Gets the default values for the auditory-brainstem processing used in
% the models by Verhulst et al. (2015) and (2018). 
%       
%   References:
%     S. Verhulst, H. Bharadwaj, G. Mehraei, C. Shera, and
%     B. Shinn-Cunningham. Functional modeling of the human auditory
%     brainstem response to broadband stimulation. jasa, 138(3):1637--1659,
%     2015.
%     
%
%   #Author: Alejandro Osses (2020)
%   #Author: Piotr Majdak (2021): clean up towards removing this file
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_auditorybrainstem.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

 
%definput.keyvals.model_version = 1.20;
%definput.keyvals.version_year = 2018;
%%% Parameters for Verhulst2018, model version 1.20:
%definput.keyvals.M1 = 4.2767e-14;
%definput.keyvals.M3 = 5.1435e-14;
%definput.keyvals.M5 = 13.3093e-14;

definput.keyvals.Tex_cn = 0.5e-3; % Tau excitation
definput.keyvals.Tin_cn = 2e-3; % Tau inhibition
definput.keyvals.dly_cn = 1e-3; % dly ms
definput.keyvals.Acn    = 1.5;
definput.keyvals.Scn    = 0.6;

definput.keyvals.Tex_ic = 0.5e-3; % Tau excitation
definput.keyvals.Tin_ic = 2e-3; % Tau inhibition
definput.keyvals.dly_ic = 2e-3; % dly ms
definput.keyvals.Aic    = 1.0;
definput.keyvals.Sic    = 1.5;
definput.keyvals.subfs  = [];
definput.keyvals.BMF    = []; % to be used in zilany2014, bruce2018
definput.flags.ic_hwr = {'ic_hwr','no_ic_hwr'}; % to be used in zilany2014, bruce2018
%%% End parameters

definput.flags.ihc  = {'ihc','no_ihc'}; % placed here for convenience
definput.flags.an   = {'an','no_an'}; 
definput.flags.cn   = {'cn','no_cn'}; 
definput.flags.ic   = {'ic','no_ic'}; 
definput.flags.mfb  = {'no_mfb','mfb'};
definput.flags.oae    = {'oae','no_oae'}; 
definput.flags.disp = {'disp','no_disp'}; 

definput.groups.bm_verhulst2012  = {'no_ihc','no_an','no_cn','no_ic','oae','disp', ...
                                    'model_version',1,'version_year',2012};
% definput.groups.abr_verhulst2015 = {'an','cn','ic','no_oae','no_disp','no_ihc', ...
%                                     'subfs',20000,'no_mfb','model_version',2015,'version_year',2015, ...
%                                     'M1',1.4792e-14,'M3',3.5456e-14,'M5',4.53e-14}; % Recalibrated for AMT toolbox
                                    % 'M1',1.845e14,'M3',93.8e-6,'M5',90.6e-6}; % as in Verhulst2015
%definput.groups.abr_verhulst2018 = {'an','cn','ic','noe','no_disp','no_ihc', ...
%                                    'subfs',20000,'no_mfb','model_version',1.2,'version_year',2018}; % default
%definput.groups.ic_verhulst2018  = {'no_ihc','no_an','no_cn','ic','no_e','no_disp', ...
%                                    'subfs',20000,'no_mfb','model_version',1.2,'version_year',2018};
%definput.groups.store_all = {'an','cn','ic','e','disp','ihc','no_mfb', ...
%                             'model_version',1.2,'version_year',2018}; % default

%%% Modulation filter bank (3 filters only), as in Osses and Verhulst 2020:
definput.groups.abr_osses2020 = {'an','cn','ic','no_oae','no_disp','no_ihc', ...
        'subfs',20000,'mfb','model_version',1.299,'version_year',2018, ...
        'M1',4.2702e-14,'M3',1.8092e-14,'M5',1.0876e-13, ... % Calibration on 01/08/2020 using g20191111_testing_2_MFB_calibration.m and three CN/IC configurations.
        'Tex_cn',[0.5  0.5  0.30]*1e-3,'Tin_cn',[2.0  2.0  0.15]*1e-3, ... % Tau excitation and inhibition CN
        'dly_cn',[1.0  1.0  1.20]*1e-3,'Acn',[1.5  1.5  1.5],'Scn',[0.6  0.6  0.6], ...
        'Tex_ic',[ 5   0.5   0.5]*1e-3,'Tin_ic',[10   2.0  0.20]*1e-3, ... % Tau excitation, inhibition IC
        'dly_ic',[2.0  2.0   2.0]*1e-3,'Aic',[1.0  1.0  0.7],'Sic',[1.5  1.5  1.0]};
    
definput.groups.ic_osses2020  = {'an','no_cn','ic','no_oae','no_disp','no_ihc', ...
        'subfs',20000,'mfb','model_version',1.299,'version_year',2018, ... % default
        'M1',4.2702e-14,'M3',1.8092e-14,'M5',1.0876e-13, ... % Calibration on 01/08/2020 using g20191111_testing_2_MFB_calibration.m and three CN/IC configurations.
        'Tex_cn',[0.5  0.5  0.30]*1e-3,'Tin_cn',[2.0  2.0  0.15]*1e-3, ... % Tau excitation and inhibition CN
        'dly_cn',[1.0  1.0  1.20]*1e-3,'Acn',[1.5  1.5  1.5],'Scn',[0.6  0.6  0.6], ...
        'Tex_ic',[ 5   0.5   0.5]*1e-3,'Tin_ic',[10   2.0  0.20]*1e-3, ... % Tau excitation, inhibition IC
        'dly_ic',[2.0  2.0   2.0]*1e-3,'Aic',[1.0  1.0  0.7],'Sic',[1.5  1.5  1.0]};
        % mf     =[22    82     301]; % Empirical centre frequencies of the modulation filter bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


