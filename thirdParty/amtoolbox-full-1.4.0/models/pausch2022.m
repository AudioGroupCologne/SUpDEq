function [itd,itd_max,itd_arg_max_phi,coef,fig] = pausch2022(featVec,varargin)
%PAUSCH2022 ITD prediction in the horizontal plane for listeners with hearing aids
%
%   Usage: [itd, itd_max, itd_arg_max_phi, coef, fig] = pausch2022(featVec, varargin)
%
%
%   Input parameters:
%     featVec  : Vector containing magnitudes in mm of the features used 
%                for the specific model implementation, see Pausch et al. (2022), Tab. 2 
%                and Sec. 3.4.
%
%                - if 'type_stf'=='hrtf': [x1, x2, x3, x4, x5]
%                - if 'type_stf'=={'hartf_front','hartf_rear'}: [x1, x2, x3, d13, d14]
%
%                - x1: head width in mm [double]
%                - x2: head height in mm [double]
%                - x3: head depth in mm [double]
%                - x4: pinna offset down in mm (only required if type_stf=='hrtf') [double]
%                - x5: pinna offset back in mm (only required if type_stf=='hrtf') [double]
%                - d13: HA-microphone offset up in mm (only required if type_stf=={'hartf_front','hartf_rear'}) [double]
%                - d14: HA-microphone offset back in mm (only required if type_stf=={'hartf_front','hartf_rear'}) [double]                           
%
%   Output parameters:
%     itd             : predicted ITDs in s for type_stf as per type_mod [double]
%     itd_max         : interaural-time-difference maximum, max{ITD} [double]
%     itd_arg_max_phi : argument of the interaural-time-difference maximum,
%                       arg max_phi{ITD} [double]
%     coef            : model coefficients as estimated after applying the 
%                       polynomial regression weights to the subset of
%                       individual features [double]
%     fig             : figure handle [matlab.ui.Figure]
%
%
%   PAUSCH2022() contains a novel hybrid model to predict the interaural time 
%   differences (ITDs) in the horizontal plane for adults with normal hearing, 
%   listening via head-related transfer functions (HRTFs), or adults fitted 
%   with behind-the-ear hearing aids (HA), listening via hearing-aid-related 
%   transfer functions (HARTFs). It also contains two previous analythical 
%   ITD models: Kuhn (1977), Woodworth and Schlosberg (1954), and  Woodworth and 
%   Schlosberg (1954) extended by Aaronson and Hartmann (2014). The ITD predictions 
%   in all models are based on invividual anthropometric features, or features 
%   describing the individual placement of the HAs, see Pausch et al. (2022) 
%   for further details.
%
%   The type_mod flag may be used to select the ITD model for the ITD 
%   predictions (only required if type_stf=='weights'):
%
%     'pausch'         hybrid ITD model by Pausch et al. (2022) (default)
%     'kuhn'           analytic ITD model by Kuhn (1977)
%     'woodworth'      analytic ITD model by Woodworth and Schlosberg (1954)
%     'woodworth_ext'  analytic ITD model by Woodworth and Schlosberg (1954) extended 
%                      by Aaronson and Hartmann (2014) (far-field assumption)
%
%
%   Note: To evaluate type_mod=='woodworth_ext', 
%   additional key/value pairs for the features 
%   x5 (only if type_stf=={'hartf_front','hartf_rear'}), 
%   d9, d10, and Theta3 have to be specified.
%
%   The type_stf flag may be used to choose between one of the following:
%
%     'hrtf'           load individual HRTF datasets
%     'hartf_front'    load invidual front HARTF datasets
%     'hartf_rear'     load invidual rear HARTF datasets
%
%   The plot flags may be:
%
%     'no_plot'        No plot (default).
%     'plot'           plot predicted ITDs over azi_min:azi_res:azi_max 
%                      (default: false) [logical]
%
%   Additional key/value pairs include:
%
%     'd9'             HA-microphones-to-ear-canal offset in mm 
%                      (default: []) [double]
%     'd10'            HA-microphones-to-scalp offset in mm
%                      (default: []) [double]
%     'Theta3'         frontal HA-microphones angle in deg
%                      (default: []) [double]
%     'azi_min'        minimum evaluation angle of azimuth range in deg 
%                      (default: 0) [double]
%     'azi_max'        maximum evaluation angle of azimuth range in deg 
%                      (default: 180) [double]
%     'azi_res'        angular resolution of the evaluated azimuth range in deg 
%                      (default: 2.5) [double]
%     'c'              speed of sound in m/s (default: 343) [double]
%
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/pausch2022.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Author: Florian Pausch (2022): integration in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   [2] G. F. Kuhn, "Model for the interaural time differences in the azimuthal 
%       plane," The Journal of the Acoustical Society of America, vol. 62, 
%       no. 1, pp. 157–167, 1977. doi: 10.1121/1.381498.
%   [3] R. S. Woodworth and H. Schlosberg, Experimental psychology, Rev. ed. 
%       Oxford, England: Holt, 1954.
%   [4] N. L. Aaronson and W. M. Hartmann, "Testing, correcting, and extending 
%       the Woodworth model for interaural time difference," The Journal of 
%       the Acoustical Society of America, vol. 135, no. 2, pp. 817–823, 2014. 
%       doi: 10.1121/1.4861243.

%% Parse flags and keyvals

definput.import={'pausch2022'};
[flags,kv] = ltfatarghelper({'azi_min','azi_max','azi_res'}, definput, varargin);

if strcmp(flags.type_mod,'woodworth_ext')...
        && ( isempty(kv.d9) || isempty(kv.d10) || isempty(kv.Theta3) )
    error('Please specify key/value pairs for ''d9'', ''d10'', and ''Theta3''.')
end

if (strcmp(flags.type_mod,'woodworth_ext') && ~strcmp(flags.type_stf,'hrtf'))...
        && ( isempty(kv.d9) || isempty(kv.d10) || isempty(kv.Theta3) || isempty(kv.x5) )
    error('Please specify key/value pair for ''x5''.')
end

%% Load the polynomials regression weights

weights = data_pausch2022('weights',flags.type_mod,flags.type_stf);
fns = fieldnames(weights);
weights = weights.(fns{1});

%% Estimate the model coefficients by applying the polynomial regression weights
%  on the feature subset magnitudes (featVec)

% construct vector/matrix alpha
M = length(featVec); % number of anthropometric features
P = (size(weights,2)-1)/M; % polynomial degree

temp_alpha = NaN(P,M); 
for mdx=1:M
    temp_alpha(:,mdx)=featVec(mdx).^(1:P);
end
alpha_mtx = [1, temp_alpha(:)'];

% estimate the model coefficients
coef = alpha_mtx*weights'; % effective head radius a_n = coef(1) in m

%% Evaluate the selected model

phi_vec = kv.azi_min:kv.azi_res:kv.azi_max;
phi_vec_rad = phi_vec*pi/180;

switch flags.type_mod

    case 'kuhn'

        itd = 3*coef/kv.c*sin(phi_vec_rad);
    
    case {'woodworth','woodworth_ext'}
    
        if strcmp(flags.type_stf,'hrtf')
            Theta_E = 90 + atand( featVec(5) / sqrt((coef*1e3)^2-featVec(5)^2) );
            if strcmp(flags.type_mod,'woodworth')
                itd = eval_woodworth_ext(coef, kv.c, phi_vec);
            else
                itd = eval_woodworth_ext(coef, kv.c, phi_vec, Theta_E);
            end
        else
            Theta_E = 90 + atand( kv.x5 / sqrt((coef*1e3)^2-kv.x5^2) );
            zeta = kv.d9*sind(kv.Theta3);
            Theta_HA = acosd( ((coef*1e3)^2 + (coef*1e3+kv.d10)^2 - zeta^2) / ...
                (2*coef*1e3*(coef*1e3+kv.d10)) ) + Theta_E;
            if strcmp(flags.type_mod,'woodworth')
                itd = eval_woodworth_ext(coef, kv.c, phi_vec);
            else
                itd = eval_woodworth_ext(coef, kv.c, phi_vec, Theta_HA);
            end
        end
    
    case 'pausch'
        
        itd = coef(1)/kv.c*( coef(2).*sin(phi_vec_rad) + ...
            coef(3).*sin(2*phi_vec_rad) + coef(4).*sin(3*phi_vec_rad) );

end

[itd_max,itd_arg_max_phi] = max(itd);

%% Optionally plot the predicted ITDs

if strcmp(flags.plot,'plot')
    
    switch flags.type_mod
        case 'kuhn'
            col = [.4 .4 .4];
        case {'woodworth','woodworth_ext'}
            col = [0, 0, 0];
        case 'pausch'
            col = [0,84,159]/255;
    end
    
    lwidth = 1.5;
    fsize = 12;
    xpos_anno = 0.97;
    ypos_anno = 0.94;

    fig = figure;
    plot(phi_vec,itd.*1e6,'color',col,'linewidth',lwidth)
    grid on
    title([flags.type_stf,', ',flags.type_mod],'interpreter','latex')
    set(gca,'xlim',[phi_vec(1) phi_vec(end)],'ticklabelinterpreter','latex','fontsize',fsize)
    set(gca,'xtick',phi_vec(1):30:phi_vec(end),'ytick',0:100:max(itd)*1e6)
    axis square
    box on
    xlabel('Azimuth (deg)','fontsize',fsize,'interpreter','latex')
    ylabel('ITD ($\mu$s)','fontsize',fsize,'interpreter','latex')

    if strcmp(flags.type_mod,'woodworth_ext')
        if strcmp(flags.type_stf,'hrtf')
            text(xpos_anno,ypos_anno,['$\Theta_{\mathrm{E}}$\,=\,',...
                num2str(round(Theta_E*10)/10),'$^\circ$'],'interpreter','latex',...
                'fontsize',fsize,'units','normalized','horizontalalignment','right')
        else
            text(xpos_anno,ypos_anno,['$\Theta_{\mathrm{HA}}$\,=\,',...
                num2str(round(Theta_HA*10)/10),'$^\circ$'],'interpreter','latex',...
                'fontsize',fsize,'units','normalized','horizontalalignment','right')
        end
    end

end

end

%% ------------------------------------------------------------------------
%  ---- INTERNAL FUNCTIONS ------------------------------------------------
%  ------------------------------------------------------------------------

function itd = eval_woodworth_ext(a, c, phi_vec, Theta_E)
%EVAL_WOODWORTH_EXT - Function to estimate interaural time differences (ITDs) 
%                     for directions in the horizontal plane based on the 
%                     Woodworth model [1] extended by specified horizontal 
%                     ear-canal angles Theta_E (or horizontal HA-microphones 
%                     angle Theta_HA), assuming an infinite source distance [2]. 
%                     The model is evaluated for source directions phi_vec(phi_vec<=pi) 
%                     in the horizontal plane.
%
%   Usage : itd = eval_woodworth_ext(a, c, phi_vec) [1]
%           itd = eval_woodworth_ext(a, c, phi_vec, Theta_E) [2]
% 
%   Input parameters (required):
%
%     a       : effective head radius (m) [double]
%     c       : speed of sound (m/s) [double]
%     phi_vec : vector of azimuth angles <= 180 (deg) [double]
%     Theta_E : horizontal ear-canal angle (deg) [double]
%
%   Output parameters:
%
%     itd : if no Theta_E is specified: predicted ITDs are based on the simple 
%           Woodworth model [1] (s) [double]
%           if Theta_E is specified: predicted ITDs are based on the extended 
%           Woodworth model [2] (s) [double]
%
%   [1] R. S. Woodworth, "Experimental Psychology," The Journal of Nervous and 
%       Mental Disease, vol. 91, no. 6, p. 811, 1940.
%   [2] N. L. Aaronson and W. M. Hartmann, "Testing, correcting, and extending 
%       the Woodworth model for interaural time difference," The Journal of 
%       the Acoustical Society of America, vol. 135, no. 2, pp. 817–823, 2014. 
%       doi: 10.1121/1.4861243.

%   Author: Florian Pausch, Institute for Hearing Technology and
%            Acoustics, RWTH Aachen University

if nargin < 4
    Theta_E = 90;
end

if any(phi_vec>180)
    error('Azimuth angles exceeding $180^\circ$ in ''phi_vec'' are not supported by this model. Please re-specify.')
end

phi_vec1 = phi_vec(phi_vec<=90);
phi_vec1 = phi_vec1(:);
phi_vec2 = phi_vec(phi_vec>90);
phi_vec2 = phi_vec2(:);

phi_vec1_rad = deg2rad(phi_vec1);
phi_vec2_rad = deg2rad(phi_vec2);

if nargin < 4 || Theta_E==90 % simple Woodworth model

    ITD_Woodworth1 = a/c * (sin(phi_vec1_rad) + phi_vec1_rad);
    ITD_Woodworth2 = a/c * (pi - phi_vec2_rad + sin(phi_vec2_rad));
    itd = [ITD_Woodworth1; ITD_Woodworth2];

else % extended Woodworth model
    
    ITD_Woodworth_ext1 = NaN(numel(phi_vec1),1);
    ITD_Woodworth_ext2 = NaN(numel(phi_vec2),1);

    Theta_E_input = Theta_E;

    if Theta_E_input ~= 90

        if Theta_E_input<90
            Theta_E = 90 + (90-Theta_E_input);
        end

        % regions for Eq. (b1)-(b5)
        line1 = Theta_E - 90;
        line2 = 180 - Theta_E;
        line3 = 270 - Theta_E;

        region_b1 = phi_vec1<line1 & phi_vec1<line2;
        region_b2 = phi_vec1>=line1 & phi_vec1<line2;
        region_b3a = phi_vec1>=line2 & phi_vec1>=line1;
        region_b3b = phi_vec2<line3;
        region_b4 = phi_vec2>=line3;
        region_b5 = phi_vec1>=line2 & phi_vec1<line1;

        Theta_E_rad = deg2rad(Theta_E);

        % evaluate ITD for the different regions depending on the ear angles
        ITD_Woodworth_ext1(region_b1)  = 2*a/c*phi_vec1_rad(region_b1);
        ITD_Woodworth_ext1(region_b2)  = a/c*(-pi/2 + phi_vec1_rad(region_b2) + Theta_E_rad + cos(phi_vec1_rad(region_b2)-Theta_E_rad));
        ITD_Woodworth_ext1(region_b3a) = a/c*(3*pi/2 - phi_vec1_rad(region_b3a) - Theta_E_rad + cos(phi_vec1_rad(region_b3a)-Theta_E_rad));
        ITD_Woodworth_ext2(region_b3b) = a/c*(3*pi/2 - phi_vec2_rad(region_b3b) - Theta_E_rad + cos(phi_vec2_rad(region_b3b)-Theta_E_rad));
        ITD_Woodworth_ext2(region_b4)  = a/c*(cos(phi_vec2_rad(region_b4)-Theta_E_rad) - cos(phi_vec2_rad(region_b4)+Theta_E_rad));
        ITD_Woodworth_ext1(region_b5)  = 2*a/c*(pi - Theta_E_rad);

        if Theta_E_input>90
            itd = [ITD_Woodworth_ext1; ITD_Woodworth_ext2];
        else % Theta_E_input<90
            itd = flipud([ITD_Woodworth_ext1; ITD_Woodworth_ext2]);
        end
   
    end

end

end


