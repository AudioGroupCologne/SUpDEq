function data = data_lyon2011(varargin)
%DATA_LYON2011   Data for demo_lyon2011CompressiveFunction, demo_lyon2011ExcitationPattern, demo_lyon2011ImpulseResponses
%
%   Usage: data = data_lyon2011(varargin)
%
%   Output parameters:
%     data     : structure
%
%   DATA_LYON2011(varargin) returns data from rhode1996, lopezpoveda2003,
%   russel1997, ren2002, glasberg1990 and deboer2000 to compare with the
%   results of the lyon2011 model.
%
%   The parameter varagin may be one of: rhode1996, lopezpoveda2003,
%   russel1997, ren2002, glasberg1990 or deboer2000.
%   See the references below for the corresponding papers.
%
%   The fields in the output data contain the following information:
%
%     'rhode1996'       physiological chinchila data at CF=500
%
%                       - L_animal_500hz           
%                       - IO_animal_norm_500hz
%
%     'lopezpoveda2003'   psychoacoustic data at CF = 500Hz, 1kHz, 2kHz, 4kHz
%
%                         - L_psych_500hz
%                         - IO_ex_norm_500hz
%                         - L_psych_1khz
%                         - IO_ex_norm_1khz
%                         - L_psych_2khz
%                         - IO_ex_norm_2khz
%                         - L_psych_4khz
%                         - IO_ex_norm_4khz
%
%     'russel1997'        physiological chinchila data
%
%                         - L_animal_4khz
%                         - IO_animal_norm_4khz
%
%     'ren2002'           physiological data (Ren,2002), Fig. 1, panel A.
%
%                         - P_ex
%                         - Ex_30dB_norm
%
%     'glasberg1990'      Experimentally-derived equation of Glasberg and Moore (1990)
%
%                         - f
%                         - QERB_exp
%
%     'deboer2000'        The QERBs derived from the impulse responses recorded by deBoer and Nuttal (2000), Fig. 1.
%
%                         - intensity_deBoer
%                         - QERB_deBoer
%
%   Examples:
%
%   To get provided data from lopezpoveda2003, use :
%
%     data = data_lyon2011('lopezpoveda2003')
%
%
%   References:
%     I. Russel and K. Nilsen. The location of the cochlear amplifier:
%     Spatial representation of a single tone on the guinea pig
%     basilar membrane. Proceedings of the National Academy of Sciences of
%     the United States of America, 94(6):2660--2664, 1997.
%     
%     E. Lopez-Poveda, C. J. Plack, and R. Meddis. Cochlear nonlinearity
%     between 500 and 8000 hz in listeners with normal hearing. J. Acoust.
%     Soc. Am., 113(951), 2003.
%     
%     W. Rhode and N. Cooper. Nonlinear mechanics in the apical turn of the
%     chinchilla cochlea in vivo,. Aud. Neurosci., 3:101--121, 1996.
%     
%     R. F. Lyon. Cascades of two-pole–two-zero asymmetric resonators are
%     good models of peripheral auditory function. J. Acoust. Soc. Am.,
%     130(6), 2011.
%     
%     B. R. Glasberg and B. Moore. Derivation of auditory filter shapes from
%     notched-noise data. Hearing Research, 47(1-2):103--138, 1990.
%     
%
%   See also: demo_lyon2011_compressivefunctions demo_lyon2011 demo_lyon2011_impulseresponses
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_lyon2011.php


%   #Author: Amin Saremi
%   #Author: Clara Hollomey (2021): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


    %% ------ Check input options -----------------------------------------
    % Define input flags
    definput.flags.type = {'missingflag','rhode1996','lopezpoveda2003',...
        'russel1997','ren2002','glasberg1990','deboer2000'};

    % Parse input options
    flags=ltfatarghelper({},definput,lower(varargin));

    if flags.do_missingflag
      flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
                 sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
      error('%s: You must specify one of the following flags: %s.',mfilename,flagnames);
    end
    
    %% ------ Get data from different papers ------------------------------
    if flags.do_rhode1996
        %% Physiological and psychoacoustic data at CF=500 Hz
        % physiological chinchila data from Rhode and Cooper(1996)
        data.L_animal_500hz=[25,31,37,42,48,54,60,66,72,78,84];
        IO_raw_500hz=[1.5,3,5,10,20,30,40,60,90,170,300];
        IO_animal_500hz=db(IO_raw_500hz);
        data.IO_animal_norm_500hz=IO_animal_500hz-IO_animal_500hz(1)+data.L_animal_500hz(1); % Normalize
    end
    
    if flags.do_lopezpoveda2003
        %% Physiological and psychoacoustic data at CF=500 Hz
        data.L_psych_500hz=[15,23,35,45,55];% from Lopez-Poveda et al., (2003)
        IO_psych_500hz=[42,55,63,65,70];
        data.IO_ex_norm_500hz=IO_psych_500hz-IO_psych_500hz(1)+data.L_psych_500hz(1); %Normalize
        %% psychoacoustic data at CF= 1kHz (Lopez-Poveda et al., 2003)
        data.L_psych_1khz=[10,20,30,40,50,60,70,80];
        IO_psych_1khz=[38,50,62,70,82,88,88,92];
        data.IO_ex_norm_1khz=IO_psych_1khz-IO_psych_1khz(1)+data.L_psych_1khz(1); %Normalize
        %% psychoacoustic data at CF= 2kHz (Lopez-Poveda et al., 2003)
        data.L_psych_2khz=10:10:100;% 
        IO_psych_2khz=[35,43,52,61,64,68,72,76,80,83];
        data.IO_ex_norm_2khz=IO_psych_2khz-IO_psych_2khz(1)+data.L_psych_2khz(1); %Normalize
        %% Physiological and psychoacoustic data at CF= 4kHz
        data.L_psych_4khz=10:10:100;% Psychoacoustic data by Lopez-Poveda et al. (2003)
        IO_psych_4khz=[31,38,48,53,60,64,65,67,68,69];
        data.IO_ex_norm_4khz=IO_psych_4khz-IO_psych_4khz(1)+data.L_psych_4khz(1); %Normalize
    end
    
    if flags.do_russel1997
        %% physiological chinchila data by Russel and Nilsen (1997)
        data.L_animal_4khz=10:10:80;
        IO_raw_4khz=[0.25,1,3,4.25,7,9,10,10.5];
        IO_animal_4khz=db(IO_raw_4khz);
        data.IO_animal_norm_4khz=IO_animal_4khz-IO_animal_4khz(1)+data.L_animal_4khz(1); %Normalize
    end
    
    if flags.do_ren2002
        %% physiological data (Ren,2002), Fig. 1, panel A. 
        data.P_ex= 0.01.*[11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16]+0.2;% in percent from base
        Ex_30dB=1e-6.*[7,10,18,30,50,80,50,30,20,10,8];Ex_30dB=20.*log10(Ex_30dB./2e-5);
        data.Ex_30dB_norm=Ex_30dB-max(Ex_30dB);
    end
    
    if flags.do_glasberg1990
        % Experimentally- derived equation of Glasberg and Moore (1990)
        data.f=[500,1000,2000,4000];
        data.QERB_exp=[6.35,7.53,8.31,8.76]; % the QERBs at 0.5, 1, 2 and 4 kHz at low intensities according to Glasberg and Moore (1990).
    end
    
    if flags.do_deboer2000
        % The QERBs derived from the impulse responses recorded by deBoer and
        % Nuttal (2000), Fig. 1.
        data.intensity_deBoer=[20,60,70,80,90,100];
        data.QERB_deBoer=[8.8923    7.7079    5.3953    4.8953    4.5704    2.9763];
    end
    
end


