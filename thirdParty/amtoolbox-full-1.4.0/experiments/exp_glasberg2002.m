function exp_glasberg2002(varargin)
%EXP_GLASBERG2002 Figures several papers by Glasberg et al. (2002)
%
%   Usage:
%     exp_glasberg2002(flags)
%
%   The following flags can be specified;
%
%     'fig1'    Reproduce Fig. 1 of glasberg2002.
%
%     'fig1b'   Similar to fig1 but using 2006 revised data for middle ear
%               filter.
%
%     'fig5'    Reproduce Fig. 5 of glasberg2002.
%
%   Examples:
%   ---------
%   To display Fig.1 of Glasberg et al. (2002) use :
%
%     exp_glasberg2002('fig1');
%
%   To display Fig.1b of Glasberg et al. (2002) use :
%
%     exp_glasberg2002('fig1b');
%
%   To display Fig.5 of Glasberg et al. (2002) use :
%
%     exp_glasberg2002('fig5');
%
%   References:
%     B. R. Glasberg and B. C. J. Moore. A Model of Loudness Applicable to
%     Time-Varying Sounds. J. Audio Eng. Soc, 50(5):331--342, 2002.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_glasberg2002.php


%   #Author: Thomas Deppisch
%   #Author: Piotr Majdak (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Retrieve and compute model paramters
    % Set flags

    definput.flags.type = {'missingflag','fig1','fig1b','fig5'};

    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
           flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
           sprintf('%s or %s',definput.flags.type{end-1},...
           definput.flags.type{end})];
           error('%s: You must specify one of the following flags: %s.',...
                 upper(mfilename),flagnames);
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig5
        fs = 32000;
        t = linspace(0,1,fs);
        inSig = sin(2*pi*1000*t).';
        results = glasberg2002(inSig,fs);
        plot(results.eLdB(500,:))
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig1
        fs = 32000;
        fVec = 20:fs/2;
        data = data_glasberg2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);
        figure
        semilogx(fVec, data.tfOuterMiddle)
        grid on
        xlim([20,16000])
        xlabel('Frequency (Hz)')
        ylabel('Relative Transmission (dB)')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 1b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig1b
        fs = 32000;
        fVec = 20:fs/2;
        data = data_glasberg2002('tfOuterMiddle2007','fieldType','free','fVec',fVec);
        figure
        semilogx(fVec, data.tfOuterMiddle)
        grid on
        xlim([20,16000])
        xlabel('Frequency (Hz)')
        ylabel('Relative Transmission (dB)')
    end


end


