function data = data_glasberg2002(varargin)
%DATA_GLASBERG2002 Filter coefficients of outer and middle ear
%   Usage: 
%     data = data_glasberg2002(varargin)
%     [data] = data_glasberg2002('tfOuterMiddle1997','fieldType','free')
%     [data] = data_glasberg2002('tfOuterMiddle1997','fieldType','diffuse')
%
%   DATA_GLASBERG2002 provides the filter coefficients of the outer and middle ear as determined in Glasberg et al. (2002)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_glasberg2002.php


%   #Author: Thomas Deppisch
%   #Author: Piotr Majdak (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

    definput.flags.type = {'missingflag','tfOuterMiddle1997','tfOuterMiddle2007','specLoud'};
                                       
    definput.keyvals.fieldType='free';
    definput.keyvals.fVec=[];

    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
           flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
           sprintf('%s or %s',definput.flags.type{end-1},...
           definput.flags.type{end})];
           error('%s: You must specify one of the following flags: %s.',...
                 upper(mfilename),flagnames);
    end
    if isempty(kv.fVec)
        error('Please specify frequency vector fVec.')
    end

    %% Returns transfer function of outer and middle ear as used in moore1997, glasberg2002
    if flags.do_tfOuterMiddle1997
        % transfer function of the outer ear
        fOuter = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 750 ...
            800 1000 1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 6300 8000 ...
            9000 10000 11200 12500 14000 15000 16000 20000];

        % values of ANSI S3.4-2007
        if kv.fieldType == 'free' % free field
            tfOuter = [0 0 0 0 0 0 0 0 0.1 0.3 0.5 0.9 1.4 1.6 1.7 2.5 2.7 2.6 2.6 3.2 5.2 ...
                6.6 12 16.8 15.3 15.2 14.2 10.7 7.1 6.4 1.8 -0.9 -1.6 1.9 4.9 2 -2 2.5 2.5];

        elseif kv.fieldType == 'diffuse' % diffuse field
            tfOuter = [0 0 0 0 0 0 0 0 0.1 0.3 0.4 0.5 1 1.6 1.7 2.2 2.7 2.9 3.8 5.3 6.8 7.2 ...
                10.2 14.9 14.5 14.4 12.7 10.8 8.9 8.7 8.5 6.2 5 4.5 4 3.3 2.6 2 2];
        else
            error('Wrong parameter for fieldType. Please use "free" or "diffuse".')
        end

        tfOuterInterp = interp1(fOuter, tfOuter, kv.fVec, 'pchip');

        % transfer function of the middle ear
        fMiddle = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 750 ...
            800 1000 1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 6300 8000 ...
            9000 10000 11200 12500 14000 15000 16000 18000 20000];

        % revised data 2006
        tfMiddle = -[39.6 32 25.85 21.4 18.5 15.9 14.1 12.4 11 9.6 8.3 7.4 6.2 4.8 ...
            3.8 3.3 2.9 2.6 2.6 3.2 4.5 5.5 8.5 10.4 7.3 7 6.6 7 9.2 10.2 12.2 ...
            10.8 10.1 12.7 15 18.2 23.8 32.3 45.5 50];

        tfMiddleInterp = interp1(fMiddle, tfMiddle, kv.fVec, 'pchip');

        data.tfOuterMiddle = tfOuterInterp + tfMiddleInterp;
        data.tfOuter = tfOuter;
        data.tfMiddle = tfMiddle;
        data.fOuter = fOuter;
        data.fMiddle = fMiddle;
    end

    %% Returns revised transfer function of outer and middle ear as in ANSI S3.4-2007
    if flags.do_tfOuterMiddle2007
        % transfer function of the outer ear
        fOuter = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 750 ...
            800 1000 1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 6300 8000 ...
            9000 10000 11200 12500 14000 15000 16000 20000];

        % values of ANSI S3.4-2007
        if kv.fieldType == 'free' % free field
            tfOuter = [0 0 0 0 0 0 0 0 0.1 0.3 0.5 0.9 1.4 1.6 1.7 2.5 2.7 2.6 2.6 3.2 5.2 ...
                6.6 12 16.8 15.3 15.2 14.2 10.7 7.1 6.4 1.8 -0.9 -1.6 1.9 4.9 2 -2 2.5 2.5];

        elseif kv.fieldType == 'diffuse' % diffuse field
            tfOuter = [0 0 0 0 0 0 0 0 0.1 0.3 0.4 0.5 1 1.6 1.7 2.2 2.7 2.9 3.8 5.3 6.8 7.2 ...
                10.2 14.9 14.5 14.4 12.7 10.8 8.9 8.7 8.5 6.2 5 4.5 4 3.3 2.6 2 2];
        else
            error('Wrong parameter for fieldType. Please use "free" or "diffuse".')
        end

        tfOuterInterp = interp1(fOuter, tfOuter, kv.fVec, 'pchip');

        % transfer function of the middle ear
        fMiddle = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 750 ...
            800 1000 1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 6300 8000 ...
            9000 10000 11200 12500 14000 15000 16000 18000 20000];

        % values of ANSI S3.4-2007
        tfMiddle = -[39.6 32 25.85 21.4 18.5 15.9 14.1 12.4 11 9.6 8.3 7.4 6.2 4.8 3.8 3.3 2.9 2.6 2.6 4.5 5.4 6.1 8.5 10.4 7.3 7 ...
                    6.6 7 9.2 10.2 12.2 10.8 10.1 12.7 15 18.2 23.8 32.3 45.5 50];

        tfMiddleInterp = interp1(fMiddle, tfMiddle, kv.fVec, 'pchip');

        data.tfOuterMiddle = tfOuterInterp + tfMiddleInterp;
    end

    %% Return data for calculating specific loudness as in ANSI S3.4-2007
    if flags.do_specLoud
        fRef = [50 63 80 100 125 160 200 250 315 400 500 630 750 800 1000];
        tQ = [28.18 23.9 19.2 15.68 12.67 10.09 8.08 6.3 5.3 4.5 3.73 3.73 3.73 3.73 3.73];
        data.tQ = interp1(fRef, tQ, kv.fVec, 'pchip');
        data.tQ500 = tQ(11);
        data.g = data.tQ500-data.tQ;    % low level gain in cochlea amplifier

        %% linearization parameter a
        g = [2454531 2378397 2278169 1978305 1621055 1123866 945902 738338 589392 497718 362882 250042 177405 ...
                 157745 124006 94596 75663 52501 35451 18750 7845 2470 0] * -1e-5;
        a = [885200 863150 835840 765258 699954 610719 582490 555322 536425 525064 510314 498203 490515 ...
                 488449 484918 481854 479890 477494 475736 474019 472899 472349 472096] * 1e-5;

        data.a = interp1(g, a, data.g, 'pchip');

        %% compressive exponent alpha
        g = [-25 -20 -15 -10 -5 -0];
        alpha = [26692 25016 23679 22228 21055 20000]*1e-5;

        data.alpha = interp1(g, alpha, data.g, 'pchip');

        data.c = 0.046871; % constant to get loudness scale to sone
    end
end


