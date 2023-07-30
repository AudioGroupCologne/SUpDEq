function varargout = exp_tabuchi2016(varargin)
%EXP_TABUCHI2016 Results from Tabuchi et al. (2016)
%   Usage: data = exp_tabuchi2016(flag) 
%
%   EXP_TABUCHI2016(flag) reproduces figures of the study from 
%   Tabuchi et al. (2016).
%
%   The following flags can be specified
%
%     'fig6'    Reproduces the lower panel of Figure 6
%
%   Examples:
%   ---------
%
%   To display Fig.6 use :
%
%     exp_tabuchi2016('fig6');
% 
%   See also: tabuchi2016
%
%   References:
%     H. Tabuchi, B. Laback, T. Necciari, and P. Majdak. The role of
%     compression in the simultaneous masker phase effect. The Journal of the
%     Acoustical Society of America, 140(4), 2016.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_tabuchi2016.php


%   #Author: Hisaaki Tabuchi (2022)
%   #Author: Clara Hollomey (2023): adaptations for AMT


definput.import = {'amt_cache'};
definput.flags.type = {'missingflag','fig6'};

[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

if flags.do_fig6
    
    [kfunc, kvals] = amt_cache('get','kfunction',flags.cachemode);
    
    if isempty(kfunc) || isempty(kvals)

        % initial parameters
        n1 = 34;
        n2 = 46;
        Cvec = -1.5:0.25:1; % This range of Cvec was tested by considering the presumed curvature, -0.5. See the footnote 2 in Tabuchi et al. (2016). 
        GmaxVec = 0:70;
        x_vec = 0:0.1:110; % target level, dB SPL with 0.1 dB step
        mlvl_vec = [60 90];
        GmaxToGetK = 34;
        gamma = 60; % for IHC 
        beta = 1; % for IHC
        sfs = 48;
        durms = 40;
        f0 = 0.1;
        ListenerStr = 'GrandMeanSd';
        CondStr_vec = {'OffFreqPre', 'OnFreqPre'};
        
        KvecMat = NaN(length(x_vec), length(Cvec), length(GmaxVec),...
        length(CondStr_vec), length(mlvl_vec));
        % Kval_Mat is computed by C vals between -1 and 1.
        CvecLen = max(find(Cvec >= -1)) - min(find(Cvec >= -1)) + 1;
        Kval_Mat = NaN(CvecLen, length(CondStr_vec), length(mlvl_vec));
        
        
        for mlvlLoop = 1:length(mlvl_vec)
            for PrecCondLoop = 1:length(CondStr_vec)

                % get the data of grand mean
                dataOut = amt_load('tabuchi2016', [ListenerStr '_' ...
                    CondStr_vec{PrecCondLoop} '_' int2str(mlvl_vec(mlvlLoop)) 'dB' '.mat']);
                S= fieldnames(dataOut);
                dataOut = eval(['dataOut.' S{1}]);
                CvecLen_flag = 0; % to sort out the if sentence below

                for CvecLoop = 1:length(Cvec)  
                    % generate a masker
                    [WaveMout, nPhase] = sig_tabuchi2016(CondStr_vec(PrecCondLoop),...
                        sfs, f0, n1, n2, mlvl_vec(mlvlLoop), durms, Cvec(CvecLoop));
                    
                    for GmaxLoop = 1:length(GmaxVec)
                    [Kfunc_vec, Kval] = tabuchi2016(WaveMout, nPhase, dataOut,'C', Cvec(CvecLoop), ...%'n1', n1, 'n2', n2,
                        'Gmax', GmaxVec(GmaxLoop), 'mlvl', mlvl_vec(mlvlLoop),...
                        'GmaxToGetK', GmaxToGetK, 'gamma', gamma, 'beta', beta);
                    if ~isempty(Kval) && (GmaxVec(GmaxLoop) == GmaxToGetK) && (Cvec(CvecLoop)>= -1)
                        CvecLen_flag = CvecLen_flag + 1; 
                        Kval_Mat(CvecLen_flag, PrecCondLoop, mlvlLoop) = Kval; % for Kave
                    end
                    KvecMat(:, CvecLoop, GmaxLoop, PrecCondLoop, mlvlLoop) = Kfunc_vec;
                    
                    end
                end
            end
        end
        kfunc = KvecMat;
        % for Kave
        kvals = Kval_Mat;
        amt_cache('set','kfunction',kfunc,kvals);

    end

    gammaVal = 60;
    %Get the predicted masked thresholds of C -1.5 and C 0, when the 
    %difference of the RMS Outputs of M and M+T (across target levels)
    %and the average of K over all the conditions is minimized.
    [Kave, ThresEstVec] = tabuchi2016_estimatethreshold(kfunc, kvals);

    %Plot a difference of the predicted thresholds for C -1.5 and C 0 
    %(i.e., MMD) as a function of Gmax.
    plot_tabuchi2016(gammaVal, Kave, ThresEstVec);
end
 
end

