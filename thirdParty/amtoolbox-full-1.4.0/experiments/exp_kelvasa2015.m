function [dataOut] = exp_kelvasa2015(varargin)
%EXP_KELVASA2015 Figures from Kelvasa and Dietz (2015)
%
%   EXP_KELVASA2015(fig) computes data to produce corresponding figure
%   number from the Kelvasa and Dietz 2015 paper.
%
%   Usage:
%     [dataOut] = exp_kelvasa2015('fig8a')
%     [dataOut] = exp_kelvasa2015('fig8a','identifier','BTE','HRTFchannels',[3,4]);
%     [dataOut] = exp_kelvasa2015('fig8a',varargin);
%
%   The following flags can be specified;
%
%     'redo'    Recomputes data for specified figure
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'no_plot'  Don't plot, only return data.
%
%     'fig5'    Reproduce Fig. 5.
%
%     'fig6'    Reproduce Fig. 6.
%
%     'fig8a'    Reproduce Fig. 8a.
%
%     'fig8b'    Reproduce Fig. 8b.
%
%     'fig9a'    Reproduce Fig. 9a.
%
%     'fig10'    Reproduce Fig. 10.
%
%     'fig12'    Reproduce Fig. 12.
%
%   Examples:
%   ---------
%
%   To display Fig. 5 use :
%
%     exp_kelvasa2015('fig5');
%
%   To display Fig. 6 use :
%
%     exp_kelvasa2015('fig6');
%
%   To display Fig. 8a use :
%
%     exp_kelvasa2015('fig8a');
%
%   To display Fig. 8b use :
%
%     exp_kelvasa2015('fig8b');
%
%   To display Fig. 9a use :
%
%     exp_kelvasa2015('fig9a');
%
%   To display Fig. 10 use :
%
%     exp_kelvasa2015('fig10');
%
%   To display Fig. 12 use :
%
%     exp_kelvasa2015('fig12');
%
%   References:
%     D. Kelvasa and M. Dietz. Auditory model-based sound direction
%     estimation with bilateral cochlear implants. Trends in Hearing,
%     19:2331216515616378, 2015.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_kelvasa2015.php


%   #Author: Daryl Kelvasa (2016)
%   #Author: Mathias Dietz (2016)
%   #Author: Thomas Deppisch (2017)
%   #Author: Piotr Majdak (2017)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Retrieve and compute model paramters
    % Set flags

      definput.flags.type = {'missingflag','fig5','fig6','fig7','fig8a',...
                                           'fig8b','fig9a','fig10','fig12'};
      definput.flags.disp = {'no_debug','debug'};
      definput.flags.plot = {'plot','no_plot'};
      definput.flags.plot_stage_fig = {'no_plot_stage_fig','plot_stage_fig'};

    % import default arguments from other functions
    definput.import={'kelvasa2015','amt_cache'}; % defaults from arg_kelvasa2015

    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
           flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
           sprintf('%s or %s',definput.flags.type{end-1},...
           definput.flags.type{end})];
           error('%s: You must specify one of the following flags: %s.',...
                 upper(mfilename),flagnames);
    end


%% Load HRTF data
HRTF = amt_load('kelvasa2015',kv.HRTFfile);
[~,ind_elev] = min(abs(HRTF.SourcePosition(:,2)-kv.HRTFelevation));
[~,ind_dist] = min(abs(HRTF.SourcePosition(:,3)-kv.HRTFsourceDistance));
ind = find(sum([HRTF.SourcePosition(:,2) == HRTF.SourcePosition(ind_elev,2),...
                HRTF.SourcePosition(:,3) == HRTF.SourcePosition(ind_dist,3)],2)...
                ==2);
HRTFnew.SourcePosition = HRTF.SourcePosition(ind,:);
HRTFnew.Data.IR = HRTF.Data.IR(ind,kv.HRTFchannels,:);
HRTFnew.Data.SamplingRate = HRTF.Data.SamplingRate;
HRTF = HRTFnew;

%% Set dB SPL offset
dboffset=71.778;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig5

    %Needed parameters
    savename = ['Figure_5_data_',kv.identifier];

    %Check for preprocessed calibration data
    [dataOut] = amt_cache('get', savename, flags.cachemode);
    if isempty(dataOut)

    signals = {'SSN.wav','FrozenSpeech.wav','PinkNoise.wav','WhiteNoise.wav'};
    levels = 35:70;
    elecPos = [3 4; 7 8; 11 12; 15 16; 19 20];

    %Main loop over all stimuli
    n = 1; numLoops = numel(levels)*numel(signals)*2;
    for ind = 1 : 8

        [signal, fs] = amt_load('kelvasa2015',signals{ceil(ind./2)});
        signal = signal(1:6*fs,:);
        signal = resample(signal,kv.FS_ACE,fs);

        if  ismember(ind,2:4:8)
            %HRTF filter signal and choose microphone channels
            [signal] = HRTFfilter(signal,0,kv,HRTF);
        end

            spikeRatePerBinPerLevel = zeros(kv.numBin,numel(levels));
            currentPerBinPerLevel = zeros(kv.numBin,numel(levels));

        for level = 1 : numel(levels)
            tic

            temp = signal(:,1)./rms(signal(:,1));
            scalor = scaletodbspl(levels(level),[],dboffset);
            scalor = rms(temp.*scalor)/rms(signal(:,1));

            signal = signal.*scalor;

            sigLengthSec = (size(signal,1)/fs);

            [electrodogram, vTime] = ...
                                kelvasa2015_ciprocessing(signal,...
                                kv.FS_ACE,'argimport',flags,kv);

            [APvec] = ...
                                kelvasa2015_anprocessing(electrodogram,...
                                vTime,'argimport',flags,kv);

           [~,spikeRatePerBinPerLevelT] = ...
                                kelvasa2015_anbinning(APvec,...
                                sigLengthSec,'argimport',flags,kv);


           spikeRatePerBinPerLevel(:,level) = ...
                            mean(spikeRatePerBinPerLevelT,2);
           currentPerBinPerLevel(:,level) =  ...
                            (sum(electrodogram(elecPos(:,2),:),2) + ...
                             sum(electrodogram(elecPos(:,1),:),2)).* ...
                                25e-6.*1e-6.*1e9;

            a = toc; timeLeft = round((a*(numLoops - n))/60);
            amt_disp(['Computing Figure 5. ' num2str(n) ' / ' num2str(numLoops) '. Time left (min):' num2str(timeLeft)],'volatile');
            n = n+1;
        end

            dataOut(ind).spikeRatePerBinPerLevel = spikeRatePerBinPerLevel;
            dataOut(ind).currentPerBinPerLevel = currentPerBinPerLevel;

    end

    amt_disp();
    amt_cache('set',savename,dataOut);
    end
    if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig6

        %Needed parameters
        azis = 0:5:90;
        levels = [45,55,65];
        savename = ['Figure_6_data_',kv.identifier];

        %Check for preprocessed calibration data
        [dataOut] = amt_cache('get', savename, flags.cachemode);
        if isempty(dataOut)

        signalName = 'SSN.wav'; 
        [signal, fs] = amt_load('kelvasa2015',signalName);
        signal = signal(1:6*fs,:);
        signal = resample(signal,kv.FS_ACE,fs);

        SpkDiffPerBin = zeros(numel(azis),kv.numBin);
        SpkSumPerBin = zeros(numel(azis),kv.numBin);

        %Main loop over all levels
        n = 1; numLoops = numel(levels)*numel(azis);

        for level = 1 : numel(levels)
            for ang = 1 : numel(azis)
                tic

                %HRTF filter signal and choose microphone channels
                [HRIR] = HRTFfilter(signal,azis(ang),kv,HRTF);

                if azis(ang) == 0
                      temp = HRIR(:,1)./rms(HRIR(:,1));
                      scalor = scaletodbspl(levels(level),[],dboffset);
                      scalor = rms(temp.*scalor)/rms(HRIR(:,1));
                end

                HRIR = HRIR .* scalor;
                sigLengthSec = (size(HRIR,1)/fs);

                spikeRatePerBin = zeros(kv.numBin,2);

                for chan = 1 : 2

                     singChanSig = HRIR(:,chan);

                     [electrodogram, vTime] = ...
                                    kelvasa2015_ciprocessing(singChanSig,...
                                    kv.FS_ACE,'argimport',flags,kv);

                      [APvec] = ...
                                    kelvasa2015_anprocessing(electrodogram,...
                                    vTime,'argimport',flags,kv);

                      [~,spikeRatePerBinT] = ...
                                    kelvasa2015_anbinning(APvec,...
                                    sigLengthSec,'argimport',flags,kv);

                      spikeRatePerBin(:,chan) =  mean(spikeRatePerBinT,2);
                end

                SpkDiffPerBin(ang,:) = spikeRatePerBin(:,2)- spikeRatePerBin(:,1);
                SpkSumPerBin(ang,:) = spikeRatePerBin(:,2) + spikeRatePerBin(:,1);

                a = toc; timeLeft = round((a*(numLoops - n))/60);
                amt_disp(['Computing Figure 6. Time left (min):',num2str(timeLeft)],'volatile');
                n = n+1;
            end
            dataOut(level).signalName = signalName;
            dataOut(level).levelDB = levels(level);
            dataOut(level).SpkDiffPerBin = SpkDiffPerBin;
            dataOut(level).SpkSumPerBin = SpkSumPerBin;
        end
         amt_disp();
         amt_cache('set',savename,dataOut);
        end
        dataOut(1).azis=azis;
        if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig7

        %Needed parameters
        azis = 0:5:90;
        level = 65;
        savename = ['Figure_7_data_',kv.identifier];

        %Check for preprocessed calibration data
        [dataOut] = amt_cache('get', savename, flags.cachemode);
        if isempty(dataOut)

        signalName = {'SSN.wav','PinkNoise.wav'};

        %Main loop over all angles and signal
        n = 1; numLoops = 2*numel(azis);
        for sig = 1 : 2

        [signal, fs] = amt_load('kelvasa2015',signalName{sig});
        signal = signal(1:6*fs,:);
        signal = resample(signal,kv.FS_ACE,fs);

        SpkDiffPerBin = zeros(numel(azis),kv.numBin);
        SpkSumPerBin = zeros(numel(azis),kv.numBin);

        for ang = 1 : numel(azis)
            tic

            %HRTF filter signal and choose microphone channels
            [HRIR] = HRTFfilter(signal,azis(ang),kv,HRTF);

            if azis(ang) == 0
                  temp = HRIR(:,1)./rms(HRIR(:,1));
                  scalor = scaletodbspl(level,[],dboffset);
                  scalor = rms(temp.*scalor)/rms(HRIR(:,1));
            end

            HRIR = HRIR .* scalor;
            sigLengthSec = (size(HRIR,1)/fs);

            spikeRatePerBin = zeros(kv.numBin,2);

            for chan = 1 : 2

                 singChanSig = HRIR(:,chan);

                 [electrodogram, vTime] = ...
                                kelvasa2015_ciprocessing(singChanSig,...
                                kv.FS_ACE,'argimport',flags,kv);

                  [APvec] = ...
                                kelvasa2015_anprocessing(electrodogram,...
                                vTime,'argimport',flags,kv);

                  [~,spikeRatePerBinT] = ...
                                kelvasa2015_anbinning(APvec,...
                                sigLengthSec,'argimport',flags,kv);

                  spikeRatePerBin(:,chan) =  mean(spikeRatePerBinT,2);
            end

            SpkDiffPerBin(ang,:) = spikeRatePerBin(:,2)- spikeRatePerBin(:,1);
            SpkSumPerBin(ang,:) = spikeRatePerBin(:,2) + spikeRatePerBin(:,1);

            a = toc; timeLeft = round((a*(numLoops - n))/60);
            amt_disp(['Computing Figure 7. Time left (min):',num2str(timeLeft)],'volatile'); 
            n = n+1;
        end
            dataOut(sig).signalName = signalName{sig};
            dataOut(sig).levelDB = (level);
            dataOut(sig).SpkDiffPerBin = SpkDiffPerBin;
            dataOut(sig).SpkSumPerBin = SpkSumPerBin;
        end
        %Save data
         amt_disp();
         amt_cache('set',savename,dataOut);
        end
        if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 8a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8a

    %Parameters
       savename = ['Figure_8a_data_',kv.identifier];
       %Check for preprocessed calibration data
        [dataOut] = amt_cache('get', savename, flags.cachemode);
        if isempty(dataOut); clear dataOut

       levels = [45,55,65];
       signalName = 'SSN.wav';
       [signal fs] = amt_load('kelvasa2015',signalName);
       signal = signal(1:6*fs,:);
       signal = resample(signal,kv.FS_ACE,fs);
       kv.AziCrvfitRange =  45;
       kv.localizationModel =  'RateLevel';

      %Run model
       for ind = 1 : numel(levels)
            tic
            dataOut(ind) = groupPredictions(signal,levels(ind),kv,flags,HRTF, dboffset);

            a = toc; timeLeft = round((a*(3 - ind))/60);
            amt_disp(['Computing Figure 8a. Time left (min):',num2str(timeLeft)],'volatile');
       end

       %Save data
       amt_disp();
       amt_cache('set',savename,dataOut);
       end
       if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 8b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8b
       %Parameters
       savename = ['Figure_8b_data_',kv.identifier];

       %Check for preprocessed calibration data
        [dataOut] = amt_cache('get', savename, flags.cachemode);
        if isempty(dataOut); clear dataOut

       levels = [45,55,65];
       signalName = 'SSN.wav';
       [signal fs] = amt_load('kelvasa2015',signalName);
       signal = signal(1:6*fs,:);
       signal = resample(signal,kv.FS_ACE,fs);
       kv.AziCrvfitRange =  45;
       kv.localizationModel = 'ResponseDifferenceAN';

       %Run model
       for ind = 1 : numel(levels)
            tic
            %Run model
            [dataOut(ind)] = groupPredictions(signal,levels(ind),kv,flags,HRTF, dboffset);

            a = toc; timeLeft = round((a*(3 - ind))/60);
            amt_disp(['Computing Figure 8b. Time left (min):',num2str(timeLeft)],'volatile');
       end

       amt_disp();

       %Save data
       amt_cache('set',savename,dataOut);
       end
       if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 9a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig9a
       %Paramters
       savename = ['Figure_9a_data_',kv.identifier];

       %Check for preprocessed calibration data
        [dataOut] = amt_cache('get', savename, flags.cachemode);
        if isempty(dataOut); clear dataOut

       levels = [45,55,65];
       signalName = 'SSN.wav';
       [signal fs] = amt_load('kelvasa2015',signalName);
       signal = signal(1:6*fs,:);
       signal = resample(signal,kv.FS_ACE,fs);
       kv.AziCrvfitRange =  90;
       kv.localizationModel = 'ResponseDifferenceAN';


       %Run model
       for ind = 1 : numel(levels)
           tic

            %Run Model
            [dataOut(ind)] = groupPredictions(signal,levels(ind),kv,flags,HRTF, dboffset);

            a = toc; timeLeft = round((a*(3 - ind))/60);
            amt_disp(['Computing Figure 9a. Time left (min):',num2str(timeLeft)],'volatile');
       end

       %Save data
       amt_disp();
       amt_cache('set',savename,dataOut);
       end
       if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig10

       %Parameters
       savename = ['Figure_10_data_',kv.identifier];

       %Check for preprocessed calibration data
        [dataOut] = amt_cache('get', savename, flags.cachemode);
        if isempty(dataOut); clear dataOut

       level = 55;
       signalName = {'MixedSpeech.wav','FrozenSpeech.wav'};
       localizationModel = {'MaxLikelihood','ResponseDifferenceAN'};
       kv.AziCrvfitRange =  90;

       %Run model
       n = 1; 
       for model = 1 : 2
           tic
           if model == 1;
               kv.localizationModelCalibWav = 'MixedSpeech.wav';end
               kv.localizationModel = localizationModel{model};
           for sig = 1 : 2
               [signal fs] = amt_load('kelvasa2015',signalName{sig});
                signal = resample(signal,kv.FS_ACE,fs);
                [dataOut(n)] = groupPredictions(signal,level,kv,flags,HRTF, dboffset);

                a = toc; timeLeft = round((a*(4 - n))/60);
                amt_disp(['Computing Figure 10. Time left (min):',num2str(timeLeft)],'volatile'); 
                n = n+1;
           end

       end

       amt_disp();
       amt_cache('set',savename,dataOut);
       end
       if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig12

       %Parameters
       savename = ['Figure_12_data_',kv.identifier];
       %Check for preprocessed calibration data
        [dataOut] = amt_cache('get', savename, flags.cachemode);
        if isempty(dataOut); clear dataOut

       level = 55;
       signalName = 'SSN.wav';
       [signal fs] = amt_load('kelvasa2015',signalName);
       signal = signal(1:6*fs,:);
       signal = resample(signal,kv.FS_ACE,fs);
       Desired_SNR_dB = 5;
       kv.localizationModel = 'ResponseDifferenceAN';
       kv.AziCrvfitRange =  90;

        %Create 360� noise interferer
        azis = -180:5:175;
        [noise fs] = amt_load('kelvasa2015','PinkNoise.wav');
        noise = resample(noise,kv.FS_ACE,fs);
        noise = noise(1:numel(signal));

        for ang = 1:numel(azis)
            %HRTF filter signal and choose microphone channels
            [HRIR] = HRTFfilter(noise,azis(ang),kv,HRTF);
            noiseL(:,ang) = HRIR(:,1);
            noiseR(:,ang) = HRIR(:,2);
        end

        %Adjust signal to desired level. Referenced to 0� azi
           rmsSig = rms(mean(HRTFfilter(signal,0,kv,HRTF),2));
           rmsNoise = rms(mean([sum(noiseL,2),sum(noiseR,2)],2));

            K = (rmsSig/rmsNoise)*10^(-Desired_SNR_dB/20);  % Scale factor
            noise = [sum(noiseL,2),sum(noiseR,2)].*K;

        %Run Model on signals
        for ind = 1 : 2

            tic

            if ind == 1
                [dataOut(ind)] = groupPredictions(signal,level,kv,flags,HRTF, dboffset);

            else
                [dataOut(ind)] = groupPredictions(signal,level,kv,flags,...
                                                   HRTF,dboffset, noise);
            end

            a = toc; timeLeft = round((a*(2 - ind))/60);
            amt_disp(['Computing Figure 12. Time left (min):',num2str(timeLeft)],'volatile');
        end
        amt_disp();

       %Save data
       amt_cache('set',savename,dataOut);
       end
       if flags.do_plot; plot_kelvasa2015(dataOut,'argimport',flags,kv); end
end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataOut] = groupPredictions(signal,level,kv,flags,HRTF,dboffset, varargin)

           azis = 0:5:90;
           groupedBinPredictions{numel(kv.binPos)} = [];
           groupedWeightedPrediction= [];

       for ang = 1 : numel(azis)

            %HRTF filter signal and choose microphone channels
            [HRIR] = HRTFfilter(signal,azis(ang),kv,HRTF);

            %Adjust signal to desired level. Referenced to 0� azi
            if azis(ang) == 0
                  temp = HRIR(:,1)./rms(HRIR(:,1));
                  scalor = scaletodbspl(level,[],dboffset);
                  scalor = rms(temp.*scalor)/rms(HRIR(:,1));
            end
            HRIR = HRIR .* scalor;

            %If 360� interferer is to be added
            if ~isempty(varargin)
                noise = varargin{1}.*scalor;
                HRIR = HRIR + noise;

            end

            %Run HRTF filtered two channel signal through models
            [results] = kelvasa2015(HRIR,kv.FS_ACE,'argimport',flags,kv);

            %Bin azimuthal predictions over all target angles
            rangeAzis = [-10:5:90];
            if ~isempty(results.ANbinPredictions)
            for bin = 1 : numel(kv.binPos)
                data = results.ANbinPredictions(bin,:);
                [counts ind] = histc(data,rangeAzis);
                groupedPredictionsT = unique(rangeAzis(ind));
                countsT = counts(counts~=0);

                newInd  = [size(groupedBinPredictions{bin},2) + 1 : ...
                           size(groupedBinPredictions{bin},2) + ...
                                numel(groupedPredictionsT)];
                groupedBinPredictions{bin}(1,newInd) = groupedPredictionsT;
                groupedBinPredictions{bin}(2,newInd) = countsT;
                groupedBinPredictions {bin}(3,newInd) = ...
                                    repmat(azis(ang),1,numel(newInd));
            end;end

                data = results.weightedPredictions;
                [counts ind] = histc(data,rangeAzis);
                groupedPredictionsT = unique(rangeAzis(ind));
                countsT = counts(counts~=0);

                newInd  = (size(groupedWeightedPrediction,2) + 1) : ...
                          (size(groupedWeightedPrediction,2) + ...
                            numel(groupedPredictionsT));
                groupedWeightedPrediction(1,newInd) = groupedPredictionsT;
                groupedWeightedPrediction(2,newInd) = countsT;
                groupedWeightedPrediction(3,newInd) = ...
                                    repmat(azis(ang),1,numel(newInd));

        end

            dataOut.groupedBinPredictions = groupedBinPredictions;
            dataOut.groupedWeightedPrediction = groupedWeightedPrediction;
            groupedBinPredictions = [];

end

%%
function [outsig] = HRTFfilter(insig,ang,kv,HRTF)

            [~,ind_ang] = min(abs(HRTF.SourcePosition(:,1)-ang));

            HRIR = resample(squeeze(HRTF.Data.IR(ind_ang,:,:))',...
                    kv.FS_ACE,HRTF.Data.SamplingRate);

            %Filter signal with HRTF frequency domain
            HRTFchan1 = ifft(fft(insig).*fft(HRIR(:,1),numel(insig)));
            HRTFchan2 = ifft(fft(insig).*fft(HRIR(:,2),numel(insig)));
            outsig = [HRTFchan1,HRTFchan2];
end


