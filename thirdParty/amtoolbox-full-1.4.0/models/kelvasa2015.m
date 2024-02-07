function [results] = kelvasa2015(insig,fs,varargin)
%KELVASA2015 Azimuthal localization in cochlear-implant listeners
%   Usage: [results] = kelvasa2015(insig,fs,varargin); 
%
%   Input parameters:
%     insig               : Can be either [N x 2] two channel audio signal
%                           or results structure of preProcessed data
%     fs                  : sampling rate (Hz)
%     varargin            : structure with all parameters required for model. If
%                           this is not included, default paramters are loaded.
% 
%   Output parameters:
%     results             : A structure containing the processed electrodograms,
%                           AN spike times, and model predicted azimuthal locations
% 
%   KELVASA2015(insig,fs) implements the ACE signal processing strategy
%   upon the two channel input signal to produce bilateral electrodograms.
%   This is further processed through an electrode nerve interface to generate
%   spike times of a population of AN neurons. A chosen localization model from Kelvasa
%   and Dietz 2015 is then used to map the two channel (right and left) outputs
%   to a predicted azimuthal position.
%
%
%   The output structure results has the following fields:
%
%     'electrodogramCHAN1'    [NxM] matrix of CI electrode current output in mA(???) with N = number of CI electrodes and M = time
%
%     'APvecCHAN1'            [Nx2] matrix of [Nx1] indices of spiking AN fibers [Nx2] spike times in seconds
%
%     'electrodogramCHAN2'    same but for second channel
% 
%     'APvecCHAN2'            same but for second channel
%
%     'SpkSumPerBin'          [NxM] matrix of Right and Left Spike Rate differences in spikes per second with N = number of AN frequency bands and M = time bins
%
%     'SpkDiffPerBin'         same but for Right and Left spike rate differences
%
%     'ANbinPredictions'      [NxM] matrix of azimuthal angle bin predictions in degrees with N= number of AN frequency bands and M = time bins
%
%     'weightedPredictions'   [1xM] matrix of bin weighted azimuthal angle bin predictions in degrees with M = time bins
%
%     'mappingData'           Structure containing data used to calibrate and implement the chosen localization model as detailed in kelvasa2015calibratemodels.m
% 
%
%   The steps of the binaural model to calculate the result are the following :
% 
%   1 Process two channel input signal through a CI strategy as detailed in (Hamacher, 2003) and (Fredelake and Hohmann, 2012) to produce bilateral electrodograms.
%   
%   2 Process electrodogram through an electrode nerve interface and auditory nerve model as detailed in (Fredelake and Hohmann, 2012)
%                                                                                                                                                            
%   3 Compute bilateral spike rate differences over chosen AN frequency bands and time windows as detailed in (Kelvasa and Dietz, 2015)
%   
%   4 Calibrate the chosen localization model with a chosen calibration signal. This step can take several hours so preprocessed calibration is loaded for "Speech Shaped Noise" at 55dB as detailed in (Kelvasa and Dietz, 2015)
%   
%   5 Map the spike rate differences for each AN frequency band to a predicted azimuthal angle using the chosen localization model  as detailed in (Kelvasa and Dietz, 2015)
%
%
%   WARNING: If Octave is used, the plots will differ from the original ones, because
%   spectrogram is not available in Octave. Instead, sgram from LTFAT is used.
%
%
%   References:
%     S. Fredelake and V. Hohmann. Factors affecting predicted speech
%     intelligibility with cochlear implants in an auditory model for
%     electrical stimulation. Hearing Research, 287(1):76 -- 90, 2012.
%     [1]http ]
%     
%     V. Hamacher. Signalverarbeitungsmodelle des elektrisch stimulierten
%     Gehörs; 1. Aufl. PhD thesis, RWTH Aachen, Aachen, 2004. Zugl.: Aachen,
%     Techn. Hochsch., Diss., 2003.
%     
%     D. Kelvasa and M. Dietz. Auditory model-based sound direction
%     estimation with bilateral cochlear implants. Trends in Hearing,
%     19:2331216515616378, 2015.
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S0378595512000639
%     
%  
%   See also: kelvasa2015_anprocessing,
%             kelvasa2015_ciprocessing, kelvasa2015_anbinning,
%             kelvasa2015_calibratemapping, kelvasa2015_localize hohmann2002
%             plot_kelvasa2015
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/kelvasa2015.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal M-Stats
%   #Author: Daryl Kelvasa (2016)
%   #Author: Mathias Dietz (2016)
%   #Author: Clara Hollomey (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Retrieve and compute model paramters
    % Set flags
      definput.flags.disp = {'debug','no_debug'};
      definput.flags.plot = {'no_plot','plot'};
      definput.flags.plot_stage_fig = {'plot_stage_fig','no_plot_stage_fig'};
   
    % Import default arguments from other functions
    definput.import={'kelvasa2015','amt_cache'};
                        
    [flags,kv]  = ltfatarghelper({},definput,varargin);   

    % Check Inputs
    if nargin<2
          error('%s: Too few input parameters.',upper(mfilename));
    end;

    if isstruct(insig)
      results = insig; doPreProcessing = 0;
    else  doPreProcessing = 1;
    end

    if  isnumeric(insig), if min(size(insig))~=2
        error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
    end; end

    if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model processing starts here 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doPreProcessing

    %Signal resampling to CI processing
    if kv.FS_ACE ~= fs
      insig = resample(insig,kv.FS_ACE,fs);
      amt_disp('Resampling stimulus to ACE sampling frequency', flags.disp); 
    end

    results.sigLengthSec = (size(insig,1)/kv.FS_ACE);
    results.signal = insig;
    results.fs = kv.FS_ACE;
    results.modelParameters = kv;        
    
        for chan = 1 : 2
             
            signal = insig(:,chan);
            
            %CI Model adapted from (Hamacher, 2003) and 
            %(Fredelake and Hohmann, 2012)
            amt_disp('Running ACE signal processing strategy', flags.disp);    
            [electrodogram, vTime] = ...
                                kelvasa2015_ciprocessing(signal,...
                                kv.FS_ACE,'argimport',flags,kv);
            
            %CI/Auditory Nerve interface and  AN model adapted from 
            %(Fredelake and Hohmann, 2012)
            amt_disp('Running Electrode/Nerve interface and AN model', flags.disp);    
            [APvec] = ...
                                kelvasa2015_anprocessing(electrodogram,...
                                vTime,'argimport',flags,kv);

            if chan == 1
                results.vTime = vTime;
                results.electrodogramCHAN1 = electrodogram; 
                results.APvecCHAN1 = APvec; else
                results.electrodogramCHAN2 = electrodogram; 
                results.APvecCHAN2 = APvec; end  
          
         end

else
  amt_disp('Using preProcessed AN and CI data from struct input', flags.disp);
end
%Processing starts here if preProcessing has already been performed
          
            %AN frequency band binning of AN spike rate outputs as
            %described in Kelvasa and Dietz 2015
           amt_disp('Binning AN outputs into AN frequency bands along cochlea', flags.disp);
            
           [~,spikeRatePerBin1] = ...
                              kelvasa2015_anbinning(results.APvecCHAN1,...
                                results.sigLengthSec,'argimport',flags,kv);
                            
           [~,spikeRatePerBin2] = ...
                              kelvasa2015_anbinning(results.APvecCHAN2,...
                                results.sigLengthSec,'argimport',flags,kv);
                            
           results.SpkSumPerBin = spikeRatePerBin2+spikeRatePerBin1;
           results.SpkDiffPerBin = spikeRatePerBin2-spikeRatePerBin1;
           
           
           %Calibration of the linear and statistical localization models as
           %described in Kelvasa and Dietz 2015. 
           calibDataFilename = [kv.localizationModelCalibWav(1:end-4),'_',...
                   num2str(kv.localizationModelCalibStimulusLevelDB),...
                   'dB_calibration_', kv.identifier];
           %Check for preprocessed calibration data  
           [mappingData] = amt_cache('get', calibDataFilename, flags.cachemode);
           if isempty(mappingData)          
                amt_disp('Recomputing localization model calibration', flags.disp);
                [mappingData] =    kelvasa2015_calibratemapping('argimport',...
                                                            flags,kv);
                amt_cache('set',calibDataFilename,mappingData);
           end

       
               
          %Mapping of Right-Left spike rate for each bin to a predicted
          %azimuthal angle as described in Kelvasa and Dietz 2015
          amt_disp(['Implementing ',kv.localizationModel,...
                ' localization model to predict azimuthal source locations.'],flags.disp);
          [ANbinPredictions,...
               weightedPredictions,...
                     mappingData] = ...
                                      kelvasa2015_localize(mappingData, ...
                                                results.SpkDiffPerBin,...
                                                    results.SpkSumPerBin,...
                                                    'argimport',...
                                                            flags,kv); 
                                     
          results.ANbinPredictions = ANbinPredictions;
          results.weightedPredictions = weightedPredictions;
          results.mappingData = mappingData;
            
          if flags.do_plot_stage_fig
             plot_kelvasa2015mainFigure(results);
          end
end

%%
function plot_kelvasa2015mainFigure(results)
%% Plot main results
  
%Generate Main Figure
mainFig = figure('units','centimeters',...
                        'position',[0.5 0.5 21 17.5]);
 
rows = 4;
columns = 2;
xIndent = 3;
yIndent = 1.5;
ySpace = 1;
axisWidth = 8;
axisHeight = 3;
xTix = linspace(0,results.sigLengthSec,4);
xTixLab = [num2str(round(100.*(xTix'))./100),repmat(' s',numel(xTix),1)];

for axisInd = 2 : rows*columns
    
%Generate axis
    mapInd = [0 1 2 3 0 1 2 3];
    cornerX = xIndent + axisWidth*(ceil(axisInd/rows)-1) + ...
                        ceil(axisInd/rows)-1;
    cornerY = yIndent + mapInd(axisInd)*(axisHeight+ySpace);
    if axisInd == 5; cornerX = cornerX - axisWidth- 1; end
    
    ax(axisInd) = axes('Parent',mainFig,...
                       'Units','centimeters',...
                       'Fontsize', 14,...
                       'XTickLabel',[],...
                       'YTickLabel',[],...
                       'Position',...
                            [cornerX cornerY axisWidth axisHeight]);
                       hold(ax(axisInd),'on')           
                       axis(ax(axisInd),'manual') ;
   
%Plot data                            
    switch axisInd
      
        case {2,6}
            if axisInd == 2; 
                data = results.APvecCHAN1;           
                set(ax(axisInd),'YTick',[1,500],...
                            'YTickLabel',round([1,500].*(35/500)));
                if isoctave
                  ylabel({'Distance in','cochlea from Apex,','(mm)'});
                  xlabel({'Sorry: in Octave not available.'});        
                else 
                  ylabel({'Distance in','cochlea from Apex,','(mm)'});
                end
           
            else data = results.APvecCHAN2;
            end
          
            a = scatter(data(:,2),data(:,1));
if ~isoctave 
            set(a,          'MarkerFaceColor','black',...
                            'MarkerEdgeColor','black',... 
                            'SizeData',1);
                        
            set(ax(axisInd),'ylim',[1 500],...
                            'xlim',[0 max(results.sigLengthSec)],...
                            'XTick',xTix,...
                            'XTickLabel', xTixLab);
else
            amt_disp('Sorry, only limited display capabilities in Octave.');
            amt_disp('Cochlear distance not available.');
end                          
        case {3,7}
            elecCHAN1 = results.electrodogramCHAN1;
            elecCHAN2 = results.electrodogramCHAN2;
            elecRange = [min(min(elecCHAN1(:)), min(elecCHAN2(:))),...
                            max(max(elecCHAN1(:)), max(elecCHAN2(:)))];

            if axisInd == 3; 
                data = elecCHAN1;
                ylabel({'CI Electrode','Current,','(\mu A)'});
                set(ax(axisInd),...
                            'YTick',[1,6,11,17,22],...
                            'YTickLabel',[1,6,11,17,22]);
            else data = elecCHAN2;
            end
      
            imagesc(results.vTime,[],data,elecRange);
        
            set(ax(axisInd),'ylim',[1 22],...
                            'xlim',[0 max(results.vTime)],...
                            'XTickLabel',[])   
                        
        case {4,8}
            WINDOW = round(0.012 * results.fs);
            NOVERLAP = round(0.5*WINDOW); NFFT = 10000;
            
            if isoctave
              %originally, spectrogram was used, this is however, not
              %implemented in Octave, therefore, LTFATs sgram is used
              [S1,~,~] = sgram(results.signal(:,1),results.fs);
              [S2,~ ,~] = sgram(results.signal(:,2),results.fs); 

              S = [S1;S2];
              amt_disp('Octave: using LTFATs sgram instead of spectrogram.');
              amt_disp('Therefore, the output plot might not exactly match that from the publication.');
              F = linspace(0, results.fs/2, NFFT/2 + 1); 
              T = NOVERLAP/results.fs : NOVERLAP/results.fs: length(results.signal(:,2))/results.fs - NOVERLAP/results.fs;

            else
            
              [S1,~,~] = spectrogram(results.signal(:,1),WINDOW,...
                            NOVERLAP,NFFT,results.fs);
              [S2 F T] = spectrogram(results.signal(:,2),WINDOW,...
                    NOVERLAP,NFFT,results.fs); S = [S1;S2]; 
            end
                  
            colorRange = [min(20.*log10(abs(S(:)))),...
                                max(20.*log10(abs(S(:))))];
            
            if axisInd == 4; 
                data = S1;
                title('Chan 1 (Left)')
                ylabel({'Frequency,','(kHz)',''});
                set(ax(axisInd),...
                            'YTick',0:2000:8000,...
                            'YTickLabel',(0:2000:8000)./1000);
            else data = S2;
                title('Chan 2 (Right)')
            end  
            
            imagesc(T,F,20.*log10(abs(data)),colorRange)
          
            set(ax(axisInd),'ylim',[0 results.fs/2],...
                            'xlim',[0 max(T)],...
                            'XTickLabel',[])         
                    

        case 5
            %Set colors for each bin plot
            col = hsv(numel(results.modelParameters.binPos));    
            col(2,:) = col(2,:) + [0.2 -0.5 0];
            
            timeWin = results.modelParameters.timeWindowSec;
            timeBins = (1:numel(results.weightedPredictions)).*timeWin;
            numBins = numel(results.modelParameters.binPos);
            for bin = 1 : numBins + 1
                if ~isempty(results.ANbinPredictions) && bin <= numBins
                    a = scatter(timeBins,...
                        results.ANbinPredictions(bin,:),100);
                    
                    set(a,'MarkerFaceColor',col(bin,:),...
                          'Marker','o');
                else
                    a = scatter(timeBins,...
                        results.weightedPredictions,100);
                    
                    set(a,'MarkerFaceColor','black',...
                        'Marker','o');
                end
           end
            
           set(ax(axisInd),'ylim',[-10 100],...
                           'xlim',[0 results.sigLengthSec+timeWin],...
                           'YTickLabel',[],...
                           'YTick',[0 30 60 90],...
                           'XTick',timeBins,...
                           'YTickLabel',[0 30 60 90],...
                           'YGrid','on',...
                           'XGrid','on')
                       
          ylabel({'Predicted','Azimuthal','Angle, (�)'})
          xlabel(['Time Bins, ',num2str(timeWin),'s'])
          leg = cellstr([num2str(results.modelParameters.binPos'),...
                            repmat(' mm',numBins,1)]);
          leg{end+1} = 'weighted';
          a = legend(leg); b = get(a,'position'); b(1) = b(1)+4.5;
          b(2) = b(2)+0.25; set(a,'Position',b)
    end
end
end


