function exp_mckenzie2022(varargin)
%EXP_mckenzie2022 experiments from McKenzie et al
%
%   Usage:
%     exp_mckenzie2022('fig11d'); % (to plot Figure 11 d).
%     exp_mckenzie2022('fig11a'); % (to plot Figure 11 a).
%
%   Reproduces the plots in the listening test section (Figure 11a-d)
%   McKenzie, T.; Armstrong, C.; Ward, L.; Murphy, D.T.; Kearney, G.
%   "Predicting the Colouration between Binaural Signals". Appl. Sci. 2022,
%   12(2441). https://doi.org/10.3390/app12052441
%
%   PEAQ and CLL are commented out in this script so that it runs without any
%   additional files necessary. To produce the respective data, download the
%   PEAQ and CLL code:
%   https://github.com/NikolajAndersson/PEAQ and
%   www.acoustics.hut.fi/-301ville/software/auditorymodel/.
%   The test sounds (ts.... etc) must be saved as wav files to work with the
%   PEAQ model.
%
%   Examples:
%   ---------
%
%   To display Figure 11a use :
%
%     exp_mckenzie2022('fig11a');
%
%   To display Figure 11d use :
%
%     exp_mckenzie2022('fig11d');
%
%   See also: mckenzie2022
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_mckenzie2022.php


%   #Author: Thomas McKenzie (2022)
%   #Author: Cal Armstrong (2022)
%   #Author: Lauren Ward (2022)
%   #Author: Damian Murphy (2022)
%   #Author: Gavin Kearney (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.flags.type={'missingflag','fig11a', 'fig11b', 'fig11c', 'fig11d'};
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

%% Read in listening test stimuli
data = amt_load('mckenzie2022', 'sig_mckenzie2022.mat');
% data = load('sig_mckenzie2022.mat');

fs = data.fs;
rs1 = data.rs1;
testDirections = data.testDirections;
ts1 = data.ts1;
ts2 = data.ts2;
ts3 = data.ts3;
ts4 = data.ts4;
ts5 = data.ts5;
ts6 = data.ts6;
ts7 = data.ts7;

% combine stimuli into one matrix
ts = cat(3,ts1,ts2,ts3,ts4,ts5,ts6,ts7);
rs = cat(3,rs1,rs1,rs1,rs1,rs1,rs1,rs1);
tsP = permute(ts,[1 3 2]);
rsP = permute(rs,[1 3 2]);

%% Read in listening test results
data = amt_load('mckenzie2022', 'data_mckenzie2022.mat');
% data = load('data_mckenzie2022.mat');

resultsRaw = data.resultsRaw;

% arrange results from repeated stimuli
resultsMean = round(squeeze(mean(resultsRaw))');

% without anchors
listTestResults = [resultsMean(:,1,:);resultsMean(:,2,:);resultsMean(:,3,:);resultsMean(:,4,:);resultsMean(:,5,:);resultsMean(:,6,:);resultsMean(:,7,:)];

%% Run spectral difference calculations
% parameters
freqRange = [20 20000]; nfft = length(rs(:,1,1));

%switch figureNumber
if flags.do_fig11a
        % BASIC SPECTRAL DIFFERENCE
        tsF = fftMatrix(tsP, fs, nfft, freqRange);
        rsF = fftMatrix(rsP, fs, nfft, freqRange);
        
        % get single values of spectral difference for all stimuli
        BSpecDiff = squeeze(mean(abs(tsF-rsF)));
        BavgSpecDiffS = mean(BSpecDiff,2);
        
elseif flags.do_fig11b
          amt_disp('perceptual evaluation of audio quality not implemented.');
%         % PERCEPTUAL EVALUATION OF AUDIO QUALITY
%         for i = 1:length(testDirectory)
%             [odg_tsA1(i), movb_tsA1(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'Ambi_1__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsA3(i), movb_tsA3(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'Ambi_3__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsA5(i), movb_tsA5(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'Ambi_5__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsD1(i), movb_tsD1(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'DFC_1__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsD3(i), movb_tsD3(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'DFC_3__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%             [odg_tsD5(i), movb_tsD5(:,i)] = PQevalAudio_fn(strcat(testDirectory,'HRIR__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'), strcat(testDirectory,'DFC_5__',num2str(testDirections(i,1)),'az_',num2str(testDirections(i,2)),'el_',num2str(i),'.wav'));
%         end
%         
%         % get single values of spectral difference for all stimuli
%         QavgSpecDiffS = [odg_tsA1 odg_tsA3 odg_tsA5 odg_tsD1 odg_tsD3 odg_tsD5]';
%         
elseif flags.do_fig11c
         amt_disp('composite loudness level not implemented.');
%         % COMPOSITE LOUDNESS LEVEL
%         for i = 1:length(tsP(1,:,1))
%             CLL_difference(i,:) = CLL(tsP(:,i,:),rsP(:,i,:),fs);
%         end
%         % get single values of spectral difference for all stimuli
%         CavgSpecDiffS = mean(abs(CLL_difference),2);
%         
elseif flags.do_fig11d
        % PREDICTED BINAURAL COLOURATION
        f.fs = fs;f.nfft = nfft;f.minFreq = freqRange(1); f.maxFreq = freqRange(2);
        [~,PSpecDiff] = mckenzie2022(tsP,rsP,0,f,0); %no fft pre model
        PSpecDiff = squeeze(PSpecDiff);
        
        % get single values of spectral difference for all stimuli
        PavgSpecDiffS = mean(PSpecDiff,2);
end

%% Calculate Pearson's Correlation Coefficient
% between spectral difference values and perceptual listening test results
%switch figureNumber
if flags.do_fig11a
        [rbsd, pbsd] = corrcoef(BavgSpecDiffS,listTestResults);
        disp(strcat('BSD correlation=',num2str(rbsd(2,1)),', p=',num2str(pbsd(2,1))));
        
% elseif flags.do_fig11b
%         [rqsd, pqsd] = corrcoef(QavgSpecDiffS,listTestResults);
%         disp(strcat('PEAQ correlation = ',num2str(rqsd(2,1)),', p = ',num2str(pqsd(2,1))));
%         
% elseif flags.do_fig11c
%         [rcsd, pcsd] = corrcoef(CavgSpecDiffS,listTestResults);
%         disp(strcat('CLL correlation = ',num2str(rcsd(2,1)),', p = ',num2str(pcsd(2,1))));
        
elseif flags.do_fig11d
        [rpsd, ppsd] = corrcoef(PavgSpecDiffS,listTestResults);
        disp(strcat('PBC correlation=',num2str(rpsd(2,1)),', p=',num2str(ppsd(2,1))));
end

%% Plot PCC vs Listening Test Results

plotColour = get(gca,'colororder');
plotColour(8,:) = [0.9,0.1,0.9];
xFit = linspace(min(listTestResults), max(listTestResults), 1000);

%switch figureNumber
if flags.do_fig11a
        j = 1; hold on;
        for i = 1:length(listTestResults)
            scatter(listTestResults(i),BavgSpecDiffS(i),60,'x','LineWidth',1.5,'MarkerEdgeColor',plotColour(j,:));
            j = j+1;
            if j == 9, j = 1; end
        end
        coefficients = polyfit(listTestResults, BavgSpecDiffS, 1);
        yFit = polyval(coefficients, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 1);
        ylabel('BSD (dB)'); xlabel('MUSHRA Test Results');
        set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
        pbaspect([1.7 1 1]); grid on; box on;
        
% elseif flags.do_fig11b
%         h = figure; j = 1; hold on;
%         for i = 1:length(listTestResults)
%             scatter(listTestResults(i),QavgSpecDiffS(i),60,'x','LineWidth',1.5,'MarkerEdgeColor',plotColour(j,:));
%             j = j+1;
%             if j == 9, j = 1; end
%         end
%         coefficients = polyfit(listTestResults, QavgSpecDiffS, 1);
%         yFit = polyval(coefficients, xFit);
%         plot(xFit, yFit, 'k-', 'LineWidth', 1);
%         ylabel('PEAQ (ODG)'); xlabel('MUSHRA Test Results');
%         set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
%         pbaspect([1.7 1 1]); grid on; box on;
%         
% elseif flags.do_fig11c
%         h = figure; j = 1; hold on;
%         for i = 1:length(listTestResults)
%             scatter(listTestResults(i),CavgSpecDiffS(i),60,'x','LineWidth',1.5,'MarkerEdgeColor',plotColour(j,:));
%             j = j+1;
%             if j == 9, j = 1; end
%         end
%         coefficients = polyfit(listTestResults, CavgSpecDiffS, 1);
%         yFit = polyval(coefficients, xFit);
%         plot(xFit, yFit, 'k-', 'LineWidth', 1);
%         ylabel('CLL (dB)'); xlabel('MUSHRA Test Results');
%         set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
%         pbaspect([1.7 1 1]); grid on; box on;
%         
%         legend('\theta = 180°, \phi = 64°',...
%             '\theta = 50°, \phi = 46°',...
%             '\theta = 118°, \phi = 16°',...
%             '\theta = 0°, \phi = 0°',...
%             '\theta = 180°, \phi = 0°',...
%             '\theta = 62°, \phi = −16°',...
%             '\theta = 130°, \phi = −46°',...
%             '\theta = 0°, \phi = −64°',...
%             'Location','NorthEast');
        
elseif flags.do_fig11d
        j = 1; hold on;
        for i = 1:length(listTestResults)
            scatter(listTestResults(i),PavgSpecDiffS(i),60,'x','LineWidth',1.5,'MarkerEdgeColor',plotColour(j,:));
            j = j+1;
            if j == 9, j = 1; end
        end
        coefficients = polyfit(listTestResults, PavgSpecDiffS, 1);
        yFit = polyval(coefficients, xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 1);
        ylabel('PBC (sones)'); xlabel('MUSHRA Test Results');
        set(gca,'FontSize', 14); set(gcf, 'Color', 'w');
        pbaspect([1.7 1 1]); grid on; box on;
end

%% Extra functions
    function [matrix_output_fft, freq_vector_fft,fft_abs_matrix_input] = fftMatrix(matrix_input, Fs, Nfft, freq_range)
        % Function to calculate the single sided frequency spectrum of two matrices
        % of HRIRs for a specified frequency range. Returns FFT of input matrix as
        % the absolute FFT in dB for the specified frequency range with the
        % associated frequency vector.
        
        % Take FFT of matrices
        fft_matrix_input = fft(matrix_input, Nfft); % Get Fast Fourier transform
        
        % Compute freq bins for x-axis limits
        fr_low = round(freq_range(1)*Nfft/Fs);
        fr_high = round(freq_range(2)*Nfft/Fs);
        
        % Get absolute values for frequency bins
        fft_abs_matrix_input = abs(fft_matrix_input(fr_low:fr_high,:,:));
        
        % Get values in dB
        matrix_output_fft = 20*log10(fft_abs_matrix_input);
        
        % Frequency vector for plotting
        f = 0:Fs/Nfft:Fs-(Fs/Nfft);
        freq_vector_fft = f(fr_low:fr_high);
    end

% Composite loudness level difference
    function [CLL_difference,freqs] = CLL(input,ref,Fs)
        if length(input(:,1)) < length(input(1,:))
            input = input';
        end
        if length(ref(:,1)) < length(ref(1,:))
            ref = ref';
        end
        
        % requires Karjalainen Auditory Toolbox.
        [~, ~, ~, CLL_input, freqs] = simuspatcues_KarAudMod(input(:,1),input(:,2),Fs);
        [~, ~, ~, CLL_ref] = simuspatcues_KarAudMod(ref(:,1),ref(:,2),Fs);
        
        CLL_difference = CLL_input-CLL_ref;
    end

end


