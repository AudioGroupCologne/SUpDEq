function [out,mpar] = eurich2022_processing(stimulus,mpar)
%eurich2022_processing  Binaural processing stage of the Eurich et al. (2022) model 
%   Usage: [out,mpar] = eurich2022_processing(stimulus,mpar)
%
%   Input parameters:
%     stim: Matrix containing stimulus (columns 1 and 2) and template (columns 3 and 4)
%     
%     mpar      : Structure with the following model parameters:
%
%                 - fs:                  sampling frequency of model
%                 - Filters_per_ERB_aud: spacing of peripheral filter central frequencies in ERB
%                 - flow:                lower bound of gammatone filterbank
%                 - fhigh:               higher bound of gammatone filterbank
%                 - gtorder:             filter order of gammatone filters (hohmann2002)
%                 - IPDnoise:            Internal noise affecting the IPD extraction
%                 - Dnoise:              Internal noise affecting the Decision stage
%                 - GaussSigma:          std of Gaussian across-channel smoothing
%                 - iKernelThresh:       treshold above which a value of the Gaussian filter window is used (below := 0)
%                 - n_hcbins:            number of bins to be used in IPD distribution computation via histcounts
%   
%   
%   Output parameters:
%     out:                       Processed output structure containing ???
%     
%     mpar:                      Updated structure with model parameters
%
%   See also: eurich2022 exp_eurich2022
%
%   References:
%     B. Eurich, J. Encke, S. D. Ewert, and M. Dietz. Lower interaural
%     coherence in off-signal bands impairs binaural detection. The Journal
%     of the Acoustical Society of America, 151(6):3927--3936, 06 2022.
%     [1]arXiv | [2]http ]
%     
%     References
%     
%     1. http://arxiv.org/abs/https://pubs.aip.org/asa/jasa/article-pdf/151/6/3927/16528275/3927\_1\_online.pdf
%     2. https://doi.org/10.1121/10.0011673
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/eurich2022_processing.php


%   #Author: Bernhard Eurich (2022): original implementation
%   #Author: Piotr Majdak (2023): adaptation to AMT 1.4

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% ==== Gammatone filtering ====

% filterbank object
sFB = hohmann2002(mpar.fs,mpar.GT_lowest_center_frequency,mpar.GT_fix_center_frequency,...
    mpar.GT_highest_center_frequency, mpar.GT_filters_per_ERBaud,'bandwidth_factor',mpar.GT_bandwidth_factor);

mpar.idx_of_500Hz_channel = find(round(sFB.center_frequencies_hz) == mpar.GT_fix_center_frequency);


% apply filterbank on stimulus
filtered_signal(:,:,1) = hohmann2002_process(sFB,stimulus(:,1));
filtered_signal(:,:,2) = hohmann2002_process(sFB,stimulus(:,2));

 

% ==== Binaural Stage ====
% complex correlation coefficient gamma for every filter channel
gamma = mean(filtered_signal(:,:,1) .* conj(filtered_signal(:,:,2))) ...
./ sqrt(mean(abs(filtered_signal(:,:,1)).^2) .* mean(abs(filtered_signal(:,:,2)).^2));



%  across-channel interference 
channel_index = 1:length(sFB.center_frequencies_hz); % to consider all channels

% for-loop represents convolution with window:
% * limit coherence with respect to one channel (index g)
% * create weighting window
% * apply window to coherence where channel g is center
% * --> for every channel --> same effect as convolution but now with coherence limitation applied with respect to each channel

% need odd number for interference window
flS = length(sFB.center_frequencies_hz);
flidx = mod(flS,2)<1;
mpar.KernelLen = flS - flidx;

if mpar.KernelLen == 1 | mpar.interference_sigma <= 0.11 % somewhat arbitrary --> decision whether interference happens
    interference_window = 1;
    interfered_coherence = abs(gamma);
    interfered_coherence_atanh_rho_max =  atanh(mpar.rho_max * interfered_coherence);
            
    
elseif find(mpar.idx_of_target_channel)
    Window_range = [-floor(mpar.KernelLen/2):ceil(mpar.KernelLen/2)-1];
    vExpWin0 = exp(-abs(Window_range)/(mpar.GT_filters_per_ERBaud*mpar.interference_sigma));
    vExpWin_ge   = vExpWin0(ge(vExpWin0,mpar.iKernelThresh));
    vExpWin_norm = vExpWin_ge ./ sum(vExpWin_ge);
    
    % interference of IPD fluctuations means only impact if off-frequency coherence lower than on-frequency coherence 
    coherence_lower = min(abs(gamma),abs(gamma(mpar.idx_of_500Hz_channel)));
    
    interference_range = mpar.idx_of_500Hz_channel - floor(length(vExpWin_norm)/2):mpar.idx_of_500Hz_channel + floor(length(vExpWin_norm)/2);
    interfered_coherence = sum(coherence_lower(interference_range) .* vExpWin_norm);

    interfered_coherence_atanh_rho_max = atanh(mpar.rho_max .* interfered_coherence);
    
    
else
    for g = channel_index
        
        % only lower off-frequency coherence is considered (compare with on-freq and take lower)
        coherence_lower = min(abs(gamma),abs(gamma(g)));
        
        zeros_padd = NaN(1, floor(length(coherence_lower)/2));
        coherence_lower_zeropadd = [zeros_padd coherence_lower zeros_padd];
        
        

        
        %     window for channel interaction
        Window_range = [-floor(mpar.KernelLen/2):ceil(mpar.KernelLen/2)-1];
        
        interference_window = NaN(1,2*length(Window_range)-1);
        window_idx1 = floor(length(interference_window)/2)-floor(length(Window_range)/2);
        window_idx2 = floor(length(interference_window)/2)+floor(length(Window_range)/2);
        
        interference_window(window_idx1:window_idx2) = ...
            exp(-abs(Window_range)/(mpar.GT_filters_per_ERBaud*mpar.interference_sigma));
        
        
        
        
        % place the window according to g (= temporary center channel)
        
        [~,center_of_window] = max(interference_window);
        interference_window_temp = interference_window(g:end-g);
        interference_window_areanorm = interference_window./nansum(interference_window);
        
        % actual interference happening here
        interfered_coherence = nansum(coherence_lower .* interference_window_areanorm(window_idx1:window_idx1+length(coherence_lower)-1));
        
        % select current interfered channel and apply internal noise (define maximum coherence)
        interfered_coherence_atanh_rho_max(g) =  atanh(mpar.rho_max * interfered_coherence);
    end
end



% write new IPD at f while the first 19 entries are empty
IPD = angle(gamma);

% binaural decision variable: complex correlation coefficient consisting of interfered coherence and (non-interfered) IPD
binaural_feature =  interfered_coherence_atanh_rho_max .* exp(1i*IPD);%



%% ==== MONAURAL PATHWAY ====

% Long-term DC of central channel: squared mean
monaural_feature = squeeze( mean( abs(filtered_signal(:,channel_index,1)),1) .^2 / 2);




% ==== ROUTE DECISION VALUES TO OUTPUT ====
out = cat(1,binaural_feature,monaural_feature);





end


