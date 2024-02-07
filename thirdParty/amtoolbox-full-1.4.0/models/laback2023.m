function output = laback2023(insig, fs, varargin)
%LABACK2023   Contextual lateralization based on interaural level differences
%
%   Usage: [outsig, fc] = laback2023(insig,fs);
%          [outsig, fc] = laback2023(insig,fs,...);
%          [outsig, fc, params] = laback2023(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%
%   Output parameters:
%     outsig  : output acoustic signal.
%     fc      : filter bank center frequencies.
%  
%   LABACK2023(insig,fs) computes the internal representation of the
%   signal insig sampled with a frequency of fs Hz.
%  
%   [outsig,fc,mfc]=LABACK2023(...) additionally returns the center
%   frequencies of the filter bank and the center frequencies of the
%   modulation filterbank.
%  
%   References:
%     B. Laback. Contextual lateralization based on interaural level
%     differences is preshaped by the auditory periphery and predominantly
%     immune against sequential segregation. Trends in Hearing,
%     27:23312165231171988, 2023. PMID: 37161352. [1]arXiv | [2]www: ]
%     
%     C. J. Smalt, M. G. Heinz, and E. A. Strickland. Modeling the
%     Time-Varying and Level-Dependent Effects of the Medial Olivocochlear
%     Reflex in Auditory Nerve Responses. Journal of the Association for
%     Research in Otolaryngology, 15(2):159--173, Apr. 2014. [3]http ]
%     
%     M. S. A. Zilany, I. C. Bruce, and L. H. Carney. Updated parameters and
%     expanded simulation options for a model of the auditory periphery. The
%     Journal of the Acoustical Society of America, 135(1):283--286, Jan.
%     2014.
%     
%     References
%     
%     1. http://arxiv.org/abs/ https://doi.org/10.1177/23312165231171988
%     2. https://doi.org/10.1177/23312165231171988
%     3. https://doi.org/10.1007/s10162-013-0430-z
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/laback2023.php


%   #StatusDoc: Unknown
%   #StatusCode: Unknown
%   #Verification: Unknown
%   #Requirements: MATLAB MEX M-Signal
%   #Author: Bernhard Laback (2023): Original version 
%   #Author: Alejandro Osses (2023): adapted to AMT 1.4

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% ------ Checking of input parameters ------------

% load defaults
definput.import={'laback2023'}; 
% definput.importdefaults={''}; 
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

model     = keyvals.model;
precmode  = keyvals.precmode;
AdaptMode = keyvals.AdaptMode;
Modeldelay= keyvals.Modeldelay;
Onset     = keyvals.Onset;

N_samples = size(insig,1); % length of the time dimension
noiseILDvec = keyvals.noiseILDvec;
gapvec = keyvals.gapvec;

length_stimTarg = keyvals.length_stimTarg;
if isempty(length_stimTarg)
    if precmode == 0
        length_stimTarg = length(insig);
    else
        error('The length of the original stimTarg (length_stimTarg) needs to be specified');
    end
end

% ------ do the computation -------------------------
if flags.do_MOC
    MOCstat = 2;
end
if flags.do_no_MOC
    MOCstat = 1;
end
fiberType = keyvals.fiberType;

switch model
    case 1
        nrep = 20; % 20 repetitions
        dt = 1/fs;
        dur_sim = N_samples*1.5/fs;
        spont = keyvals.spontvec(fiberType);
        
        %%% For smalt2014, cohc is always 1
        keyvals.cohc=1; % maximum OHC activity 
        switch MOCstat
            case 1
                % MOC off
                mocr_max = 0; % minimum gain reduction (relative to COHC value)
            case 2
                % MOC on
                mocr_max=keyvals.mocr_max_def;   % maximum gain reduction (relative to COHC value)                            
        end
        if keyvals.BinMOC % binaural efferent model fit 2:1 dB ipsi/contra ratio
            mocr_threshold = -130.5076;  % 2:1
            mocr_slope = .0163840; % 2:1
            mocr_binauralratio = [0.4654    0.2332] ;% 2:1 => differs slightly from Tab. 1
        else %contralateral MOC only
            warning('Not validated yet...')
            mocr_threshold = -130.5076;  % 2:1
            mocr_slope = .0163840; % 2:1
            mocr_binauralratio = [0 0.2332] ;% 2:1 => differs slightly from Tab. 1
        end
        
        flags_z_fixed = {keyvals.CF_Model, nrep, dt, dur_sim, keyvals.cohc, keyvals.cihc, ...
            spont,mocr_max,mocr_threshold,mocr_slope,mocr_binauralratio,keyvals.shocks};
            
        %%% Part 1 ------------------------------------
        insig = transpose(insig); % needs to be a row array
        try
            [timeout,meout,mocr,c1filterout,c2filterout, ...
                c1vihc,c2vihc,vihc,synout,psth500k] = smalt2014(insig,flags_z_fixed{:});
        catch
            error('smalt2014 failed')
        end
        
        % numears = size(stim,1);   % 1=monaural, 2=binaural
        % psthbins = round(kv.psthbinwidth*fs);  % number of psth500k bins per psth bin
        %psthtime = timeout(1:psthbins:end); % time vector for psth
        %pr = squeeze(sum(reshape(psth500k,numears,psthbins,length(psth500k)/psthbins),2))/kv.nrep; % pr of spike in each bin
        %psth = pr/kv.psthbinwidth; % psth in units of spikes/s        

        %%% Part 2 ------------------------------------
        %   Obtains 'diffvec'

        %use synapse for determining internal ILD
        if precmode == 0
            wstart=Modeldelay;   
            wend = Modeldelay + length_stimTarg;

        elseif precmode == 1

            %Optimal
            wstart = Modeldelay + length(noiseILDvec) + length(gapvec);
            wend = Modeldelay + length(noiseILDvec) + length(gapvec) + length_stimTarg; 

        end
        MedWin=400; %window length for running average (smooting of response)

        diffpsth = medfilt1(synout(1,wstart:wend),MedWin) - medfilt1(synout(2,wstart:wend), MedWin);

        %use synapse output
        %Separate into rectangular windows, weight output in each
        %window  and integrate across windows
        WL = fs*0.002; %window duration (in ms)
        Overlap = 50; %overlap of windows in percent
        Shift = WL*Overlap/100; %shift between successive windows

        %number of windows (trailing portion (WNumb*(Shift)+Shift) ignored, if no integer number fits)
        WNumb = floor(length(synout(1,wstart:wend))/Shift - 1);
        WBeg = wstart;

        Wvec = zeros(2, WNumb);
        for i_wind = 1:WNumb             
            WEnd = WBeg+WL-1;
            W_L = mean(synout(1,WBeg:WEnd));  %mean rate in current L window
            W_R = mean(synout(2,WBeg:WEnd));  %mean rate in current R window
            Wvec(1, i_wind) = W_L;            %for analyzing monaural responses
            Wvec(2, i_wind) = W_R;            %for analyzing monaural responses
            MeanRate = (W_L + W_R)/2;         %needed for determining weight
            WBeg = WBeg+Shift;
        end
        %Without windowing and MED-filtering
        RespL = mean(synout(1,wstart:wend)); %wend
        RespR = mean(synout(2,wstart:wend)); %wend
        diffvec = RespL - RespR;
        diffvec_description = 'diffvec';
                                    
    case 2
        % definput.import = {'zilany2014'};
        % [fg_zilany, kv_zilany] = ltfatarghelper({},definput,varargin);
        
        %%% Part 1 ------------------------------------
        switch MOCstat
            case 1
                % MOC off
                keyvals.cohc=0; % no OHC activity %minimum gain reduction (relative to COHC value)
            case 2
                % MOC on
                keyvals.cohc=1; % maximum OHC activity %maximum gain reduction (relative to COHC value)
        end

        flags_z_common = {keyvals.CF_Model,'fiberType',fiberType, 'cohc', keyvals.cohc};
        
        % Left: Ipsilateral with 20 fibers, contralateral with 8 fibers:
        [synout_I_L,psth_I_L] = zilany2014(insig(:,1),fs, flags_z_common{:}, 'nrep', 20);   
        [synout_C_L,psth_C_L] = zilany2014(insig(:,1),fs, flags_z_common{:}, 'nrep', 8);   
        
        % Right: Ipsilateral with 8 fibers, contralateral with 20 fibers:
        [synout_I_R,psth_I_R] = zilany2014(insig(:,2),fs, flags_z_common{:}, 'nrep', 8);
        [synout_C_R,psth_C_R] = zilany2014(insig(:,2),fs, flags_z_common{:}, 'nrep', 20);

        disp('')
        % Contra LSO: Left(8), right(20);
        
        %%% Part 2 ------------------------------------

        if precmode == 0   
            switch AdaptMode
                case 0
                  wstart=Modeldelay + Onset*fs;   
                  wend = Modeldelay + N_samples;
                  
                case 1
                  wstart=Modeldelay + length(PrecPrecVec) + length(kv.gapvecLong) + Onset*Fs;   
                  wend = Modeldelay + length(PrecPrecVec) + length(kv.gapvecLong) + length_stimTarg;
            end
        elseif precmode == 1
            
            switch AdaptMode
                case 0
                  wstart = Modeldelay + length(noiseILDvec) + length(gapvec) + Onset*fs; %240 corresponding to processing delay => %adding 15000 samples reverts contra/ipsi effect (as does subtractin) => how possible.
                  wend   = Modeldelay + length(noiseILDvec) + length(gapvec) + length_stimTarg; 
                  
                case 1
                  wstart = Modeldelay + length(noiseILDvec)*2 + length(gapvec) + length(kv.gapvecLong) + Onset*fs; %240 corresponding to processing delay => %adding 15000 samples reverts contra/ipsi effect (as does subtractin) => how possible.
                  wend   = Modeldelay + length(noiseILDvec)*2 + length(gapvec) + length(kv.gapvecLong) + length_stimTarg; 
            end
        end
        
        MedWin=400; %window length for running average (smooting of response)        
        
        % Mean of the 'moving average':
        RespL_tmp = mean(medfilt1(synout_I_L(wstart:wend),MedWin));
        RespR_tmp = mean(medfilt1(synout_C_R(wstart:wend),MedWin));
        
        %Left LSO
        I_L=psth_I_L(wstart:wend);
        I_R=psth_I_R(wstart:wend);

        %Right LSO
        C_L=psth_C_L(wstart:wend);
        C_R=psth_C_R(wstart:wend);
        
        switch keyvals.ANmode 
        case 1
          RespL = RespL_tmp;
          RespR = RespR_tmp;
          
          % as used to predict LSO model based on Zilany2014
          % LSO model (Ashida, 2016)
          [spOut_I, ~] = ashida2016_LSOmodelCOC(I_L, I_R, 1/fs);
          [spOut_C, ~] = ashida2016_LSOmodelCOC(C_R, C_L, 1/fs);
          spOutVec_I = numel(find(spOut_I));
          spOutVec_C = numel(find(spOut_C));
          diffvec = spOutVec_I - spOutVec_C;
          diffvec_description = 'diffvec_LSO (Ashida 2016)';
          
        case 2 
          %   %used to predict spike-based AN model
          error('Not validated yey by AO')
          %   RespL(MOCstat, i_ILD) = mean(psth_I_L(wstart:wend));
          %   RespR(MOCstat, i_ILD) = mean(psth_I_R(wstart:wend));
          %   diffvec_LSO(MOCstat, i_ILD) = RespL(MOCstat, i_ILD) - RespR(MOCstat, i_ILD); 
           
        case 3
          % used to predict synapse-based AN model
          RespL = mean(psth_I_L(wstart:wend));
          RespR = mean(psth_I_R(wstart:wend));
          
          diffvec = RespL_tmp - RespR_tmp; 
          diffvec_description = 'diffvec_LSO';
        end
        
end

output.RespL = RespL;
output.RespR = RespR;
output.diffvec = diffvec;
output.diffvec_description = diffvec_description;
