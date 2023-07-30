function out_struct = hauth2020(MixedSig, fs, varargin)
%HAUTH2020 Blind equalization cancellation model
%
%   Usage: out_struct = hauth2020(MixedSig, fs)
%          out_struct = hauth2020('RequiredSignal',mixed_input,'OptionalSignal',...[OptionalSignal1 OptionalSignal2]...,OptionalSignalN,'model_params',model_params)
%
%   Input parameters:
%     MixedSig       : two-channel matrix containing signal + noise
%     fs             : the sampling frequency of MixedSignal [Hz]
%
%   Output parameters:
%     out_struct     : Structure containing the processed signals, in the order [Sel Min Max]
%
%
%   The SII used in the study by Hauth et al. (2020) 
%   requires band specific levels. If optional signals are present,
%   the levels of all optional signals after EC processing and better ear
%   selection are part of the output structure.
%
%   Additional input parameters:
%
%     'OptSigs',OptSig           matrix containing one or more optional two channel signals
%
%     'fmin',fmin                center frequency of the lowest filter of the gammatone filterbank [Hz]
%
%     'fmax',fmax                center frequency of the highest filter of the gammatone filterbank [Hz]
%
%     'f_target',f_target        specified frequency that will have a matched gammatone filter [Hz]
%
%     'ERB_factor',ERB_factor    gammatone filter bandwidth
%
%
%   HAUTH2020 accepts the following flags:
%
%     'binauralinaccuracies'     enable binaural processing inaccuracies
%
%     'no_binauralinaccuracies'  disable binaural processing inaccuracies
%
%     'longterm'                 SRMR-decision making considers whole signal
%
%     'shortterm'                blockwise SRMR-decision making based on same time constants as EC and BE processing
%
%
%   Parameter description:
%
%     'RequiredSignal'   is the mixture of speech and noise, which is a two
%                        channel (left,right) matrix (mixed_input)
%
%     'OptionalSignal'   One or more optional two channel signals, which are
%                        processed in the same way as the 'RequiredSignal'. They have to be
%                        arranged in a single matrix such that always two columns correspond to
%                        the left and right ear signals [left1 right1 left2 right2 ... leftN rightN]
%                        (e.g. if you need clean speech and noise for your back-end)
%
%     'model_params'     Structure containing model parameters like frequency range etc.
%
%     'out_struct'         Structure containing the processed signals, in the order [Sel Min Max]
%                          e.g. out_struct.signals.Mixsigsynfin contains the EC/BE processed mixed
%                          signal, where the best strategy was chosen blindly
%
%
%   See also: demo_hauth2020 hauth2020_fftcon hauth2020_sii hauth2020_ecprocess4optsigs
%             hauth2020_srmr hohmann2002
%
%   References:
%     C. F. Hauth, S. C. Berning, B. Kollmeier, and T. Brand. Modeling
%     binaural unmasking of speech using a blind binaural processing stage.
%     Trends in Hearing, 2020.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/hauth2020.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Christopher F. Hauth (2020)
%   #Author: Dr. Thomas Brand (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%--------------------------------------------------------------------------

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end

definput.import = {'hauth2020'}; % load defaults from arg_hauth2020
[flags,kv]  = ltfatarghelper({},definput,varargin);

fs_model            = 44100; % model sampling rate (must be 44100 Hz)

iNumofOptsigs = size(kv.OptSigs,2)/2;
frange = [kv.fmin kv.fmax];
BMtype = kv.Filterbank;
ERB_factor = kv.ERB_factor;

if flags.do_longterm
  batch_flag = 1;
else
  batch_flag = 0;
end
if flags.do_binauralinaccuracies
  bin_inaccuracy  = 1;
else
  bin_inaccuracy  = 0;
end  

if frange(1)==frange(2)
    ftarget         = frange(1);
else 
    ftarget         = kv.f_target; %[Hz]    
end
 % resmaple input if necessary
if fs_model~=fs
    MixedSig        = resample(MixedSig,fs_model,fs);
    if ~isempty(kv.OptSigs)
        kv.OptSigs     = resample(kv.OptSigs,fs_model,fs);
    end
end

% define temporal windows for EC processing and better ear selection
% batch processing means that the whole signal is considered as a single time frame
if batch_flag 
    iBlockLen_EC    = length(MixedSig); 
    iBlockLen_BE    = length(MixedSig);
    iNumofBlocks_EC = 1;
    iNumofBlocks_BE = 1;
else
    % If short time processing is performed, EC process and better ear
    % processing operate on different time constants.
    % The EC process is conducted in
    % time frames of 300ms (according to Hauth&Brand, 2018)
    % Using a frame shift of 50% leads to an effective temporal window of
    % 150 ms.
    iBlockLen_EC    = round(0.300.*fs_model); 
    % Better ear processing is performed on time frames of 20ms
    % Using a frame shift of 50% leads to an effective temporal window of
    % 10 ms.
    iBlockLen_BE    = round(0.020.*fs_model);
    
    % make sure the block length is even
    if mod(iBlockLen_EC,2)
        iBlockLen_EC = iBlockLen_EC+1;
    end
    iNumofBlocks_EC = 2.*floor(length(MixedSig)./iBlockLen_EC)-1;
    iNumofBlocks_BE = 2.*floor(length(MixedSig)./iBlockLen_BE)-1;
    iBlockshift_EC  = floor(iBlockLen_EC/2); % 50 percent overlap
    iBlockshift_BE  = floor(iBlockLen_BE/2);
    window_EC       = hann(iBlockLen_EC,'periodic');
    window_BE       = hann(iBlockLen_BE,'periodic');
end
% Binaural Processing Inaccuracies defined by vom H�vel (1984)
% For details see supplementary material of Hauth et al. (2020) (not yet published)
if bin_inaccuracy
    % error in level 
    sigma_epsilon0 = 1.5; % [dB]
    alpha0         = 15;  % [dB]
    p              = 1.6; % no unit as it used as exponent
    
    % error in time
    sigmadelta0    = 65.*10^-6; % [sec]
    Delta0         = 1.6.*10^-3;% [sec] 
else
    sigmadelta0    = 0;% [sec]
    Delta0         = 1;
end
% Parameters for the gammatone filterbank
order              = 4;
filter_per_ERB     = 1;

%parameter defining the cut-off frequency of EC processing. Above this
%frequency, better-ear processing is assumed to play the dominant role.
f_cut = 1500; % 1500 Hz

% Gammatone Filterbank from described in Hohmann (2002)
if strcmp(kv.Filterbank,'GT')
    analyzer1 = hohmann2002(fs_model,frange(1),ftarget,frange(2),filter_per_ERB,order,ERB_factor);
     %% Apply the Gammatone filterbank
     % LEFT EAR
     [MixsigLbatch, ~] = hohmann2002_process(analyzer1, MixedSig(:,1));
     fc = analyzer1.center_frequencies_hz;
     MixsigLbatch = real(MixsigLbatch);
     SomeImagSig = imag(MixsigLbatch);
     % RIGHT EAR
     [MixsigRbatch, ~] = hohmann2002_process(analyzer1, MixedSig(:,2));
     MixsigRbatch = real(MixsigRbatch);
     % If optional signals are present, process in the same way as mixed
     % signals

  if ~isempty(kv.OptSigs)
     OptSigs = kv.OptSigs.';
     if ~isempty(kv.OptSigs)
          for mm=1:iNumofOptsigs
                eval(sprintf('[OptSig%dLbatch, ~] = hohmann2002_process(analyzer1, OptSigs(2*(mm-1)+1,:)'');',mm));
                eval(sprintf('OptSig%dLbatch = (real(OptSig%dLbatch));',mm,mm));
                eval(sprintf('[OptSig%dRbatch, ~] = hohmann2002_process(analyzer1, OptSigs(2*mm,:)'');',mm));
                eval(sprintf('OptSig%dRbatch = (real(OptSig%dRbatch));',mm,mm));
                %eval(sprintf('[OptSig%dLbatchH, ~] = gfb_analyzer_process(analyzer1, OptSigs(:,2*(mm-1)+1));',mm));
                %eval(sprintf('OptSig%dLbatch = transpose(real(OptSig%dLbatchH));',mm,mm));               
                %eval(sprintf('[OptSig%dRbatchH, ~] = gfb_analyzer_process(analyzer1, OptSigs(:,2*mm));',mm));
                %eval(sprintf('OptSig%dRbatch = transpose(real(OptSig%dRbatchH));',mm,mm));
          end
                          OptSigs = OptSigs.';
     end
   end
end
% Allocate buffer for the final resynthesis of signals
% the ending fin denotes final representation of the signal
MixsigLfin    = zeros(size(MixsigLbatch));
MixsigRfin    = MixsigLfin;
Mixsigminfin  = MixsigLfin;
Mixsigmaxfin  = MixsigLfin;
Mixsigsynfin  = MixsigLfin;
% Allocate buffer for optional signals if present:
%if exist('OptSigs','var')
if ~isempty(kv.OptSigs)
    for ll=1:iNumofOptsigs
        eval(sprintf('OptSig%dLfin  = MixsigLfin;',ll));
        eval(sprintf('OptSig%dRfin  = MixsigLfin;',ll));
        eval(sprintf('OptSig%dminfin = MixsigLfin;',ll));
        eval(sprintf('OptSig%dmaxfin = MixsigLfin;',ll));
        eval(sprintf('OptSig%dsynfin = MixsigLfin;',ll));
    end
end
%% Binaural stage processing starts here:
% iterate through frequencies
for kk = 1:length(analyzer1.center_frequencies_hz)
    % Check, if frequency is in the range of frequencies where binaural
    % processing is assumed
    if analyzer1.center_frequencies_hz(kk)<=f_cut
        % Iterate through time blocks for EC processing
        for mm = 1:iNumofBlocks_EC  
            if iNumofBlocks_EC==1 % i.e. considering the signal as a single time frame
                % save signal portion 
                MixsigL = MixsigLbatch(:,kk);
                MixsigR = MixsigRbatch(:,kk);
                %if exist('OptSigs','var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch(:,kk);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch(:,kk);',ll,ll));
                    end
                end       
            else
                % Do block-wise processing 
                MixsigL = MixsigLbatch((mm-1)*iBlockshift_EC+1:...
                        iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);
                MixsigR = MixsigRbatch((mm-1)*iBlockshift_EC+1:...
                        iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);
                    % Do the same for optional signals if present
                    %if exist('OptSigs','Var')
                    if ~isempty(kv.OptSigs)
                        for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk).*sqrt(window_EC);',ll,ll));
                        end
                    end  
            end
            % Get levels of left and right ear signal
            rmsL = sqrt(mean(MixsigL.^2));
            rmsR = sqrt(mean(MixsigR.^2));
            % calculate mean between ears. This will be the target level
            % for the re-synthesis
            rms_target = mean([rmsL' rmsR'],2);
            
            % Level difference between left and right ear.
            alpha = 20.*log10(rmsL./rmsR);
            
            % Calculate level error (depending on the level difference)
            if bin_inaccuracy
                 sigma_epsilon = sigma_epsilon0.*(1+(abs(alpha)/alpha0).^p);
                 levErrorL = mean(10.^(sigma_epsilon.*randn(1,1)./20));
                 levErrorR = mean(10.^(sigma_epsilon.*randn(1,1)./20));
            else
                % if processing is assumed to be perfect:
                levErrorL = 1;
                levErrorR = 1;
            end
            % Calculate factors for level equalization:
            facR = sqrt((rmsL./rmsR).*levErrorR*levErrorL);
            facL = 1./facR;
           
            %% Apply equalization in Level (First step of EC Process)
            MixsigLproc = facL.*MixsigL;
            MixsigRproc = facR.*MixsigR;
            % Apply the same equalization to the Optional Signals:
            %if exist('OptSigs','var')
            if ~isempty(kv.OptSigs)
                for ll=1:iNumofOptsigs 
                    eval(sprintf('OptSig%dLproc=facL.*OptSig%dL;',ll,ll));
                    eval(sprintf('OptSig%dRproc=facR.*OptSig%dR;',ll,ll));
                end
            end
    %% Calculate delays to maximize and minimize the output of the EC
     % mechanism
           % run EC processing for the mixed signals:
           [Mixsigmin, Mixsigmax, ECparam4OptSigs] = hauth2020_fftcon(MixsigRproc,MixsigLproc,...
           analyzer1.center_frequencies_hz(kk),fs_model,sigmadelta0,Delta0,bin_inaccuracy);     
           
           %if exist('OptSigs','var')
           if ~isempty(kv.OptSigs)
               for ll=1:iNumofOptsigs
                    eval(sprintf('[OptSig%dmin,OptSig%dmax] = hauth2020_ecprocess4optsigs(OptSig%dRproc,OptSig%dLproc,ECparam4OptSigs,fs_model,bin_inaccuracy);',ll,ll,ll,ll));
               end
           end
 
        rms_max   = rms(Mixsigmax);
        Mixsigmax = rms_target./rms_max.*Mixsigmax;
        
        rms_min   = rms(Mixsigmin);
        Mixsigmin = rms_target./rms_min.*Mixsigmin;
        
        %if exist('OptSigs','var')
        if ~isempty(kv.OptSigs)
            for ll=1:iNumofOptsigs
                    eval(sprintf('OptSig%dmin = rms_target./rms_min.*OptSig%dmin;',ll,ll));
                    eval(sprintf('OptSig%dmax = rms_target./rms_max.*OptSig%dmax;',ll,ll));
            end
        end
        % Use the modified SRMR and apply it to both EC processed signals
        [ratiomax(mm,kk), ~] = hauth2020_srmr(Mixsigmax, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_EC);
        [ratiomin(mm,kk), ~] = hauth2020_srmr(Mixsigmin, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_EC);
           
        if mm>2
             if mean([ratiomax(mm,kk) ratiomax(mm-1,kk)])>mean([ratiomin(mm,kk) ratiomin(mm-1,kk)]);
                Mixsigsyn = Mixsigmax; 
                %if exist('OptSigs','Var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmax;',ll,ll));
                    end
                end
                amt_disp(['Max used in band ',num2str(kk)],flags.disp);
            else
                Mixsigsyn = Mixsigmin;
                %if exist('OptSigs','Var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmin;',ll,ll));
                    end
                end
                  amt_disp(['Min used in band ',num2str(kk)],flags.disp);
             end
         else
             
            if ratiomax(mm,kk)>ratiomin(mm,kk)
                Mixsigsyn = Mixsigmax; 
                %if exist('OptSigs','var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmax;',ll,ll));
                    end
                end
                  amt_disp(['Max used in band ',num2str(kk)],flags.disp);
            else
                Mixsigsyn = Mixsigmin;
                %if exist('OptSigs','var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dsyn = OptSig%dmin;',ll,ll));
                    end
                end
                  amt_disp(['Min used in band ',num2str(kk)],flags.disp);
            end
        end
         
        % As the SII was used in the study of Hauth et al.(2020), the levels of the optional
        % signals (if present) are saved
        %if exist('OptSigs','var')
        if ~isempty(kv.OptSigs)
            for ll=1:iNumofOptsigs
                eval(sprintf('LevelOptSig%dmin(kk) = 20.*log10(rms(OptSig%dmin));',ll,ll));
                eval(sprintf('LevelOptSig%dmax(kk) = 20.*log10(rms(OptSig%dmax));',ll,ll));
                eval(sprintf('LevelOptSig%dsyn(kk) = 20.*log10(rms(OptSig%dsyn));',ll,ll));
                eval(sprintf('LevelOptSig%dL(kk) = 20.*log10(rms(OptSig%dL));',ll,ll));
                eval(sprintf('LevelOptSig%dR(kk) = 20.*log10(rms(OptSig%dR));',ll,ll));
            end
        end
         
        % Construct final signals for resynthesis
         if iNumofBlocks_EC==1
            Mixsigminfin(:,kk) = Mixsigminfin(:,kk) + Mixsigmin;
            Mixsigmaxfin(:,kk) = Mixsigmaxfin(:,kk) + Mixsigmax;
            Mixsigsynfin(:,kk) = Mixsigsynfin(:,kk)+ Mixsigsyn;
            MixsigLfin(:,kk)   = MixsigLfin(:,kk) + (rms_target./rmsL).* MixsigL;
            MixsigRfin(:,kk)   = MixsigRfin(:,kk) + (rms_target./rmsR).*MixsigR;
            
            %if exist('OptSigs','var')
            if ~isempty(kv.OptSigs)
                for ll=1:iNumofOptsigs
                    eval(sprintf('OptSig%dminfin(:,kk) = OptSig%dminfin(:,kk) + OptSig%dmin;',ll,ll,ll));
                    eval(sprintf('OptSig%dmaxfin(:,kk) = OptSig%dmaxfin(:,kk) + OptSig%dmax;',ll,ll,ll));
                    eval(sprintf('OptSig%dsynfin(:,kk) = OptSig%dsynfin(:,kk) + OptSig%dsyn;',ll,ll,ll));
                    eval(sprintf('OptSig%dLfin(:,kk) = OptSig%dLfin(:,kk) +(rms_target./rmsL).* OptSig%dL;',ll,ll,ll));
                    eval(sprintf('OptSig%dRfin(:,kk) = OptSig%dRfin(:,kk) +(rms_target./rmsR).* OptSig%dR;',ll,ll,ll));
                end
            end
     
         else
                     
            Mixsigminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = Mixsigminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              Mixsigmin.*sqrt(window_EC);
        
            Mixsigmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = Mixsigmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              Mixsigmax.*sqrt(window_EC);
          
            Mixsigsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
           = Mixsigsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
             Mixsigsyn.*sqrt(window_EC);
           
            MixsigLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = MixsigLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              MixsigL.*sqrt(window_EC);
         
            MixsigRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)...
            = MixsigRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+...
              MixsigR.*sqrt(window_EC);
          
        % Do the same for optional signals if present
            %if exist('OptSigs','var')
            if ~isempty(kv.OptSigs)
                for ll=1:iNumofOptsigs
                    eval(sprintf('OptSig%dminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dminfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dmin.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dmaxfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dmax.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dsynfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dsyn.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dLfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dL.*sqrt(window_EC);',ll,ll,ll));
                    eval(sprintf('OptSig%dRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk) = OptSig%dRfin((mm-1)*iBlockshift_EC+1:iBlockLen_EC+(mm-1)*iBlockshift_EC,kk)+OptSig%dR.*sqrt(window_EC);',ll,ll,ll));
                end
            end
         end 
        end   
    else
        % Better ear processing:
        for mm = 1:iNumofBlocks_BE
            % Save left and right ear signal to temporary processing
            % variables.
            % if batch processing is used (whole signal is a single time
            % frame), save the whole signal.
            if iNumofBlocks_BE==1
                MixsigL = MixsigLbatch(:,kk);
                MixsigR = MixsigRbatch(:,kk);
                % Do the same for the optional signals.
                %if exist('OptSigs','var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch(:,kk);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch(:,kk);',ll,ll));
                    end
                end
            else
            % if short time processing is used, save only a signal portion.
                MixsigL = MixsigLbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                          (mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);
                MixsigR = MixsigRbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                          (mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);
                 % Do the same for the optional signals.
                %if exist('OptSigs','var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dL = OptSig%dLbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);',ll,ll));
                        eval(sprintf('OptSig%dR = OptSig%dRbatch((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk).*sqrt(window_BE);',ll,ll));
                    end
                end  
            end
            
            rmsL = sqrt(mean(MixsigL.^2));
            rmsR = sqrt(mean(MixsigR.^2));
            rms_target = mean([rmsL' rmsR'],2);
            % Use the modified SRMR and apply it to the left and right ear
            % signal to identify the better ear.
            [ratioL, ~] = hauth2020_srmr(rms_target./rmsL.*MixsigL, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_BE);
            [ratioR, ~] = hauth2020_srmr(rms_target./rmsR.*MixsigR, fs_model, 'fast', 0, 'norm', 1, 'minCF', 4, 'maxCF', 32, 'single',0,'window',iBlockLen_BE);
           
            % Independent of the EC processing strategy, the better ear is
            % always selected, and the worse ear is always neglected. 
            % Therefore, the better ear channels are combined with all EC
            % outputs. (Min+BE, Max+BE, Syn+BE).
            if ratioL>=ratioR
                  % if the left ear channel has stronger modulation cues, 
                  % select it set it to the target level.
                  Mixsigmin  = rms_target./rmsL.*MixsigL;
                  Mixsigmax  = Mixsigmin;
                  Mixsigsyn = Mixsigmin;
                  % Do the same for optional signals
                  %if exist('OptSigs','var')
                  if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs  
                        eval(sprintf('OptSig%dsyn = rms_target./rmsL.*OptSig%dL;',ll,ll));
                        eval(sprintf('OptSig%dmin = rms_target./rmsL.*OptSig%dL;',ll,ll));
                        eval(sprintf('OptSig%dmax = rms_target./rmsL.*OptSig%dL;',ll,ll));
                    end
                  end
                    amt_disp('Left Ear used', flags.disp);
            else
                  % if the right ear channel has stronger modulation cues, 
                  % select it and set it to the target level.
                  Mixsigmin  = rms_target./rmsR.*MixsigR;
                  Mixsigmax  = Mixsigmin;
                  Mixsigsyn = Mixsigmin;
                  % Do the same for optional signals
                  %if exist('OptSigs','var')
                  if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs 
                        eval(sprintf('OptSig%dsyn = rms_target./rmsR.*OptSig%dR;',ll,ll));
                        eval(sprintf('OptSig%dmin = rms_target./rmsR.*OptSig%dR;',ll,ll));
                        eval(sprintf('OptSig%dmax = rms_target./rmsR.*OptSig%dR;',ll,ll));
                    end
                  end
                    amt_disp('Right Ear used', flags.disp);
            end
            % Because the SII was used in the study by Hauth et al.(2020),
            % the frequency specific levels of the optional signals are
            % calculated. They can later be used by the SII.
            % If no optional signals are present, the levels are not
            % calculated.
            %if exist('OptSigs','var')
            if ~isempty(kv.OptSigs)
                for ll=1:iNumofOptsigs
                    eval(sprintf('LevelOptSig%dmin(kk) = 20.*log10(rms(OptSig%dmin));',ll,ll));
                    eval(sprintf('LevelOptSig%dmax(kk) = 20.*log10(rms(OptSig%dmax));',ll,ll));
                    eval(sprintf('LevelOptSig%dsyn(kk) = 20.*log10(rms(OptSig%dsyn));',ll,ll));
                    eval(sprintf('LevelOptSig%dL(kk)   = 20.*log10(rms(OptSig%dL));',ll,ll));
                    eval(sprintf('LevelOptSig%dR(kk)   = 20.*log10(rms(OptSig%dR));',ll,ll));
                end
            end
             % save the channels in the final output matrices.
            if iNumofBlocks_BE==1
                % Minimized output
                Mixsigminfin(:,kk)=  Mixsigminfin(:,kk) + Mixsigmin;
                % Maximized output
                Mixsigmaxfin(:,kk)  = Mixsigmaxfin(:,kk)+ Mixsigmax;
                % Synthesized output
                Mixsigsynfin(:,kk) = Mixsigsynfin(:,kk) + Mixsigsyn;
                % Left ear output
                MixsigLfin(:,kk)  =  MixsigLfin(:,kk)   +(rms_target./rmsL).*MixsigL;
                % Right ear output
                MixsigRfin(:,kk)  =  MixsigRfin(:,kk)   +(rms_target./rmsR).*MixsigR;
                % Do the same for optional signals.
                %if exist('OptSigs','var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dminfin(:,kk)  = OptSig%dminfin(:,kk) + OptSig%dmin;',ll,ll,ll));
                        eval(sprintf('OptSig%dmaxfin(:,kk) = OptSig%dmaxfin(:,kk) + OptSig%dmax;',ll,ll,ll));
                        eval(sprintf('OptSig%dsynfin(:,kk) = OptSig%dsynfin(:,kk) + OptSig%dsyn;',ll,ll,ll));
                        eval(sprintf('OptSig%dLfin(:,kk)    = OptSig%dLfin(:,kk) +(rms_target./rmsL).* OptSig%dL;',ll,ll,ll));
                        eval(sprintf('OptSig%dRfin(:,kk)    = OptSig%dRfin(:,kk) +(rms_target./rmsR).* OptSig%dR;',ll,ll,ll));
                    end
                end
                % save channels to final output matrices ( if short time processing was used)
            else 
                %Mixed Signals
                % Mimizized
                Mixsigminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = Mixsigminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+Mixsigmin.*sqrt(window_BE);
                % Maximimized
                Mixsigmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = Mixsigmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+ Mixsigmax.*sqrt(window_BE);
                % Synthesized
                Mixsigsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = Mixsigsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+ Mixsigsyn.*sqrt(window_BE);              
                 %Left
                MixsigLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = MixsigLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+MixsigL.*sqrt(window_BE);
                %Right
                MixsigRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)...
                = MixsigRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+...
                (mm-1)*iBlockshift_BE,kk)+MixsigR.*sqrt(window_BE);
           
                % Process Optional Signals accordingly:
                %if exist('OptSigs','var')
                if ~isempty(kv.OptSigs)
                    for ll=1:iNumofOptsigs
                        eval(sprintf('OptSig%dminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dminfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dmin.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dmaxfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dmax.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dsynfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dsyn.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dLfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dL.*sqrt(window_BE);',ll,ll,ll));
                        eval(sprintf('OptSig%dRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk) = OptSig%dRfin((mm-1)*iBlockshift_BE+1:iBlockLen_BE+(mm-1)*iBlockshift_BE,kk)+OptSig%dR.*sqrt(window_BE);',ll,ll,ll));
                    end
                end
            end
        end
    end
end
 %% Resynthesize time domain signals using gammatone filters and a delay of 10ms
synthesizer = hohmann2002_synth(analyzer1,0.01);%1/fs_model);%1/fs_model);
Mixsigminfin = hohmann2002_process(synthesizer, Mixsigminfin+SomeImagSig)';
Mixsigmaxfin= hohmann2002_process(synthesizer, Mixsigmaxfin)';
Mixsigsynfin = hohmann2002_process(synthesizer, Mixsigsynfin)';
MixsigLfin = hohmann2002_process(synthesizer,MixsigLfin+SomeImagSig)';
MixsigRfin = hohmann2002_process(synthesizer,MixsigRfin+SomeImagSig)';
% Do the same for optional signals
%if exist('OptSigs','var')
if ~isempty(kv.OptSigs)
    for ll=1:iNumofOptsigs
        eval(sprintf('OptSig%dminfin = hohmann2002_process(synthesizer, (OptSig%dminfin)+SomeImagSig)'';',ll,ll));
        eval(sprintf('OptSig%dmaxfin = hohmann2002_process(synthesizer, (OptSig%dmaxfin))'';',ll,ll));
        eval(sprintf('OptSig%dsynfin = hohmann2002_process(synthesizer, (OptSig%dsynfin))'';',ll,ll));
        eval(sprintf('OptSig%dLfin = hohmann2002_process(synthesizer, (OptSig%dLfin)+SomeImagSig)'';',ll,ll));
        eval(sprintf('OptSig%dRfin = hohmann2002_process(synthesizer, (OptSig%dRfin)+SomeImagSig)'';',ll,ll));
    end  
end

%% Resample signals if necessary
if fs_model~= fs
    Mixsigminfin  = resample(Mixsigminfin,fs,fs_model);
    Mixsigmaxfin  = resample(Mixsigmaxfin,fs,fs_model);
    Mixsigsynfin  = resample(Mixsigsynfin,fs,fs_model);
    MixsigLfin    = resample(MixsigLfin,fs,fs_model);
    MixsigRfin    = resample(MixsigRfin,fs,fs_model);
    
    %if exist('OptSigs','var')
    if ~isempty(kv.OptSigs)
        for ll=1:iNumofOptsigs
            eval(sprintf('[OptSig%dminfin, ~] = resample(OptSig%dminfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dmaxfin, ~] = resample(OptSig%dmaxfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dsynfin, ~] = resample(OptSig%dsynfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dLfin, ~]   = resample(OptSig%dLfin,fs,fs_model);',ll,ll));
            eval(sprintf('[OptSig%dRfin, ~]   = resample(OptSig%dRfinfs,fs_model);',ll,ll));
        end
    
    end
end
%% Save signals to output structure
out_struct.signals.MixsigSyn = Mixsigsynfin;
out_struct.signals.MixsigMin = Mixsigminfin;
out_struct.signals.MixsigMax = Mixsigmaxfin;
out_struct.signals.MixsigL = MixsigLfin;
out_struct.signals.MixsigR = MixsigRfin;
% Save also the optional signals to the output structure
    %if exist('OptSigs','var')
    if ~isempty(kv.OptSigs)
        for ll=1:iNumofOptsigs
            eval(sprintf('out_struct.signals.OptSig%dsyn = OptSig%dsynfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dsyn = LevelOptSig%dsyn;',ll,ll)); 
            eval(sprintf('out_struct.signals.OptSig%dmin = OptSig%dminfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dmin = LevelOptSig%dmin;',ll,ll));
            eval(sprintf('out_struct.signals.OptSig%dmax = OptSig%dmaxfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dmax = LevelOptSig%dmax;',ll,ll));
            eval(sprintf('out_struct.signals.OptSig%dL = OptSig%dLfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dL = LevelOptSig%dL;',ll,ll));
            eval(sprintf('out_struct.signals.OptSig%dR = OptSig%dRfin;',ll,ll));
            eval(sprintf('out_struct.levels.LevelOptSig%dR = LevelOptSig%dR;',ll,ll));
        end
    end
%--------------------Licence ---------------------------------------------
% Copyright (c) <2020> Christopher F. Hauth
% Dept. Medical Physics and Acoustics
% Carl von Ossietzky University Oldenburg 
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to 
% permit persons to whom the Software is furnished to do so, subject 
% to the following conditions:
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% END OF FILE


