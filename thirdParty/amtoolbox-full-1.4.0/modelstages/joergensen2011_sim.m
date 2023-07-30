function [dSRT] = joergensen2011_sim(NSpeechsamples,varargin)
% JOERGENSEN2011_SIM Simulate the experiments shown in figure 5 and 6 of Jørgensen and Dau (2011)
%   Usage: [dSRT]  = joergensen2011_simFig6(NSpeechsamples,varargin)
%
%   Output parameters:
%     dSRT   : the delta-SRT between the reference condition and the
%              conditions with spectral subtraction processing
% 
%   The flag may be one of:
%     'fig5'       Simulate experiment with reverberation shown in fig5 of Jørgensen and Dau (2011)
%
%     'fig6'        Simulate experiment with spectral subtraction shown in fig6 of Jørgensen and Dau (2011)
%                   
%   
%
%   Please cite Jørgensen and Dau (2011) if you use
%   this model.
%
%   See also: joergensen2011, plot_joergensen2011
%
%   References:
%     S. Joergensen and T. Dau. Predicting speech intelligibility based on
%     the signal-to-noise envelope power ratio after modulation-frequency
%     selective processing. J. Acoust. Soc. Am., 130(3):1475--1487, 2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/joergensen2011_sim.php


%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Untrusted
%   #Requirements: M-Signal M-Stats
%   #Author: Peter L. Sondergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.flags.type = {'missingflag','fig5','fig6'};
definput.import={'amt_cache'};

[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% ------- simulate experiment with reverberation
if flags.do_fig5
  
    x=amt_load('joergensen2011','Danish_CLUE_10sentence_samples_22kHz.mat');
    sentenceArray=x.sentenceArray;
    
    result=amt_cache('get',['fig5 NSpeechsamples=' num2str(NSpeechsamples)],flags.cachemode);
    
    if isempty(result)

      amt_disp(['joergensen2011 started: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
      for q = 1:NSpeechsamples

          x = sentenceArray{q}';
          fs = 22050;

          sentenceFileLevel = -26.00; % The RMS level of all CLUE sentence files (in dB relative to 1)
          SPL = 65; % speech presentation level
          x = x*10^((SPL-sentenceFileLevel)/20);

          % the level of the sentence is set such that the long-term RMS of all
          % sentences are presented at a level of 65 dB SPL.

          Ts = 1/fs;
          T = length(x)/fs;
          t = 0:Ts:T;
          t = t(1:end-1);
          N = length(t);

          % load the noise file
          noise_glob = amt_load('joergensen2011','SSN_CLUE_22kHz.wav');

          Nsegments = floor(length(noise_glob)/N);
          % pick a random segment from the noise file
          startIdx = randi(Nsegments-2 ,1)*N;
          noise_glob = noise_glob(startIdx:startIdx+N -1)';

          SNRs = [ -9 -6 -3 0 3 6 9 ];
          conditions = [0 0.4 0.7 1.3 2.3];

          for n = 1:length(conditions)

              for k = 1:length(SNRs)

                  noise = noise_glob/rms(noise_glob)*10^((SPL-SNRs(k))/20);
                  if size(noise) ~= size(x)
                      noise = noise';
                  end
                  test = noise + x;

                  %     ------------ applying reverberation
                  if n>1

                    [tmp,Fs] = sig_joergensen2011(conditions(n));

                    tmp = resample(tmp,fs,Fs); % downsampling to 22.05 kHz
                    tmp_test = fconv(tmp',test); %
                    tmp_noise = fconv(tmp',noise);
                    test_env = abs(hilbert(tmp_test));

                    [bb, aa] = butter(4, 20*2/fs);
                    test_env = filter(bb,aa,test_env);
                    cut =  0.05;
                    threshold = floor(max(test_env)*cut);
                    %         finding the index in the env vector corresponding to the
                    %         threshold:
                    idx = find(test_env(floor(length(test_env)/4):end) > threshold, 1,'last' );

                    idx_end = idx + floor(length(tmp_test)/4);
                    test_env(idx_end);
                    idx_start = floor(0.047*fs);

                    test = tmp_test(idx_start:idx_end);
                    noise = tmp_noise(idx_start:idx_end);
                  end
                  % ---------------------------------

                  IOparameters = [0.82 0.5 8000 0.6]; % parameters from Jørgensen and Dau (2011).
                  tmp = joergensen2011(test,noise,fs,IOparameters);
                  SNRenvs(k,n,q) = tmp.SNRenv;
                  Pcorrect(k,n,q) = tmp.P_correct;


              end

          end
          amt_disp(['sentence nr: ' num2str(q) '/' num2str(NSpeechsamples)]);

      end

      result.Pcorrect = Pcorrect;
      result.SNRenvs = SNRenvs;
      result.conditions =conditions;
      result.SNRs =SNRs;
      amt_cache('set',['fig5 NSpeechsamples=' num2str(NSpeechsamples)],result);
    end
    
    %% Average across speech samples
    Pc_est_mean = mean(result.Pcorrect,3);
    
    
    % ------------------  Estimating changes in SRTs based on the mean Pcorrect -------------
    %   The first column of Pc_est_mean should always be the reference
    selection = 1:5;
    dSRT = joergensen2011_pctodsrt(Pc_est_mean,result.SNRs,selection);
end

if flags.do_fig6
    % Loads a cell array with 10 sentences from the CLUE material
    x=amt_load('joergensen2011','Danish_CLUE_10sentence_samples_22kHz.mat');
    sentenceArray=x.sentenceArray;
    
    amt_disp(['start: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
    
    for q = 1:NSpeechsamples
        
        x = sentenceArray{q}';
        fs = 22050;
        
        sentenceFileLevel = -26.00; % The RMS level of all CLUE sentence files (in dB relative to 1)
        SPL = 65; % speech presentation level
        x = x*10^((SPL-sentenceFileLevel)/20);
        
        % the level of the sentence is set such that the long-term RMS of all
        % sentences are presented at a level of 65 dB SPL.
        
        Ts = 1/fs;
        T = length(x)/fs;
        t = 0:Ts:T;
        t = t(1:end-1);
        N = length(t);
        
        % load the noise file
        noise_glob = amt_load('joergensen2011','SSN_CLUE_22kHz.wav');
        
        Nsegments = floor(length(noise_glob)/N);
        % pick a random segment from the noise file
        startIdx = randi(Nsegments-2 ,1)*N;
        noise_glob = noise_glob(startIdx:startIdx+N -1)';
        
        SNRs = [ -9 -6 -3 0 3 6 9 ];
        %
        conditions = [0 .5 1 2 4 8 ];
        
        for n = 1:length(conditions)
            
            for k = 1:length(SNRs)
                
                noise = noise_glob/rms(noise_glob)*10^((SPL-SNRs(k))/20);
                if size(noise) ~= size(x)
                    noise = noise';
                end
                test = noise + x;
                
                %        %% --------- spec sub parameters
                W=1024/2; % frame length
                padz=1024/2;
                %zero padding (pad with padz/2 from the left and padz/2 from the right )
                % Note that (W+padz) is the final frame window and hence the
                % fft length (it is normally chose as a power of 2)
                SP=0.5; %Shift percentage is 50%
                
                %     all this goes into a function:
                %% --------- spec sub
                %               Spectral subtraction is applied to both the noisy speech and the noise alone
                factor = conditions(n);
                [test2, noise2] = joergensen2011_sepsub(test,noise,W,padz,SP,factor);
                
                
                %%         ---------------------------------------------------------
                % It is important that the input to model is not zero-padded
                % of include silence (zeros) in the beginning and/or end of
                % the signal, since this leads to a large low-frequency
                % modulation, which will bias the model.  The SpecSub -
                % processing includes windowing which adds small periods of
                % silence in the beginning and end of the output.  The
                % output from the SpecSub is therefore truncated in the
                % beginning and the end. The truncation is inperceivable and
                % does not change the intelligibility of the sentence.
                test2_adj = test2(1.5*W:N);
                noise2_adj = noise2(1.5*W:N);
                
                
                
                IOparameters = [0.82 0.5 8000 0.6]; % parameters from Jørgensen et al., (2013).
                tmp = joergensen2011(test2_adj,noise2_adj,fs,IOparameters);
                SNRenvs(k,n,q) = tmp.SNRenv;
                Pcorrect(k,n,q) = tmp.P_correct;
                
                
            end
            
        end
        amt_disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
        
    end
    
%     s a v e Name = ['tmp_SNRenvs_specSub_' num2str(q) 'sent_' datestr(now, 'dd-mm-yyyy') '_1'];
    result.Pcorrect = Pcorrect;
    result.SNRenvs = SNRenvs;
    result.conditions =conditions;
    result.SNRs =SNRs;
%     s a v e (s a v e Name, 'result') % the mulktichannel SNRenv is stored.
        
    %% Average across speech samples
    Pc_est_mean = mean(result.Pcorrect,3);
    
    
    % ------------------  Estimating changes in SRTs based on the mean Pcorrect -------------
    %   The first column of Pc_est_mean should always be the reference
    selection = 1:6;
    dSRT = joergensen2011_pctodsrt(Pc_est_mean,result.SNRs,selection);
end



function [y]=fconv(h,x)
%FCONV Fast Convolution
%   [y] = FCONV(h, x) convolves x and h.  The output of this 
%         function is scaled.
%
%      x = input vector
%      h = input vector
% 
%      See also CONV
%
%Version 2.0
%2003-2004, Stephen G. McGovern

Ly=length(x)+length(h)-1;   % 
Ly2=pow2(nextpow2(Ly));     % Find smallest power of 2 that is > Ly
m= max(abs(x));      
X=fft(x, Ly2);		    % Fast Fourier transform
H=fft(h, Ly2);		    % Fast Fourier transform
Y=X.*H;        	           
y=real(ifft(Y, Ly2));       % Inverse fast Fourier transform
y=y(1:1:Ly);                % Take just the first N elements

m=m/max(abs(y));
y=m*y;


