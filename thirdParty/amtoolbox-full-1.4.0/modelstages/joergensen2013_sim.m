function [SRTs conditions] = joergensen2013_sim(NSpeechsamples,varargin)
%JOERGENSEN2013_SIM Simulate the experiments shown in figure 2 of Jørgensen, Ewert and Dau (2013)
%   Usage: [SRTs conditions]  = joergensen2013_sim(NSpeechsamples,varargin)
%
%   Output parameters:
%     SRT   : the SRT as a function of the condition used in the experiment
%     conditions : the conditions used for a given experiment
%
%   The flag may be one of:
% 
%     'Jetal2013'           Simulate the conditions with SSN, SAM and ISTS using the
%                           CLUE speech materual, shown in fig2 of
%                           Jørgensen, Ewert and Dau (2013)
% 
%     'JandD2011reverb'     Simulate the conditions with reverberation
%                           shown in fig2 of Jørgensen, Ewert and Dau (2013)
%
%     'JandD2011specsub'    Simulate conditions with spectral subtraction
%                           shown in fig2 of Jørgensen, Ewert and Dau (2013)
%
%     'Kjems2009'           Simulate the data from Kjems et al. (2009)
%                           shown in fig2 of Jørgensen, Ewert and Dau (2013)
% 
%     'FP1990'              Simulate the data from Festen and Plomp
%                           (1990) shown in fig2 of Jørgensen, Ewert and Dau (2013)
%
%   
%
%   Please cite Joergensen et al. (2013) if you use
%   this model.
%
%   See also: joergensen2013, plot_joergensen2013
%
%   References:
%     S. Jørgensen, S. D. Ewert, and T. Dau. A multi-resolution envelope
%     power based model for speech intelligibility. J. Acoust. Soc. Am.,
%     134(1):436--446, 2013.
%     
%     U. Kjems, J. B. Boldt, M. S. Pedersen, T. Lunner, and D. Wang. Role of
%     mask pattern in intelligibility of ideal binary-masked noisy speech. J.
%     Acoust. Soc. Am., 126:1415--1426, 2009.
%     
%     J. Festen and R. Plomp. Effects of fluctuating noise and interfering
%     speech on the speech-reception threshold for impaired and normal
%     hearing. J. Acoust. Soc. Am., 88(4):1725--1736, 1990.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/joergensen2013_sim.php


%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Untrusted
%   #Requirements: MATLAB M-Signal M-Stats
%   #Author: Peter L. Sondergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

disp('**********************************************************');
disp('AMT WARNING: incomplete file --- missing: original stimuli');
disp('**********************************************************');
 
definput.flags.type = {'JandD2011specsub','JandD2011reverb','FP1990','Kjems2009','Jetal2013'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


% ------- simulate experiment with reverberation
if flags.do_Jetal2013
    %in the absence of the original stimuli - using those from joergensen2011 instead
    data = amt_load('relanoiborra2019', 'single_150_SentArray22kHz_varLength.mat');
    sentenceArray = data.sentenceArray;
    
    amt_disp(['start Jetal2013: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
    sentenceFileLevel = -26;
    %noise_names = {'SSN_CLUE_22kHz.wav','SSN_MOD_CLUE','ISTS_eq'};stimuli
    %not available
    noise_names = 'SSN_CLUE_22kHz.wav';
    speechSPL = 65;
    fixedNoiseLevel = 0;
    SNRs = -27:3:3;
    conditions = {'SSN'  'SAM' 'ISTS'};
    
    IOparameters = [0.61 0.5 8000 0.6]; %
    
    for q = 1:NSpeechsamples
        
        x = sentenceArray{q}';
        fs = 22050;
        x = x*10^((speechSPL-sentenceFileLevel)/20);
        
        % the level of the sentence is set such that the long-term RMS of all
        % sentences are presented at a level of 65 dB SPL.
        
        Ts = 1/fs;
        T = length(x)/fs;
        t = 0:Ts:T;
        t = t(1:end-1);
        N = length(t);
        
        % load the noise files
        for k = 1:length(noise_names)
            
            clear tmp tmp2

            %[tmp, fs] = audioread (noise_names{k});
            [tmp, fs] = amt_load('relanoiborra2019', noise_names);
            %if fs_tmp ~= fs
            %    noise_scaled{k} = resample(tmp,fs,fs_tmp);
            %else
                noise_scaled{k} = tmp;
            %end
            
            Nsegments = floor(length(noise_scaled{k})/N);
            % pick a random segment
            startIdx = randi(Nsegments-2 ,1)*N;
            
            noise_glob{k} = noise_scaled{k}(startIdx:startIdx+N -1);
        end
        
        for n = 1:length(conditions)
            
            for k = 1:length(SNRs)
                noise = noise_glob{n};
                noise = noise/rms(noise)*10^((speechSPL-SNRs(k))/20);
                if size(noise) ~= size(x)
                    noise = noise';
                end
                test = noise + x;
                
                tmp = joergensen2013(test,noise,fs,IOparameters);
                SNRenvs(k,n,q) = tmp.SNRenv;
                Pcorrect(k,n,q) = tmp.P_correct;
                
            end
            
        end
        amt_disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
        
    end
  
    result.Pcorrect = Pcorrect;
    result.SNRenvs = SNRenvs;
    result.conditions =conditions;
    result.SNRs =SNRs;

    
    
    
    %% Average across speech samples
    Pc_est_mean = mean(result.Pcorrect,3);
    
    
    % ------------------  Estimating changes in SRTs based on the mean Pcorrect -------------
    %   The first column of Pc_est_mean should always be the reference
    selection = 1:length(conditions);
    [dSRT SRTs] = joergensen2011_pctodsrt(Pc_est_mean,result.SNRs,selection);
end

if flags.do_JandD2011specsub
    %in the absence of the original stimuli - using those from joergensen2011 instead
    data = amt_load('relanoiborra2019', 'single_150_SentArray22kHz_varLength.mat');
    sentenceArray = data.sentenceArray;
    
    amt_disp(['start JandD2011specsub: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
    sentenceFileLevel = -26;
    noise_names = 'SSN_CLUE_22kHz.wav';
    speechSPL = 65;
    SNRs = -9:3:9;
    conditions = [0 0.5 1 2 4 8];
    
    IOparameters = [0.61 0.5 8000 0.6]; %
    
    for q = 1:NSpeechsamples
        k = 1;
        x = sentenceArray{q}';
        fs = 22050;
        x = x*10^((speechSPL-sentenceFileLevel)/20);
        
        % the level of the sentence is set such that the long-term RMS of all
        % sentences are presented at a level of 65 dB SPL.
        
        Ts = 1/fs;
        T = length(x)/fs;
        t = 0:Ts:T;
        t = t(1:end-1);
        N = length(t);
        
        % load the noise file
        
        clear tmp
        %[tmp, fs] = audioread (noise_names);
        [tmp, fs] = amt_load('relanoiborra2019', noise_names);
        %if fs_tmp ~= fs
        %    noise_scaled = resample(tmp,fs,fs_tmp);
        %else
            noise_scaled = tmp;
        %end
        
        Nsegments = floor(length(noise_scaled)/N);
        % pick a random segment
        startIdx = randi(Nsegments-2 ,1)*N;
        noise_glob = noise_scaled(startIdx:startIdx+N -1);
        
        
        for n = 1:length(conditions)
            
            for k = 1:length(SNRs)
                noise = noise_glob;
                noise = noise/rms(noise)*10^((speechSPL-SNRs(k))/20);
                if size(noise) ~= size(x)
                    noise = noise';
                end
                test = noise + x;
                
                %                     --------- spec sub -----------------
                W=1024/2; % frame length
                padz=1024/2; %zero padding (pad with padz/2 from the left and padz/2 from the right )
                % % Note that (W+padz) is the final frame window and hence the fft length (it is normally chose as a power of 2)
                SP=0.5; %Shift percentage is 50%
                
                factor = conditions(n);
                ProcMix  = joergensen2011_specsub(test,noise,W,padz,SP,factor);
                
                ProcMixInv = joergensen2011_specsub(x-noise,-noise,W,padz,SP,factor);
%                 Estimating the noise alone using the approach by Hagerman and Olofsson (2004)
                NoiseEst = (ProcMix  - ProcMixInv)/2;
                
                test = ProcMix(1.5*W:N); %cutting off first and last bit
                
                noise = NoiseEst(1.5*W:N);
                %                      -----------------------------
                
                tmp = joergensen2013(test,noise,fs,IOparameters);
                SNRenvs(k,n,q) = tmp.SNRenv;
                Pcorrect(k,n,q) = tmp.P_correct;
                
            end
            
        end
        amt_disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
        
    end
    
    result.Pcorrect = Pcorrect;
    result.SNRenvs = SNRenvs;
    result.conditions =conditions;
    result.SNRs =SNRs;

    
    
    
    %% Average across speech samples
    Pc_est_mean = mean(result.Pcorrect,3);
    
    
    % ------------------  Estimating changes in SRTs based on the mean Pcorrect -------------
    %   The first column of Pc_est_mean should always be the reference
    selection = 1:length(conditions);
    [dSRT SRTs] = joergensen2011_pctodsrt(Pc_est_mean,result.SNRs,selection);
end

if flags.do_JandD2011reverb
    %in the absence of the original stimuli - using those from joergensen2011 instead
    data = amt_load('relanoiborra2019', 'single_150_SentArray22kHz_varLength.mat');
    sentenceArray = data.sentenceArray;    
    
    amt_disp(['start JandD2011reverb: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
    sentenceFileLevel = -26;
    noise_name = 'SSN_CLUE_22kHz';
    speechSPL = 65;
    SNRs = -9:3:9;
    conditions = [0 0.4 0.7 1.3 2.3];
    
    IOparameters = [0.61 0.5 8000 0.6]; %
    
    for q = 1:NSpeechsamples
        
        x = sentenceArray{q}';
        fs = 22050;
        x = x*10^((speechSPL-sentenceFileLevel)/20);
        
        % the level of the sentence is set such that the long-term RMS of all
        % sentences are presented at a level of 65 dB SPL.
        
        Ts = 1/fs;
        T = length(x)/fs;
        t = 0:Ts:T;
        t = t(1:end-1);
        N = length(t);
        
        % load the noise file
        noise_names = 'SSN_CLUE_22kHz.wav';
        clear tmp
        %[tmp, fs] = audioread (noise_names{k});
        [tmp, fs] = amt_load('relanoiborra2019', noise_names);
%        if fs_tmp ~= fs
%            noise_scaled = resample(tmp,fs,fs_tmp);
%        else
            noise_scaled = tmp;
%        end
        
        Nsegments = floor(length(noise_scaled)/N);
        % pick a random segment
        startIdx = randi(Nsegments-2 ,1)*N;
        noise_glob = noise_scaled(startIdx:startIdx+N -1);
        
        
        for n = 1:length(conditions)
            
            for k = 1:length(SNRs)
                noise = noise_glob;
                noise = noise/rms(noise)*10^((speechSPL-SNRs(k))/20);
                if size(noise) ~= size(x)
                    noise = noise';
                end
                test = noise + x;
                
                %              %     ------------ applying reverberation
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
                
                tmp = joergensen2013(test,noise,fs,IOparameters);
                SNRenvs(k,n,q) = tmp.SNRenv;
                Pcorrect(k,n,q) = tmp.P_correct;
                
            end
            
        end
        amt_disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
        
    end
    
    result.Pcorrect = Pcorrect;
    result.SNRenvs = SNRenvs;
    result.conditions =conditions;
    result.SNRs =SNRs;

    
    
    
    %% Average across speech samples
    Pc_est_mean = mean(result.Pcorrect,3);
    
    
    % ------------------  Estimating changes in SRTs based on the mean Pcorrect -------------
    %   The first column of Pc_est_mean should always be the reference
    selection = 1:length(conditions);
    [dSRT SRTs] = joergensen2011_pctodsrt(Pc_est_mean,result.SNRs,selection);
end
%AMT: Stimuli for these last two not available
 if flags.do_FP1990
     amt_disp('amt warning: not enough data to calculate this.');
          SRTs = []; 
     conditions = [];
%     load PlompMimpen_130_SentArray22kHz
%     amt_disp(['start FP1990: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
%     sentenceFileLevel = 10^((-17.93)/20);
%     noise_names = {'MaleS00.wav','MaleF00.wav','P&MRunningSpeech_FemaleTR.wav'};
%     conditions = {'SSN','SMN','RT'};
%     SNRs = [  -21 -18 -12 -9 -6  -3 0 3];
%     noisedBA = 80;
%     fs = 22050;
%     IOparameters = [0.36 0.5 8000 0.6]; %
%     
%     clear noise_scaled
%     for k = 1:length(noise_names)
%         clear tmp tmp2
%         [tmp, fs] = audioread (noise_names{k});
%         tmp = resample(tmp,fs,fs_tmp);
%         tmp = tmp/rms(tmp);
%         noise_scaled{k} = Leq2dBA(tmp,fs,noisedBA); %here, all noises have same level dBA
%         noiseSPL_dBA(k) = LeqdBA(noise_scaled{k},fs);
%     end
%     
%    
%     
%     for q = 1:NSpeechsamples
%         
%         x_unscaled = sentenceArray{q}';
%         
%         
%         Ts = 1/fs;
%         T = length(x_unscaled)/fs;
%         t = 0:Ts:T;
%         t = t(1:end-1);
%         N = length(t);
%         
%         % load the noise files
%         for k = 1:length(noise_names)
%             
%             Nsegments = floor(length(noise_scaled{k})/N);
%             % pick a random segment
%             startIdx = randi(Nsegments-2 ,1)*N;
%             
%             noise_glob{k} = noise_scaled{k}(startIdx:startIdx+N -1);
%         end
%         
%         
%         
%         for n = 1:length(conditions)
%             
%             for k = 1:length(SNRs)
%                 noise = noise_glob{n};
%                 speech_dBA = (noisedBA+SNRs(k));
%                 SPL_sent(k) = 20*log10(rms(Leq2dBA(noise_scaled{1},fs,speech_dBA))); %#check!
%                 SpeechGain = 10^(SPL_sent(k)/20);
%                 x = x_unscaled ./ sentenceFileLevel *SpeechGain;
%                 SpeechSPL_dBA(k) = LeqdBA(x,fs);
%                 
%                 if size(noise) ~= size(x)
%                     noise = noise';
%                 end
%                 test = noise + x;
%                 
%                 
%                 tmp = joergensen2013(test,noise,fs,IOparameters);
%                 SNRenvs(k,n,q) = tmp.SNRenv;
%                 Pcorrect(k,n,q) = tmp.P_correct;
%                 
%             end
%             
%         end
%         amt_disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
%         
%     end
%     
%     result.Pcorrect = Pcorrect;
%     result.SNRenvs = SNRenvs;
%     result.conditions =conditions;
%     result.SNRs =SNRs;
% 
%     
%     
%     
%     %% Average across speech samples
%     Pc_est_mean = mean(result.Pcorrect,3);
%     
%     
%     % ------------------  Estimating changes in SRTs based on the mean Pcorrect -------------
%     %   The first column of Pc_est_mean should always be the reference
%     selection = 1:length(conditions);
%     [dSRT SRTs] = joergensen2011_pctodsrt(Pc_est_mean,result.SNRs,selection);
%     
 end

 if flags.do_Kjems2009
     amt_disp('amt warning: not enough data to calculate this.');
     SRTs = []; 
     conditions = [];
%      load DANTALE2_144single_SentArray44kHz_varLength
%      amt_disp(['start Kjems2009: ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
%             sentenceFileLevel = 10^((-19.79)/20);
%             noise_names = {'ssn_noise_20k','cafe_noise_20k','car_noise_20k','Bottle_noise_20k'};
%             fs = 22050;
%             noisedBA = 65;
%             SNRs = -30:3:10;
%             conditions = {'SSN','Cafe','Car','Bottle'};
%             
%             IOparameters = [0.42 0.5 50 0.9]; %
%             clear noise_scaled
%             
%             for k = 1:length(noise_names)
%                 [tmp, fs] = audioread (noise_names{k});
%                 tmp = resample(tmp,fs,fs_tmp);
%                 tmp = scaletodbspl(tmp,0);
%                 if k == 1
%                     [tmp GaindB]= Leq2dBA(tmp,fs,noisedBA); % the level in dBA is set for the SSN to 65
%                 end
%                 noise_scaled{k} = tmp/rms(tmp).*10^(GaindB/20); %here, all noises have same level in dB, thus different levels dBA!
%                 noiseSPL_dBA(k) = LeqdBA(noise_scaled{k},fs);
%             end
%       
%     for q = 1:NSpeechsamples
%         
%         x_unscaled = sentenceArray{q}';
%         x_unscaled = resample(x_unscaled,fs,44100);
%         
%         Ts = 1/fs;
%         T = length(x_unscaled)/fs;
%         t = 0:Ts:T;
%         t = t(1:end-1);
%         N = length(t);
%         
%         % load the noise files
%         for k = 1:length(noise_names)
%             
%             Nsegments = floor(length(noise_scaled{k})/N);
%             % pick a random segment
%             startIdx = randi(Nsegments-2 ,1)*N;
%             
%             noise_glob{k} = noise_scaled{k}(startIdx:startIdx+N -1);
%         end
%    
%         for n = 1:length(conditions)
%             
%             for k = 1:length(SNRs)
%                 noise = noise_glob{n};
%                 SpeechGain = 10^((GaindB + SNRs(k))/20);
%                  
%                 x = x_unscaled ./ sentenceFileLevel *SpeechGain;
%                 SpeechSPL_dBA(k) = LeqdBA(x,fs);
%                 
%                 if size(noise) ~= size(x)
%                     noise = noise';
%                 end
%                 test = noise + x;
%              
%                 tmp = joergensen2013(test,noise,fs,IOparameters);
%                 SNRenvs(k,n,q) = tmp.SNRenv;
%                 Pcorrect(k,n,q) = tmp.P_correct;
%                 
%             end
%             
%         end
%         amt_disp(['sentence nr: ' num2str(q) ' ' datestr(now, 'dd-mm-yyyy HH:MM:SS')]);
%         
%     end
%     
%     result.Pcorrect = Pcorrect;
%     result.SNRenvs = SNRenvs;
%     result.conditions =conditions;
%     result.SNRs =SNRs;
% 
%         
%     
%     %% Average across speech samples
%     Pc_est_mean = mean(result.Pcorrect,3);
%     
%     % ------------------  Estimating changes in SRTs based on the mean Pcorrect -------------
%     %   The first column of Pc_est_mean should always be the reference
%     selection = 1:length(conditions);
%     [dSRT SRTs] = joergensen2011_pctodsrt(Pc_est_mean,result.SNRs,selection);
%     
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


function [ret Signalscale]=Leq2dBA(sig,Fs,dBA)


InSig=sig/rms(sig).*10^(dBA/20);           % scale to dB SPL
x2=InSig;
[B,A] = adsgn(Fs);                    %
r=filter(B,A,InSig);


Laeq=20*log10(rms(r));       %dB(A)
Lleq=20*log10(rms(x2));      %dB
cor=Lleq-Laeq;
Signalscale=dBA+cor;
ret=InSig/rms(InSig).*10^(Signalscale/20);

function ret=LeqdBA(sig,Fs)

InSig=sig;           % scale to dB SPL

[B,A] = adsgn(Fs);
r=filter(B,A,InSig);
ret=20*log10(rms(r));       %dB(A)


function [B,A] = adsgn(Fs)
% ADSGN  Design of a A-weighting filter.
%    [B,A] = ADSGN(Fs) designs a digital A-weighting filter for
%    sampling frequency Fs. Usage: Y = FILTER(B,A,X).
%    Warning: Fs should normally be higher than 20 kHz. For example,
%    Fs = 48000 yields a class 1-compliant filter.
%
%    Requires the Signal Processing Toolbox.
%
%    See also ASPEC, CDSGN, CSPEC.

% Author: Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
%         couvreur@thor.fpms.ac.be
% Last modification: Aug. 20, 1997, 10:00am.

% References:
%    [1] IEC/CD 1672: Electroacoustics-Sound Level Meters, Nov. 1996.


% Definition of analog A-weighting filter according to IEC/CD 1672.
f1 = 20.598997;
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;
pi = 3.14159265358979;
NUMs = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]);
DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]);

% Use the bilinear transformation to get the digital filter.
[B,A] = bilinear(NUMs,DENs,Fs);


