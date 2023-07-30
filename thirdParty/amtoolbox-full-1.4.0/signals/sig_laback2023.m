function [stimTarg, stimPrec, noiseILDvec, TotDur, B_P, A_P, PrecNoise] = sig_laback2023(precmode,tvec, nvec, ILDNoiseLin, PrecITD_side, PrecITD_t, model, varargin)
%SIG_LABACK2023 generates stimuli used in exp_laback2023
%
%   Usage: [stimTarg, stimPrec, noiseILDvec, TotDur] = sig_laback2023(precmode, frozennoise, Experiment, tvec, nvec, CF, Fs, F0cue, RampDur, stimdb, DurT, DurN, DurG, factor, ILDNoiseLin, PrecITD_side, PrecITD_t)
%          [stimTarg, stimPrec, noiseILDvec, TotDur, B_P, A_P, PrecNoise] = sig_laback2023(precmode, frozennoise, Experiment, tvec, nvec, CF, Fs, F0cue, RampDur, stimdb, DurT, DurN, DurG, factor, ILDNoiseLin, PrecITD_side, PrecITD_t)
%
%   Input parameters:
%     precmode     : if there is a precursor or not, can be 0 or 1
%     frozennoise  : if noise is constant or not, can be 0 or 1
%     Experiment   : experiment index
%     tvec         : target vector
%     nvec         : noise vector
%     ILDNoiseLin  : linear noise ILD
%     PrecITD_side : precursor ITD side
%     PrecITD_t    : precursor ITD target
%
%
%   Output parameters:
%     stimTarg     : target stimulus
%     stimPrec     : precursor stimulus
%     noiseILDvec  : noise ILD vector
%     TotDur       : total duration
%
%   SIG_LABACK2023 generates a noise, shaped by a raised cosine, with or
%   without a precursor.
%
%   Optional parametersL
%     CF           center frequency [Hz]
%     Fs           sampling frequency [Hz]
%     F0cue        if 0: T and P have same envelope rate, if 1: T has double envelope rate as P
%     RampDur      duration of onset ramp
%     stimdb       stimulus level [dB SPL]
%     DurT         duration of target
%     DurN         duration of noise
%     DurG         gap duration
%     factor       amplification factor for transposed stimuli
%     ILDNoiseLin  linear noise ILD
%
%   See also: exp_laback2023 data_laback2023
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_laback2023.php


% #Author: Bernhard Laback (2023)
% #Author: Alejandro Osses (2023) : Integration in AMT 1.4

definput.keyvals.SubExp = [];
definput.import = {'laback2023'};
if model == 1
    definput.importdefaults={'exp1'};
else
    definput.importdefaults={'exp2'};
end

[flags,kv]  = ltfatarghelper({},definput,varargin);

PrecNoise = [];

fs = kv.fs;
[B_P, A_P, B_T, A_T, C, D] = local_getfiltercoefs(tvec, nvec, fs, kv.Experiment, kv.F0cue, kv.SubExp, kv.CF);

if precmode == 0

    if kv.frozennoise == 0 %random noise sample for each trial
        TargNoise=randn(length(tvec),1);
    end
    if kv.Experiment == 1
        stimTarg = filter(B_T, A_T, TargNoise);
    elseif kv.Experiment >= 2                        
        LFband=filter(B_T, A_T, TargNoise);
        index=find(LFband >= 0);
        LFband_env=zeros(1, length(tvec));
        LFband_env(index)=LFband(index);
        LFband_filt=filter(D,C,LFband_env);
        TargNoiseTrans = LFband_filt .* sin(2*pi*tvec*kv.CF);
        stimTarg = TargNoiseTrans'*kv.factor;
    end
    stimTarg = local_raisedcosine(stimTarg,fs,kv.RampDur);
    stimTarg= stimTarg * sqrt(2)* 20e-6*10^(kv.stimdb/20) / rms(stimTarg); 
    stimPrec = [];
    noiseILDvec = [];

    TotDur=kv.DurT; %Duration of stimulus
elseif  precmode == 1 
    if kv.frozennoise == 0 %random noise sample for each trial
        %Generate noise vectors
        TargNoise=randn(length(tvec),1);
        PrecNoise=randn(length(nvec),1);
    end

    if kv.Experiment == 1
        stimTarg = filter(B_T, A_T, TargNoise);
        stimPrec = filter(B_P, A_P, PrecNoise);
    elseif kv.Experiment >= 2   
        %Target
        LFband=filter(B_T, A_T, TargNoise);
        index=find(LFband >= 0);
        LFband_env=zeros(1, length(tvec));
        LFband_env(index)=LFband(index);
        LFband_filt=filter(D,C,LFband_env);
        TargNoiseTrans = LFband_filt .* sin(2*pi*tvec*kv.CF);
        stimTarg = TargNoiseTrans'*kv.factor;

        %Precursor
        LFband=filter(B_P, A_P, PrecNoise);
        index=find(LFband >= 0);
        LFband_env=zeros(1, length(nvec));
        LFband_env(index)=LFband(index);
        LFband_filt=filter(D,C,LFband_env);
        PrecNoiseTrans = LFband_filt .* sin(2*pi*nvec*kv.CF);
        stimPrec = PrecNoiseTrans'*kv.factor;
    end


    stimTarg = local_raisedcosine(stimTarg,fs,kv.RampDur);
    stimTarg= stimTarg * sqrt(2)* 20e-6*10^(kv.stimdb/20) / rms(stimTarg);
    targetlevel=rms(stimTarg);    

    stimPrec = local_raisedcosine(stimPrec,fs,kv.RampDur);
    preclevel=rms(stimPrec);

    %Add ILD to precursor
    stimPrec = stimPrec * targetlevel/preclevel;

    noiseILDvec_(:,1) = stimPrec * sqrt(ILDNoiseLin);
    noiseILDvec_(:,2) = stimPrec / sqrt(ILDNoiseLin);

    %Apply ITD to precursor
    noiseILDvec = local_ITD(noiseILDvec_,PrecITD_side,PrecITD_t, fs);

    TotDur=kv.DurT+kv.DurN+kv.DurG; %Duration of stimulus
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B_P, A_P, B_T, A_T, C, D] = local_getfiltercoefs(tvec, nvec, Fs, Experiment, F0cue, SubExp, CF)

%Define time vectors for P and T
%nvec = 0: 1/Fs : DurN-1/Fs; 
%tvec = 0: 1/Fs : DurT-1/Fs;
C = [];
D = [];
%Generate noises for P and T
PrecNoise=randn(length(nvec),1);
TargNoise=randn(length(tvec),1);

if Experiment == 1
    %P filter
    fl = CF/sqrt(2); 
    fu = fl*(2^1.2); 
    [B_P,A_P]=butter(4,[ fl/(Fs/2) fu/(Fs/2) ]);

    %T filter
    fu = CF*sqrt(2);
    [B_T,A_T]=butter(4,[ fl/(Fs/2) fu/(Fs/2) ]);

elseif Experiment >= 2
    if SubExp == 0 %Exp 2A
        %P filter        
        fu = 138; %CF: 125 Hz; BW: 25 Hz;  
        fl = 113; 
        [B_P,A_P]=butter(3,[ fl/(Fs/2) fu/(Fs/2) ]);
        [D,C]=butter(12,2000/(Fs/2) , 'low'); %low-pass filter for filtering envelope signal 
        %T filter
        if F0cue == 0
            B_T = B_P;
            A_T = A_P;
        elseif F0cue == 1
            fu = 263; %CF: 250 Hz; BW: 25 Hz; 
            fl = 238; 
            [B_T,A_T]=butter(3,[ fl/(Fs/2) fu/(Fs/2) ]);
        end  
    elseif SubExp == 1 %Exp 2B
        %P filter        
        fu = 263;  %CF: 250 Hz; BW: 25 Hz;   
        fl = 238; 
        [B_P,A_P]=butter(3,[ fl/(Fs/2) fu/(Fs/2) ]);
        [D,C]=butter(12,2000/(Fs/2) , 'low'); %low-pass filter for filtering envelope signal 
        %T filter
        if F0cue == 0 %Target has same (high) envelope  rate
            B_T = B_P;
            A_T = A_P;
        elseif F0cue == 1 %Target has different (lower) envelope  rate
            fu = 138; %CF: 125 Hz; BW: 25 Hz; 
            fl = 113; 
            [B_T,A_T]=butter(3,[ fl/(Fs/2) fu/(Fs/2) ]);
        end  
    end 
end

% function [stimTarg, stimPrec, noiseILDvec, TotDur] = local_signal(precmode,...
%     frozennoise, Experiment, tvec, nvec, CF, Fs, B_T, A_T, B_P, A_P, D, C, RampDur, stimdb, DurT, DurN, DurG, factor, ILDNoiseLin, PrecITD_side, PrecITD_t)
% 
%     if precmode == 0
% 
%         if frozennoise == 0 %random noise sample for each trial
%             TargNoise=randn(length(tvec),1);
%         end
%         if Experiment == 1
%             stimTarg = filter(B_T, A_T, TargNoise);
%         elseif Experiment >= 2                        
%             LFband=filter(B_T, A_T, TargNoise);
%             index=find(LFband >= 0);
%             LFband_env=zeros(1, length(tvec));
%             LFband_env(index)=LFband(index);
%             LFband_filt=filter(D,C,LFband_env);
%             TargNoiseTrans = LFband_filt .* sin(2*pi*tvec*CF);
%             stimTarg = TargNoiseTrans'*factor;
%         end
%         stimTarg = local_raisedcosine(stimTarg,Fs,RampDur);
%         stimTarg= stimTarg * sqrt(2)* 20e-6*10^(stimdb/20) / rms(stimTarg); 
%         stimPrec = [];
%         noiseILDvec = [];
%         
%         TotDur=DurT; %Duration of stimulus
%     elseif  precmode == 1 
%         if frozennoise == 0 %random noise sample for each trial
%             %Generate noise vectors
%             TargNoise=randn(length(tvec),1);
%             PrecNoise=randn(length(nvec),1);
%         end
% 
%         if Experiment == 1
%             stimTarg = filter(B_T, A_T, TargNoise);
%             stimPrec = filter(B_P, A_P, PrecNoise);
%         elseif Experiment >= 2   
%             %Target
%             LFband=filter(B_T, A_T, TargNoise);
%             index=find(LFband >= 0);
%             LFband_env=zeros(1, length(tvec));
%             LFband_env(index)=LFband(index);
%             LFband_filt=filter(D,C,LFband_env);
%             TargNoiseTrans = LFband_filt .* sin(2*pi*tvec*CF);
%             stimTarg = TargNoiseTrans'*factor;
% 
%             %Precursor
%             LFband=filter(B_P, A_P, PrecNoise);
%             index=find(LFband >= 0);
%             LFband_env=zeros(1, length(nvec));
%             LFband_env(index)=LFband(index);
%             LFband_filt=filter(D,C,LFband_env);
%             PrecNoiseTrans = LFband_filt .* sin(2*pi*nvec*CF);
%             stimPrec = PrecNoiseTrans'*factor;
%         end
% 
% 
%         stimTarg = local_raisedcosine(stimTarg,Fs,RampDur);
%         stimTarg= stimTarg * sqrt(2)* 20e-6*10^(stimdb/20) / rms(stimTarg);
%         targetlevel=rms(stimTarg);    
% 
%         stimPrec = local_raisedcosine(stimPrec,Fs,RampDur);
%         preclevel=rms(stimPrec);
% 
%         %Add ILD to precursor
%         stimPrec = stimPrec * targetlevel/preclevel;
% 
%         noiseILDvec_(1,:) = stimPrec * sqrt(ILDNoiseLin);
%         noiseILDvec_(2,:) = stimPrec / sqrt(ILDNoiseLin);
% 
%         %Apply ITD to precursor
%         noiseILDvec = local_ITD(noiseILDvec_,PrecITD_side,PrecITD_t, Fs);
% 
%         TotDur=DurT+DurN+DurG; %Duration of stimulus
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal = local_raisedcosine(signal,fs,timeS)
%signal = local_raisedcosine(signal,fs,timeS)
%applies raised-cosine Onset and Offset ramp
%Inputs:
%signal...input signal mono
%fs...sampling rate [Hz]
%timeS is ramp time [s]
%AUTHOR: Bernhard Laback, 2021

size = length(signal);

time = timeS/1000;
rs = round(time*fs); 
for i=1:(rs-1)
   ampfac = 0.5 - 0.5*cos(pi*i/rs);
   signal(i) = ampfac*signal(i);
end

for i=(size-rs+1):size 
   ampfac = 0.5 - 0.5*cos(pi*(size-i)/rs);
    signal(i) = ampfac*signal(i);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = local_ITD(insig,side,ITD, fs)
%Function to apply ITD to input signal 
%B. Laback, 03.06.2016  
%Call: [out] = ITD(in,side,ITD, fs)
%Paramters
%in         binaural signal (2-D vector) to which ITD shall be applied
%side       -1: leading left; 0:zero ITD (zeros added at end of signal at both ears; +1: leading right 
%ITD        ITD in us
%fs         sampling rate (in Hz)
%out        output signal containing ITD

inputL = insig(:,1);
inputR = insig(:,2);

ITDsamples = round(ITD/(1000000/fs));
Lag = zeros(ITDsamples,1);

%Introduce ITD
if side < 0
    outsig(:,1) = [inputL; Lag];
    outsig(:,2) = [Lag; inputR];
elseif side == 0
    outsig(:,1)=[inputL; Lag];
    outsig(:,2)=[inputR; Lag];
elseif side > 0
    outsig(:,1) = [Lag; inputL];
    outsig(:,2) = [inputR; Lag];
end

