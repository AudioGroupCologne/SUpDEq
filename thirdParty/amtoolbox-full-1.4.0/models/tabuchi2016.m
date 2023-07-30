function [Kfunc_vec, Kval] = tabuchi2016(insig, nPhase, data, varargin)
%TABUCHI2016 calculates difference of predicted masker thresholds
%   Usage: [kfunc_acrCvalGmax, kvals] = tabuchi2016(insig, nPhase, dataOut)
%
%   Input parameters:
%     insig :            input signal
%     nPhase :           phase index
%     data :             input data
%
%   Output parameters:
%     kfunc_acrCvalGmax :    matrix (target level, Cvec, Gmax, precursor condition, masker level)
%     kvals :                values of the k function
%
%
%    TABUCHI2016 takes a a roex-weighted Schroeder phase complex as an input.
%    To simplify the code, the frequency weighting is already saved in the 
%    data file (e.g., Attenu_Roex_GrandMeanSd_OffFreq_60dB.mat). Then, the
%    I/O function processing is performed to get the RMS Outputs of Masker
%    and Masker + Target from C-1.5 to 1. Then, the average of K based on 
%    Gmax 34 and C from are retrieved, which involves the processing of 
%    the I/O function.
%
%   Optional parameters:
%
%     'C',C                      C (curvature) value
%
%     'Gmax',Gmax                maximum gain value
%
%     'mlvl',mlvl                sound pressure level
%
%     'GmaxToGetK',GmaxToGetK    conversion parameter to get K value
%
%     'gamma',gamma              gamma value
%
%     'beta',beta                beta value
%
%     'freq',freq                masker frequency [Hz]
%
%     See also: exp_tabuchi2016 sig_tabuchi2016
%
%   References:
%     H. Tabuchi, B. Laback, T. Necciari, and P. Majdak. The role of
%     compression in the simultaneous masker phase effect. The Journal of the
%     Acoustical Society of America, 140(4), 2016.
%     
%     B. Glasberg and B. Moore. Frequency selectivity as a function of level
%     and frequency measured with uniformly exciting notched noise. The
%     Journal of the Acoustical Society of America, 108:2318--28, 12 2000.
%     
%     N. P. Cooper. Harmonic distortion on the basilar membrane in the basal
%     turn of the guinea‐pig cochlea. The Journal of Physiology, 509, 1998.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/tabuchi2016.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Verified
%   #Requirements: 
%   #Author: Hisaaki Tabuchi
%   #Author: Clara Hollomey (adaptations for AMT)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.


% % initial parameters
definput.keyvals.C = 0; %-1.5:0.25:1; % This range of Cvec was tested by considering the presumed curvature, -0.5. See the footnote 2 in Tabuchi et al. (2016). 
definput.keyvals.Gmax = 60; %0:70;
definput.keyvals.mlvl = 60; %[60 90];
definput.keyvals.GmaxToGetK = 34;
definput.keyvals.gamma = 60; % for IHC 
definput.keyvals.beta = 1; % for IHC
definput.keyvals.freq = 4; %masker frequency [Hz]
% 
[~,kv]  = ltfatarghelper({},definput,varargin);


 CvecLen_flag = 0; % to sort out the if sentence below
% 
    Slvl_in = 0:0.1:110;
    Model_dBOut_MplsT_vec = NaN(length(Slvl_in),1);
    
    %model outputs for masker plus target
    for Slvl_in_loop = 1:length(Slvl_in)

        [model_wavout] = local_outdiffmaskertargetcore(Slvl_in(Slvl_in_loop),...
            insig, kv.freq, nPhase, kv.Gmax, kv.beta, kv.gamma);

        Model_dBOut_MplsT_vec(Slvl_in_loop,1) = dbspl(local_rmsamp(model_wavout),...
            'dboffset', 93.9794);

    end
    %model outputs for masker
    [model_wavout] = local_outdiffmaskercore(insig, kv.Gmax, kv.beta, kv.gamma);
    Model_dBOut_M_vec = ones(length(Slvl_in),1)*dbspl(local_rmsamp(model_wavout),...
        'dboffset', 93.9794);


    Kfunc_vec = Model_dBOut_MplsT_vec - Model_dBOut_M_vec; 



    if (kv.Gmax == kv.GmaxToGetK) && (kv.C >= -1)
    % get the K value based on data (Cvec = -1:0.25:1) and GmaxToGetK = 34.

        CvecLen_flag = CvecLen_flag + 1; % to sort out the outer loop
        Slvl = data(CvecLen_flag,2); % get a target threshold

        %Combine a masker and target, extract the envelope, and call the I/O function
        [model_wavout_k_mt] = local_outdiffmaskertargetcore(Slvl,...
            insig, kv.freq, nPhase, kv.Gmax, kv.beta, kv.gamma);                 
        %Extract the envelope, and call the I/O function
        [model_wavout_k_m] = local_outdiffmaskercore(insig,...
            kv.Gmax, kv.beta, kv.gamma); 

        Kval = dbspl(local_rmsamp(model_wavout_k_mt), 'dboffset', 93.9794) -...
            dbspl(local_rmsamp(model_wavout_k_m), 'dboffset', 93.9794);

    else
        Kval = [];
    end

end




function [model_wavout] = local_outdiffmaskertargetcore(spl, WaveMout, frq, nPhase, Gmax, beta, gamma)
% Purpose: Combine a masker and target, extract the envelope, and call the I/O function                    
% Command line:
%        [model_wavout] = OutDiffFunc_MplsT_Env_CoreFunc(Slvl_in, WaveMout, nPhaseOut, Gmax)
% Example inputs:
% Slvl_in = 0:0.1:110; % 0.1 is the step required to get threshold.
% WaveMout; % a masker waveform
% nPhaseOut % the vector of masker's phase
% Gmax = 70;
%

%(sfs, frq, spl, leadms, trailms, durms, riz, nPhase)
sfs = 48;
%frq = 4;
leadms = 5;
trailms = 5;
durms = 30;
%%% generate a target tone
%[wavOut] = local_tonephase(48, 4, spl, 5, 5, 30, 0, nPhase); % generate a target tone
% synthesize
smpl = 1/sfs;
tt = (smpl:smpl:durms)'; % waveform should be a column vector %

wav1 = cos( 2*pi*frq*tt + nPhase(round(length(nPhase)/2),1)); % a sinusoid starts with the phase in the center component of masker %      

% leading zeros
nz1 = round(leadms/smpl);
wav1 = cat(1,zeros(nz1,1),wav1);

% trailing zeros
nz2 = round(trailms/smpl);
wav1 = cat(1,wav1,zeros(nz2,1));

% scale amplitude 
wavOut = local_amp2press(wav1,spl);

%%% add the target and masker in phase
[wavMix] = local_mixer(wavOut,WaveMout);

%%% get the hilbert envelope
xhi = hilbert(wavMix);
xEnv = sqrt(xhi.*conj(xhi)); 

%%% the level of Envelope decreased by 3.0103 dB SPL is equivalent to the level of input waveform.
xEnvIn = xEnv*10^(-3.0103/20);

%%% get the model response from Glasberg's I/O
[model_wavout_glasberg] = local_iogmaxpascal(xEnvIn, Gmax);

%%% get the model response from Cooper's Boltzmann function
alpha = 100; % the maximum y of function
% beta = 1; % a default value of Boltzmann function
% gamma = 20; % slope in dB. the smaller gamma value, the steeper slope 
Xshift = -20; % the horizontal shift of function
[model_wavout] = local_ioboltzmannpascal(model_wavout_glasberg,...
    alpha, beta, gamma, Xshift);
end


function [model_wavout] = local_outdiffmaskercore(WaveMout, Gmax, beta, gamma)
% Purpose: Extract the envelope, and call the I/O function                    
% Command line:
%        [model_wavout] = OutDiffFunc_M_Env_CoreFunc(WaveMout, Gmax)
% Example inputs:
% WaveMout; % a masker waveform
% Gmax = 70;
%


%%% get the hilbert envelope
xhi = hilbert(WaveMout);
xEnv = sqrt(xhi.*conj(xhi)); 

%%% the level of Envelope decreased by 3.0103 dB SPL is equivalent to the level of input waveform.
xEnvIn = xEnv*10^(-3.0103/20);
% xEnvIn_dB = dbspl(RMSamp(xEnvIn));

%%% get the model response from Glasberg's I/O
[model_wavout_glasberg] = local_iogmaxpascal(xEnvIn, Gmax);

%%% get the model response from Cooper's Boltzmann function
alpha = 100; % the maximum y of function
% beta = 1; % a default value of Boltzmann function
% gamma = 20; % slope in dB. the smaller gamma value, the steeper slope 
Xshift = -20; % the horizontal shift of function
[model_wavout] = local_ioboltzmannpascal(model_wavout_glasberg,...
    alpha, beta, gamma, Xshift);
end



function [wavout] = local_iogmaxpascal(wavin_pascals, gmax)
% IOfunc_InstLvl_Glasberg2000 - get a time waveform after processed by 
% the input-output function (Glasberg et al, 2000)
% command line:
%        [wavout] = IOfunc_InstLvl_Gmax_Pascal(wavin, dbval, gmax)
% example inputs:
% wavin_pascals % the vector of a time waveform with pascals
% gmax = 34.6; % The estimated value at 4 kHz (Glasberg et al, 2000).
%


ref = 2*10^-5;

wav_db = 20*log10(abs(wavin_pascals)/ref);

vecpos = ones(length(wav_db), 1);

ind = find(wav_db <= 0); % the small instantaneous dB (<0) is not processed by I/O function

if isempty(ind) ~= 1  % replace the small dB (<0) with NaN
    for ind_n = 1:length(ind)
        vecpos(ind(ind_n)) = NaN;
    end
end

%try
wav_db_in = wav_db.*vecpos;
 %   dbstop if error
%catch
%    size(wav_db)
%    size(vecpos)

%end

a = -0.0894*gmax + 10.894;

b = 1.1789*gmax - 11.789;

dBspl_vec = 0.9*wav_db_in + a + b * (1 - (1./(1 + exp(-0.05*(wav_db_in-50))))); % dB SPL at a time point is processed by I/O

io_out = ref * 10.^(dBspl_vec/20); % convert dB SPL to pascals

wavout = sign(wavin_pascals).*io_out; % undo the full-wave rectification

if isempty(ind) ~= 1 % replace NaN with the input instantaneous value
    for ind_n = 1:length(ind)
        wavout(ind(ind_n)) = wavin_pascals(ind(ind_n));
    end
end

end


function [wavout] = local_ioboltzmannpascal(wavin_pascals, alpha, beta, gamma, Xshift)
% IOfunc_InstLvl_Boltzmann_Pascal - get a time waveform after processed by 
% the Boltzmann function (Cooper, 1998)
% command line:
%        [wavout] = IOfunc_InstLvl_Boltzmann_Pascal(wavin_pascals, alpha, beta, gamma, Xshift)
% example inputs:
% wavin_pascals % the vector of a time waveform with pascals
% alpha = 100; % the maximum y of function
% beta = 1; % a default value of Boltzmann function
% gamma = 20; % slope in dB. the smaller gamma value, the steeper slope 
% Xshift = -20; % the horizontal shift of function
%


ref = 2*10^-5;

wav_db = 20*log10(abs(wavin_pascals)/ref);

vecpos = ones(length(wav_db),1);

ind = find(wav_db <= 0); % the small instantaneous dB (<0) is not processed by I/O function

if isempty(ind) ~= 1  % replace the small dB (<0) with NaN
    for ind_n = 1:length(ind)
        vecpos(ind(ind_n)) = NaN;
    end
end

wav_db_in = wav_db.*vecpos;

[dBspl_vec] = local_boltzmann1storder(wav_db_in, alpha, beta, gamma, Xshift); % Boltzmann function

io_out = ref * 10.^(dBspl_vec/20); % convert dB SPL to pascals

wavout = sign(wavin_pascals).*io_out; % undo the full-wave rectification

if isempty(ind) ~= 1 % replace NaN with the input instantaneous value
    for ind_n = 1:length(ind)
        wavout(ind(ind_n)) = wavin_pascals(ind(ind_n));
    end
end

end


function [out] = local_boltzmann1storder(xvec, alpha, beta, gamma, Xshift)
% Purpose:
% - Return the output value of 1st Boltzmann function in dB (Cooper, 1998).
%
% Input examples:
% xvec = 0:110; % the value of x-axis in dB
% alpha = 100; % the maximum y of function
% beta = 1; % a default value of Boltzmann function
% gamma = 20; % slope in dB. the smaller gamma value, the steeper slope 
% Xshift = -20; % the horizontal shift of function
% 
% Notes:
% - Beta determines a slope and a saturating point of function.
% - Gammma also determines a slope of function. 
% - The term, (alpha/(alpha-Yshift)) gurantees a maximum of the y value with alpha. 
% - The function goes through the origin (x,y)=(0,0) when Xshift = 0.
%


term = beta * exp((gamma^-1)*(Xshift - xvec)); 
Yshift = (alpha * (1 + beta)^-1);
out = (alpha/(alpha-Yshift))*((alpha * (1 + term).^-1) - Yshift);

end



function [wavM]=local_mixer( wavS,wavN )
% Mixer - check length, then sum two waveforms
% command line:
%        [wavM]=Mixer( wavS,wavN )

% if lengths differ, truncate the longest one
if length(wavS)>length(wavN)
    nnn=length(wavS)-length(wavN);
%     fprintf(1,'%i training zeros added to WavN...\n',nnn);
    wavN=cat(1,wavN,zeros(nnn,1));
elseif length(wavS)<length(wavN)
    nnn=length(wavS);
%     fprintf(1,'WavN truncated from %i to %i samples...\n',length(wavN),nnn);
    wavN=wavN(1:nnn);
end
% mix
wavM=wavS+wavN;
end


function [ww2,sfact]=local_amp2press( www,SPL )
% =======================================
% arbitrary amplitude to Pascals
% command line:
%        [ww2,sfact]=amp2press( www,SPL )
% constant

pref=0.00002;

% check amplitude, excluding leading/trailing silence
mxx=max(abs(www));
n1=find( www>(0.02*mxx),1,'first' );
n2=find( www>(0.02*mxx),1,'last' );
rrr=local_rmsamp(www(n1:n2));
%ss=sum(www(n1:n2).*www(n1:n2));
%rrr = sqrt( ss/length(www(n1:n2)) );
if rrr>0
    dbtmp=20*log10( rrr/pref );
    dbdif=SPL-dbtmp;
    sfact=10^(dbdif/20);
    ww2=sfact*www;
else
    ww2=www;
    sfact=1;
end

end


function [rrr] = local_rmsamp( dd )
% ===================================
% calculate rms amplitude of a matrix
% command line:
%        [rrr] = RMSamp( dd )
% rrr is a scalar if dd is a vector,
% rrr is a vector of length nc if rrr is an nr x nc  matrix
%   (e.g., operates on columns)
%
ss=sum(dd.*dd);
rrr = sqrt( ss/length(dd) );

end

