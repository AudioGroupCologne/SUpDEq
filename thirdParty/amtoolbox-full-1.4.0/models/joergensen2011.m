function output = joergensen2011(x,y,fs_input,IO_param)
%JOERGENSEN2011  Speech-based envelope power spectrum (EPSM)
%   Usage: output = joergensen2011(x,y,fs_input,IO_param)
%
% 
%   Input parameters:
%     x          : noisy speech mixture 
%     y          : noise alone
%     fs         : sample rate in Hz
%     IO_param   : (optional) vector with parameters for the ideal observer 
%                  that converts the SNRenv to probability of correct, assuming a
%                  given speech material. It contains four parameters of the ideal observer 
%                  formatted as [k q m sigma_s].
%
%   Output parameters:
%     output  :     struct containing the SNRenv and the probability of correct given the
%                   SNRenv. This field is only included if IO_param is specified.
%                   Its calculation requires the Statistics ToolBox. 
%
%   output = JOERGENSEN2011(x,y,fs_input,IO_param) returns the output of
%   signal-to-noise envelope-power (SNRenv) ratio using the 
%   multi-resolution speech-based envelope spectrum model (mr-sEPSM)
%   described in Joergensen et al. (2013)
%
%   The model consists of the following stages:
%
%   1 A gammatone bandpass filterbank to simulate the auditory filters
%
%   2 An envelope extraction stage via the Hilbert Transform
%
%   3 A modulation filterbank
%  
%   4 Computation of the long-term envelope power (*output.SNRenv*)
%
%   5 A decision mechanism based on a statistically ideal observer (*output.P_correct*)
%
%   See also: joergensen2013 demo_joergensen2013 data_joergensen2011 sig_joergensen2011
%             plot_joergensen2011 joergensen2011_overlapadd3
%             joergensen2011_sim joergensen2011_pctodsrt joergensen2011_sepsub
%             joergensen2011_multchansnrenv joergensen2011_combineinformation
%             joergensen2013_sim exp_joergensen2011
%
%   References:
%     S. Joergensen and T. Dau. Predicting speech intelligibility based on
%     the signal-to-noise envelope power ratio after modulation-frequency
%     selective processing. J. Acoust. Soc. Am., 130(3):1475--1487, 2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/joergensen2011.php


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

if nargin < 3
    error('Too few input arguments.');
end

if length(x)~=length(y)
    error('x and y should have the same length');
end

% initialization
x           = x(:)';                             % Noisy speech row vector
y           = y(:)';                             % Noise alone speech row vector
fs          = 22050;
cf_aud = [63    80   100  125  160  200    250    315  400   500  630 800  1000  1250  1600  2000 2500 3150 4000 5000 6300  8000]; % centerfrequencies of the gammatone filters
cf_mod = [1 2 4 8 16 32 64 ];% centerfrequencies of the modulation filters
HT_diffuse = [37.5 31.5  26.5  22.1 17.9 14.4 11.4 8.4 5.8  3.8  2.1  1.0  0.8  1.9  0.5 -1.5 -3.1 -4.0 -3.8 -1.8 2.5   6.8 ]; % Diffuse field hearing threshold in quiet: ISO 389-7:2005


if fs_input ~= fs
    x	= resample(x, fs, fs_input);
    y 	= resample(y, fs, fs_input);
end

N = length(x);

%% Filtering through gammatone filterbank
g = gammatonefir(cf_aud,fs,'complex');
X  = 2*real(ufilterbank(x,g,1));
Y  = 2*real(ufilterbank(y,g,1));

%% ------------------ determining which frequency bands that are above
%% the hearing threshold.
% The spectrum levels (in SPL) of the stimuli are determined from a 1/3-octave analysis, since the threshold to compare with
% is spectrum levels.

mix_rms = thirdOctRMSAnalysis(x,fs,cf_aud);
mixRMS_dB =  20*log10(mix_rms);
bands = find(mixRMS_dB>HT_diffuse(1:length(cf_aud))); %     The bands to process further are the bands where the mix has spectrum levels above the hearing threshold.

%% ------------------ calculating envelopes of temporal outputs
%  lowpass filtering at 150 Hz
[bb, aa] = butter(1, 150*2/fs);

x_env = zeros(N,length(cf_aud));
y_env = zeros(N,length(cf_aud));

x_env(:,bands) = abs(hilbert(X(:,bands))); %  envelope
y_env(:,bands) = abs(hilbert(Y(:,bands)));

x_env = filter(bb,aa,x_env);
y_env = filter(bb,aa,y_env);

%downsampling
D = 10;
fsNew = fs/D;

x_env = resample(x_env,fsNew,fs);
y_env = resample(y_env,fsNew,fs);
%
SNRenv_n_p = zeros(length(cf_mod),length(cf_aud));

%% ----------------- Analysis via modulation filterbank
    P_env_xx = zeros(length(cf_mod),length(cf_aud));
    P_env_yy = P_env_xx;
    
for p = bands% For every audio channel
        
    
    P_env_xx(:,p)  = modFbank_v2(x_env(:,p),fsNew,cf_mod);
    P_env_yy(:,p)  = modFbank_v2(y_env(:,p),fsNew,cf_mod);
    
    
    if sum(sum(isnan(P_env_xx(:,p))))
        P_env_xx(isnan(P_env_xx(:,p)),n,p) = 0;
    end
    
    if sum(sum(isnan(P_env_yy(:,p))))
        P_env_yy(isnan(P_env_yy(:,p)),n,p) = 0;
    end
    
    P_env_yy(:,p) = min(P_env_xx(:,p),P_env_yy(:,p));    %     The envelope power of the noise is the minimum of Penv of the mixture
    %     or the noise.
   
    threshold = 0.01;          % The envelope power cannot go below 0.001 (-30 dB) reflecting
    % our minimum threshold of sensitivity to modulation detection
    
    P_env_yy(:,p) = max( P_env_yy(:,p),threshold);
    P_env_xx(:,p)= max(P_env_xx(:,p),threshold);
    
    SNRenvs_tmp = (P_env_xx(:,p)-P_env_yy(:,p)) ./P_env_yy(:,p); % calculation of SNRenv
    
    SNRenv_n_p(:,p) = max(0.001,SNRenvs_tmp); % Truncated at -30 dB for numerical reasons.
 
    
end

% ----------------- Integrating the SNRenv across audio and modulation bands---------------------------------------------------

SNRenv_p = sqrt(sum(SNRenv_n_p.^2,1));   %   Combine across modulation filters:  %SNRenv(p) eq. (4)

SNRenv = (sqrt(sum(SNRenv_p.^2)));%         Combine across audio filters: eq. (5)
output.SNRenv = SNRenv;

if nargin < 4
elseif nargin < 5
    P_correct = IdealObserver_v1(SNRenv,IO_param);
    output.P_correct = P_correct;
end
% ---------------------------------------------------------------------------------------------

end


function y = rms(x)
L = length(x);
y=norm(x)/sqrt(L);
end


function rms_out = thirdOctRMSAnalysis(x,fs,midfreq)

if nargin<3
    midfreq=[63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 ];
end

N = length(x);
X = (fft(x));
X_mag  = abs(X) ;
X_power = X_mag.^2/N ;% power spectrum.
X_power_pos = X_power(1:fix(N/2)+1) ;
X_power_pos(2:end) = X_power_pos(2:end).* (2)  ; %take positive frequencies only and mulitply by two-squared to get the same total energy(used since the integration is only performed for positive freqiencies)

freq= linspace(0,fs/2,length(X_power_pos));

%resolution of data
resol=freq(2)-freq(1);

crossfreq(1)=midfreq(1)/(2^(1/6));
crossfreq(2:length(midfreq)+1)=midfreq*(2^(1/6));

%cross-over indicies
y=crossfreq/resol;

%rounding up
crosselem=round(y);
for n=1:length(y)
    if crosselem(n)<y(n)
        crosselem(n)=crosselem(n)+1;
    end
end

nn=1;
rms_out(1:length(midfreq)) = 0;
while crossfreq(nn+1)<=freq(end)% for nn =1:length(crossfreq)-1
 
    rms_out(nn) = sqrt(sum(X_power_pos(crosselem(nn):crosselem(nn+1)-1))/N);
    nn=nn+1;
    if 1+nn > length(crossfreq)
        break
    end
end

end


function P_env = modFbank_v2(Env,fs,cf_mod)
%
% This function is an implementation of a modulation filterbank similar to the EPSM-filterbank
% as presented by Ewert & Dau 2000. This implementation consists of a lowpass
% filter with a cutoff at 1 Hz, in parallel with 8 bandpass filters with
% octave spacing. the Center-frequencies of the bandpass filters are lower
% than the original from Ewert & Dau (2000).
%
% Inputs:
%   Env:  The envelope to be filtered
%   fs: sampling frequency of the envelope
%   cf_mod:  centerfrequencies of the modulation filters
%
% Outputs:
%   P_env:  Envelope power for each of the modulation filters
%
%
% Created by Søren Jørgensen jan 2010
% last update 02 may 2014
% Copyright Søren Jørgensen

if nargin<3
    %band center frequencies
    cf_mod=[1 2 4 8 16 32 64];
end

if size(Env,1) > 1
    Env = Env';
end

if mod(length( Env),2) == 0
    %number is even
    Env =  Env(1:end-1);
else
    %number is odd
end

Q = 1;
N = length(Env);
X = fft(Env);
X_mag  = abs(X) ;
X_power = X_mag.^2/N ;% power spectrum.

X_power_pos = X_power(1:fix(N/2)+1) ;
X_power_pos(2:end) = X_power_pos(2:end).* (2)  ; %take positive frequencies only and mulitply by two-squared to get the same total energy(used since the integration is only performed for positive freqiencies)
pos_freqs= linspace(0,fs/2,length(X_power_pos));
freqs = [pos_freqs -1*fliplr(pos_freqs(2:end))];

% Initialize transfer function
TFs = zeros(length(cf_mod),length(freqs));

% Calculating frequency-domain transferfunction for each center frequency:
for k = 2:length(cf_mod)
    TFs(k,1:end) = 1./(1+ (1j*Q*(freqs(1:end)./cf_mod(k) - cf_mod(k)./freqs(1:end)))); % p287 Hambley.
end
% squared filter magnitude transfer functions
Wcf = (abs(TFs)).^2;

fcut = 1;% cutoff frequency of lowpassfilter:
n = 3;% order:
% Lowpass filter squared transfer function:
Wcf(1,:) =  1./(1+((2*pi*freqs/(2*pi*fcut)).^(2*n))); % third order butterworth filter TF from: http://en.wikipedia.org/wiki/Butterworth_filter

% ------------ DC-power, --------------------------
% %here devided by two such that a fully modulated tone has an AC-power of 1.
DC_power =  X_power_pos(1)/ N /2;

X_power_pos = repmat(X_power_pos,length(cf_mod),1);

Vout = X_power_pos.*Wcf(:,1:floor(end/2)+1);
% Integration estimated as a sum from f > 0
% integrate envelope power in the
% passband of the filter. Index goes from 2:end since
P_env =  sum(Vout(:,2:end),2)/ N / DC_power;
% integration is for f>0


end


function Pcorrect  = IdealObserver_v1(SNRenv_lin,parameters)
%%
% IdealObserver: Converts the overall SNRenv to percent correct.
%
% Usage: [Pcorrect SNRenv ] = IdealObserver(SNRenv_lin,parameters)
%
% SNRenv_lin :  vector with the SNRenv values (not in dB) for each input SNR
% Parameters :  vector with the parameters for the ideal Observer formatted as [k q m sigma_s]
%
% Søren Jørgensen august 2010
% last update 10-June 2013
% Copyright Søren Jørgensen
%%

if nargin < 2
    error('You have to specify the k,q,m,sigma_s parameters for the IdealObserver')
end
k = parameters(1);
q = parameters(2);
m = parameters(3);
sigma_s = parameters(4);


% ---------- Converting from SNRenv to d_prime  --------------
d_prime = k*(SNRenv_lin).^q;

%----------- Converting from d_prime to Percent correct, Green and Birdsall (1964)----------
Un = 1*norminv(1-(1/m));
mn = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
sig_n=  1.28255/Un;
Pcorrect = normcdf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;


% Green, D. M. and Birdsall, T. G. (1964). "The effect of vocabulary size",
% In Signal Detection and Recognition by Human Observers,
% edited by John A. Swets (John Wiley & Sons, New York)
end


