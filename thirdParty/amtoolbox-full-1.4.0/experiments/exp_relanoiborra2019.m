function varargout = exp_relanoiborra2019(varargin)
%EXP_RELANOIBORRA2019 experiments from Relano-Iborra et al. 2019
%
%   Usage: exp_relanoiborra2019('fig3');
%
%   This script reproduces Figure 3 of Relano-Iborra et al. It finds the free 
%   parameters of the logistic function that relates the output of the sCASP 
%   model (correlation coefficient d) with the percentage of correct answers. 
%   The data used for this fitting is taken from Nielsen and Dau (2009). The 
%   testing conditions are replicated and fed to the sCASP model. Later a least 
%   squares analysis is used to fit the function to the data.
%
%   Examples:
%   ---------
%
%   To display Figure 3 use :
%
%     exp_relanoiborra2019('fig3');
%
%   See also: relanoiborra2019
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_relanoiborra2019.php


%   #Author: Helia Relano Iborra (2020): Code provided for the AMT
%   #Author: Piotr Majdak (2021): Various bug fixes for the AMT 1.0
%   #Author: Clara Hollomey (2021): Adaptations for the AMT 1.0

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.import = {'amt_cache'};

definput.flags.type = {'missingflag', 'fig3'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},...
             definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.', ...
      upper(mfilename),flagnames);
end

%% Initialization

Pcorrect_human = [0 8 35 71 90 100 ]; %% Human data
SNRs= -8:2:2;

%load single_150_SentArray22kHz_varLength % CLUE speech
x = amt_load('relanoiborra2019', 'single_150_SentArray22kHz_varLength.mat');
sentenceArray = x.sentenceArray;
fsSent = 22050;

noise_name = 'SSN_CLUE_22kHz.wav'; % SSN Noise
speechSPL = 65;
Nsentences=25;

fsRef = 22050; % sampling freq.

Pref= 20e-6; % Transformation to Pascals

d= zeros(length(SNRs), Nsentences);

%% Run experiment:

for q=1:Nsentences
    
    amt_disp(['Processing sentence: ' num2str(q) ' out of ' num2str(Nsentences)]);
    speech  = sentenceArray{q};
    speech = resample(speech, fsRef, fsSent );
    
    N_org= length(speech);  % Calculate length of original sentence
    speech = [speech; speech]; % Prepane the same sentence
    
    speech = Pref*speech*(1/rms(speech))*10^((speechSPL)/20); % Set speech level
    N = length(speech);  % Overall speech size
    
    
    for n=1:length(SNRs)
        
        %noise = audioread(noise_name);
        noise = amt_load('relanoiborra2019', noise_name);
        Nsegments = floor(length(noise)/N);
        startIdx = randi(Nsegments-2 ,1)*N;
        noise = noise(startIdx:startIdx+N -1)'; % random segment from the noise file
        
        noise = Pref*(noise./rms(noise)*10^((speechSPL-SNRs(n))/20)); % sets the level of the noise signal
        
        if size(noise) ~= size(speech)
            noise = noise'; % Flips noise signal if needed
        end
        
        test = noise + speech; % Gerating the noisy mixture
        
        % DRNL conf:
        flow = 100;
        fhigh = 8e3;
        
        sbj = 'NH';
        
        % Call model:
        tmp = relanoiborra2019(speech, test, fsRef, flow, fhigh, 'N_org',N_org, 'subject',sbj);
        
        d(n, q) = tmp.dfinal;  % correlation value per sentence and SNR
        
    end % End loop over SNRs
end % End loop over sentences

d_mean= mean(d, 2); % Model outputs averaged across sentences

%% Fitting

xdata = d_mean';
ydata = Pcorrect_human;

fun = @(a,xdata) 100./(1 + exp(a(1)*xdata + a(2))); %Logistic function

pguess = [0 0]; %starting guess
[fit_param,R,J,CovB,MSE] = nlinfit(xdata,ydata,fun,pguess); % non-linear least squares optimization

% Goodness of fit
ysim= fun(fit_param, xdata);
rsq2 = 1 - sum(R.^2) / sum((ydata- mean(ydata)).^2);
varagout{1} = 0;
%% Plotting

x = linspace(0, 1, 200);
fit_funct= fun(fit_param, x);

figure
scatter(d_mean, Pcorrect_human, 68, 'filled', 'r')
hold on
plot(x, fit_funct,'k', 'LineWidth', 2)
xlabel('Model output (d)'), ylabel('% correct'), title('sCASP mapping for CLUE material')
le= legend(' Nielsen & Dau (2009) data', 'f_C_L_U_E', 'Location', 'southeast');
set(le, 'box', 'off')
text(3,12,{[' R^2 = ',num2str(rsq2,2)]},'fontsize',12,'FontName', 'cambria');
set(gca, 'fontsize',12)


