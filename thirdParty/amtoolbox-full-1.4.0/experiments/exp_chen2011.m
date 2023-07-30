function varargout=exp_chen2011(varargin)
%EXP_CHEN2011 experiments from Chen et al. 2011
%
%   EXP_CHEN2011(fig) reproduce Fig. no. fig from the Chen et
%   al. 2011 paper.
%
%   The following flags can be specified;
%
%     'fig7'    Reproduce Fig. 7, yields the excitation patterns
%
%     'fig9'    Reproduce Fig. 9, contrasting loudness (sones)
%               with input level (dBSPL)
%
%     'fig10'   Reproduce Fig. 10, comparing the equal loudness contours
%               from chen2011_loudness with the ISO 2006 standard
%
%     'fig11'   Reproduce Fig. 11, displaying loudness as a function of bandwidth
%
%
%   Examples:
%   ---------
%
%   To display Fig. 7 use :
%
%     exp_chen2011('fig7');
%
%   To display Fig. 9 use :
%
%     exp_chen2011('fig9');
%
%   To display Fig. 10 use :
%
%     exp_chen2011('fig10');
%
%   To display Fig. 11 use :
%
%     exp_chen2011('fig11');
%
%   References:
%     Z. Chen, G. Hu, B. R. Glasberg, and C. Moore, Brian. A new model for
%     calculating auditory excitation patterns and loudness for cases of
%     cochlear hearing loss. J. Acoust. Soc. Am., 282(1), 2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_chen2011.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Author : Zhangli Chen: original code
%   #Author : Clara Hollomey (2020): integration in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.flags.type = {'missingflag',...
    'fig7','fig9','fig10','fig11'};
  definput.flags.plot = {'plot','no_plot'};
  
    % Parse input options
  [flags,kv]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

if flags.do_fig7 || flags.do_fig9
    
    OuterEarOpt = 'FreeField';

    CF1k = 1000;
    % LdB1k = [2 10:10:90]';
    LdB1k = [2 5 10:5:90]';

%% ANSI S3.4-2007
    phon_ref = [0 1 2 3 4 5 7.5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];

    sone_ref =  [0.0011 0.0018 0.0028 0.0044 0.0065 0.0077 0.017 0.029 0.07 0.142 0.255 0.422 0.662 0.997 1.461 ...
             2.098 2.970 4.166 5.813 8.102 11.326 15.98 22.929 33.216 48.242 70.362 103.274 152.777 227.855 341.982];

    Ldn1k_Binaural_ANSI = interp1(phon_ref,sone_ref,LdB1k,'linear','extrap');

%% model calculation
    for j = 1:length(LdB1k)
        [Ldn1k_Monaural(j),E1k(:,j),Cam,CF] = chen2011(CF1k,LdB1k(j), OuterEarOpt);
    end

    Ldn1k_Binaural_Model = 2* Ldn1k_Monaural;
    
    
    if flags.do_plot && flags.do_fig9
        figure;
semilogy(LdB1k,Ldn1k_Binaural_ANSI,'k--',LdB1k,Ldn1k_Binaural_Model,'k');
legend('ANSI S3.4-2007','This model');
xlabel('Input level of 1-kHz tone, dB SPL');ylabel('Loundess, sones');
axis([0 100 0.001 1000]);
set(findobj('FontSize',10),'FontSize',12);
set(get(gca,'XLabel'),'FontSize',12);
set(get(gca,'YLabel'),'FontSize',12);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'LineWidth',1');
    elseif flags.do_plot && flags.do_fig7
        figure;
semilogx(CF,10*log10(E1k));
xlabel('CF, Hz');ylabel('Excitation, dB');
axis([20 20000 0 100]);
set(findobj('FontSize',10),'FontSize',12);
set(get(gca,'XLabel'),'FontSize',12);
set(get(gca,'YLabel'),'FontSize',12);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'LineWidth',1');
    end
    

elseif flags.do_fig10
    OuterEarOpt = 'FreeField';

%% 1-kHz as reference
CF1k = 1000;
LdB1k = [2 10:10:90]';
for j = 1:length(LdB1k)
    Ldn1k(j) = chen2011(CF1k,LdB1k(j), OuterEarOpt);
end

%% calculating equal loudness contours
%this is the data from iso226-2003
data = amt_load('chen2011', 'data_chen2011.mat');
freq = data.freq;
spl = data.spl;

f = freq(:,1);
LdB = -10:2:140;

for i = 1:length(f)  
    for j = 1:length(LdB)
        [Ldn(j),Exitation,Cam,CF] = chen2011(f(i),LdB(j), OuterEarOpt);          
    end
    
    LdB_sameLdn(i,:) = interp1(Ldn,LdB,Ldn1k,'linear','extrap');
end

if flags.do_plot
figure;
semilogx(freq(:,1),spl(:,1)-1.6,'k--',freq(:,1),LdB_sameLdn(:,1),'k',...
    freq(:,2:end),spl(:,2:end),'k--',freq(:,2:end),LdB_sameLdn(:,2:end),'k');
legend('ISO 226-2003','This paper');
xlabel('Frequency, Hz');ylabel('LdB, dB SPL');
axis([20 20000 -10 130]);
set(findobj('FontSize',10),'FontSize',12);
set(get(gca,'XLabel'),'FontSize',12);
set(get(gca,'YLabel'),'FontSize',12);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'LineWidth',1');
text(600,44,'40 phons');
end

elseif flags.do_fig11
    %% oringinal data
%%%%%% Fig 14 in Moore et al,1997
NoiseF0 = 1420;
NoiseBw0 = 210;
NoiseLdB01 = 30;
NoiseBw11 = [30 70 150 250 500 1000 3600];
NoiseLdB11 = [30.9 31.7 30.1 30.9 31.7 34.0 34.0];
NoiseLdB02 = 50;
NoiseBw12 = [35 75 150 250 500 1125 3400];
NoiseLdB12 = [51.1 50.3 50.3 50.3 54.2 57.3 61.9];
NoiseLdB03 = 50;
NoiseBw13 = [35 75 150 250 500 1125 3400];
NoiseLdB13 = [79.0 79.8 79.8 79.3 81.4 83.1 86.2];

%% model calculation

OuterEarOpt = 'PDR10';

F0 = 1420;      % Hz

%%%%%%%%%%%%%%%% Step1: loudness of matching noise
FixNoiseBW = 210;          % Hz
MatchNoiseLdB = 20:5:100;

for k = 1:length(MatchNoiseLdB)
    
    FixNoiseSpectrumLevel(k) = MatchNoiseLdB(k) - 10*log10(FixNoiseBW);    % dB SPL/Hz

    fl2 = (-FixNoiseBW+sqrt(FixNoiseBW.^2+4*F0^2))/2;
    fu2 = fl2 + FixNoiseBW-1;
    comp_f = fl2:fu2;
    comp_dB = FixNoiseSpectrumLevel(k) * ones(1,length(comp_f));

    Ldn2(k) = chen2011(comp_f,comp_dB, OuterEarOpt);    
end

%%%%%%%%%%%%%%%% Step2: loudness of variable bandwidth noise
FlexNoiseBW = [25:50:125,150:50:300, 400:100:1000,1250,1500,2000,3000,4000];     % Hz
NoiseLevel = [30 50 80];    % dB SPL

for i = 1:length(NoiseLevel)
    for j = 1:length(FlexNoiseBW)
        FlexNoiseSpectrumLevel(i,j) =  NoiseLevel(i) - 10*log10(FlexNoiseBW(j));    % dB SPL/Hz
        
        fl1 = (-FlexNoiseBW(j)+sqrt(FlexNoiseBW(j).^2+4*F0^2))/2;
        fu1 = fl1 + FlexNoiseBW(j)-1;
        comp_f = fl1:fu1;
        comp_dB = FlexNoiseSpectrumLevel(i,j) * ones(1,length(comp_f));
        
        Ldn1(j) = chen2011(comp_f,comp_dB, OuterEarOpt);
    end
        
    LdB_samLdn(i,:) = interp1(Ldn2,MatchNoiseLdB,Ldn1,'linear','extrap');   
end

if flags.do_plot
figure;
semilogx(FlexNoiseBW,LdB_samLdn','k',NoiseBw11,NoiseLdB11,'ko',...
        NoiseBw12,NoiseLdB12,'ko',NoiseBw13,NoiseLdB13,'ko');
axis([20 5000 20 100]);
xlabel('Variable bandwidth, Hz');ylabel('Level of matching noise, dB SPL');
set(findobj('FontSize',10),'FontSize',12);
set(get(gca,'XLabel'),'FontSize',12,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',12,'Vertical','middle');
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'LineWidth',1');
end
end


