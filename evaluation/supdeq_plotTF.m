%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function supdeq_plotTF( tf1, tf2, fs, singleSided)
%
% This function plots magnitude, phase, and group delay 
% of one or two transfer functions
%
% Output:
% Figure with 3 subplots
%
% Input:
% tf1           - Transfer Function 1 - For example HRTF left
% tf2           - Transfer Function 2 - For example HRTF right
% fs            - Sampling rate
%                 Default: 48000
% singleSided   - Boolen with true/false - 
%                 Defines if spectrum is single or both sided
%                 Default: true - single sided spectrum only  
%
% Dependencies: -
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

%%
function supdeq_plotTF( tf1, tf2, fs, singleSided)

plotTF2 = true;
if nargin < 2
    plotTF2 = false;
    tf2 = [];
end

if nargin < 3 || isempty(fs)
    fs = 48000;
end

if nargin < 4 || isempty(singleSided)
    singleSided = true;
end

%Check size of tf
if size(tf1,2) > size(tf1,1)
    tf1 = tf1';
end

if plotTF2
    if size(tf2,2) > size(tf2,1)
        tf2 = tf2';
    end
end

%% Get phase and group delay

NFFTtf1 = length(tf1);
if plotTF2
    NFFTtf2 = length(tf2);
end

if singleSided
    fVec_tf1    = linspace(0,fs/2,NFFTtf1);
    fVec2_tf1   = linspace(0,fs/2,NFFTtf1-1);
    if plotTF2
        fVec_tf2    = linspace(0,fs/2,NFFTtf2);
        fVec2_tf2   = linspace(0,fs/2,NFFTtf2-1);
    end
else
    fVec_tf1    = linspace(0,fs,NFFTtf1);
    fVec2_tf1   = linspace(0,fs,NFFTtf1-1);
    if plotTF2
        fVec_tf2    = linspace(0,fs,NFFTtf2);
        fVec2_tf2   = linspace(0,fs,NFFTtf2-1);
    end
end
    
phitf1  = unwrap(angle(tf1));
if plotTF2
    phitf2  = unwrap(angle(tf2));
end
if singleSided
    grtf1   = -diff(phitf1) / (2*pi) * NFFTtf1/(fs/2);
    if plotTF2
        grtf2   = -diff(phitf2) / (2*pi) * NFFTtf2/(fs/2);
    end
else
    grtf1   = -diff(phitf1) / (2*pi) * NFFTtf1/fs;
    if plotTF2
        grtf2   = -diff(phitf2) / (2*pi) * NFFTtf2/fs;
    end
end

%% Plot

figure;
set(gcf,'Color','w');
linewidth = 1.1;

subplot(2,2,[1,2])
semilogx(fVec_tf1,20*log10(abs(tf1)),'k','Linewidth',linewidth);
if plotTF2
    hold on;
    semilogx(fVec_tf2,20*log10(abs(tf2)),'r','Linewidth',linewidth);
    legend('TF1','TF2','Location','NorthWest');
end
xlim([20, fs/2]);
title('Magnitude Spectrum');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;

subplot(2,2,3)
semilogx(fVec_tf1,phitf1,'k','Linewidth',linewidth);
if plotTF2
    hold on;
    semilogx(fVec_tf1,phitf1,'r','Linewidth',linewidth);
end
xlim([20, fs/2]);
xlabel('Frequency [Hz]');
ylabel('Phase [°]');
title('Phase');
grid on;

subplot(2,2,4)
semilogx(fVec2_tf1,grtf1*1000,'k','Linewidth',linewidth);
if plotTF2
    hold on;
    semilogx(fVec2_tf2,grtf2*1000,'r','Linewidth',linewidth);
end
xlim([20, fs/2]);
title('Group Delay');
xlabel('Frequency [Hz]');
ylabel('Group Delay [ms]');
grid on;

end

