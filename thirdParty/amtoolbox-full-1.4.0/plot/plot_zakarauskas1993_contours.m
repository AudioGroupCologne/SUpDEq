function plot_zakarauskas1993_contours(varargin)
% PLOT_ZAKARAUSKAS1993_CONTOURS Script for plotting HRTF contour plots 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_zakarauskas1993_contours.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Robert Baumgartner (2010), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% b=load('hrtf_M_KEMAR-CIPIC 22kHz.mat');
b=load('hrtf_M_KEMAR-CIPIC 44kHz.mat');

% settings
fs=b.stimPar.SamplingRate;       % sampling frequency
lat=00;         % lateral angle of sagital plane
delta=0;        % lateral variability in degree
q10=10;         % relativ bandwidth of filter bands;  default: 10
space=1.10;     % spacing factor of filter bank (distance of bands)
fstart =1000;   % start frequency of filter bank;     default: 1kHz
fend   =15000;   % end frequency of filter bank;       default: 11kHz

% sorting data and filtering
idx=find(b.posLIN(:,4)>=-(delta+0.01)/2+lat & b.posLIN(:,4)<=(delta+0.01)/2+lat); % median plane with lateral delta-variation; +0.01 because of rounding errors in posLIN
[pol,polidx]=unique(b.posLIN(idx,5));   % sorted polar angles
sagir=double(b.hM(:,idx,:));            % unsorted impulse responses on sagital plane
ir=sagir(:,polidx,:);                   % sorted
[y2n,y2n,fn]=butterfb(ir,ir,q10,fstart,fend,fs,space); % filter bank
y2n=y2n-max(max(max(y2n)))+14;

% plots
do=0;
figure;
[C,h] = contour(fn(1:end-do),pol,dif(squeeze(y2n(:,:,1)),do)');
set(gca,'XScale','log')
set(gca,'XLim',[1000 10000])
set(gca,'XTick',1000:1000:10000)
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10'})
set(gca,'YLim',[-90 270])
set(gca,'YTick',[-90 -50:50:250 270])
set(gca,'YTickLabel',{'','-50','0','50','100','150','200','250',''})
set(gca,'YMinorTick','on')
set(h,'LevelStep',2)
colormap hot
c=colormap;
colormap(flipud(c));
axis square
clabel(C,h,[-12 -4 0 4 8 12 14])
xlabel('Frequency in kHz');ylabel('Elevation in degrees');
title(sprintf('X_n  ;  Q_{10}=%d  ;  c.f. Fig. 4 of Zakarauskas et al. (1993)',q10))

do=1;
figure;
[C,h] = contour(fn(1:end-do),pol,dif(squeeze(y2n(:,:,1)),do)');
set(gca,'XScale','log')
set(gca,'XLim',[1000 10000])
set(gca,'XTick',1000:1000:10000)
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10'})
set(gca,'YLim',[-90 270])
set(gca,'YTick',[-90 -50:50:250 270])
set(gca,'YTickLabel',{'','-50','0','50','100','150','200','250',''})
set(gca,'YMinorTick','on')
set(h,'LevelStep',2)
colormap(flipud(c));
axis square
clabel(C,h,-8:4:12)
xlabel('Frequency in kHz');ylabel('Elevation in degrees');
title(sprintf('D_n  ;  Q_{10}=%d  ;  c.f. Fig. 5 of Zakarauskas et al. (1993)',q10))

do=2;
figure;
[C,h] = contour(fn(1:end-do),pol,dif(squeeze(y2n(:,:,1)),do)');
set(gca,'XScale','log')
set(gca,'XLim',[1000 10000])
set(gca,'XTick',1000:1000:10000)
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10'})
set(gca,'YLim',[-90 270])
set(gca,'YTick',[-90 -50:50:250 270])
set(gca,'YTickLabel',{'','-50','0','50','100','150','200','250',''})
set(gca,'YMinorTick','on')
set(h,'LevelStep',2)
colormap(flipud(c));
axis square
clabel(C,h,-5:5:20)
xlabel('Frequency in kHz');ylabel('Elevation in degrees');
title(sprintf('C_n  ;  Q_{10}=%d  ;  c.f. Fig. 6 of Zakarauskas et al. (1993)',q10))


