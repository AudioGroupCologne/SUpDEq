%DEMO_VERHULST2012 Demo of the cochlear model calculating otoacoustic emissions
%  
%   This script computes and plot the otoacoustic emission for a 500 Hz sinusoid,
%   This is generated using the cochlear model described in verhulst2012
%   In particular otoacoustic emissions are computed as the signal difference between the
%   sound pressure at the middle ear with model non linearities and irregularities enabled,
%   and the sound pressure at the middle ear in case of linear model.
%
%   Figure 1: Output of the cochlear model.
%
%
%   See also: verhulst2012
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_verhulst2012.php


%   #Author: Alessandro Altoe (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

fs=48000;
t=0:1.0/fs:0.05;
insig=zeros(length(t),2);
f0=500;
insig(:,1)=sin(2*pi*f0*t);
insig(:,2)=insig(:,1);
spl=[60,60];
irron=[1,0];
normalizeRms=[1 1];
subjectNo=rand();
fc='all';

[V,Y,OAE,CF]=verhulst2012(insig,fs,fc,spl,'normalize',normalizeRms,'subject',subjectNo,'irr',irron);

OtoacousticEmissionPa=OAE(:,1)-OAE(:,2);

figure; 
plot(t.*1e3,OtoacousticEmissionPa);
grid on;
xlabel('Time (ms)');
ylabel('Emission (Pa)');
title('Otoacoustic emission for a 500 Hz sinsuoid');


