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
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_verhulst2012.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

[V,Y,E,CF]=verhulst2012(insig,fs,fc,spl,normalizeRms,subjectNo,irron);

OtoacousticEmissionPa=E(:,1)-E(:,2);

figure; 
plot(t.*1e3,OtoacousticEmissionPa);
grid on;
xlabel('Time (ms)');
ylabel('Emission (Pa)');
title('Otoacoustic emission for a 500 Hz sinsuoid');
