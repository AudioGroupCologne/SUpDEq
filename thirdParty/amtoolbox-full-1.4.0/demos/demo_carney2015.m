%DEMO_CARNEY2015 peripheral and central neural responses to AM tones
%
%   DEMO_CARNEY2015 generates a neurogram for a specific hearing loss
%   and outputs the responses of the (contralateral) cochlear nucleus and
%   inferior colliculus, contrasting the band-suppressed and the band-enhanced
%   cell output, as predicted by the model devised by Carney et al. (2015)
%
%   Figure 1: Cochlear nucleus and inferior colliculus responses to the generated neurogram
%
% 
%
%   See also: carney2015 carney2015_generateneurogram carney2015_fitaudiogram
%             carney2015_getalphanorm
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_carney2015.php


%   #Author : University of Rochester (UR EAR) team
%   #Author: Clara Hollomey (2021): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

Fs = 100e3; % "model sampling frequency" [samples/sec]

%% the stimulus============================================================
path1 = 'm06iy.wav';
[wav1,Fs_wav1] = amt_load('carney2015', path1);
spl = 65;
Pref = 20*10^(-5);
if ~isnan(spl)
    wav1 = wav1*(Pref*10.^(spl/20)/rms(wav1));
end
stimulus = resample(wav1,Fs,Fs_wav1);

%% generate a neurogram for a specific hearing loss========================

% Hearing loss in dB, to determine Cohc and Cihc
% see also: carney2015_fitaudiogram2 bruce2018 zilany2014
ag_dbloss = [0 0 0 0 0 0 0]; %modify loss here (cooresponding to ag_fs below)
ag_fs = [125 250 500 1e3 2e3 4e3 8e3]; % audiometric frequencies

% determine the species to be used
species = 2; % 1=cat; 2=human

% nerve fiber settings
numCF = 10;
dur = 10;
fiber_num = 5;
CF_range = [200 3000];
fiberType = 3;% AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
  % with  carney2015_generateneurogram
% [psth,neurogram_ft] = carney2015_generateneurogram(stimulus,Fs,species,...
% 				ag_fs,ag_dbloss,CF_num,dur,n,fiber_num,CF_range,fiberType);  
% an_sout = (100000*psth)/fiber_num;

  % with  zilany2014
cf = audspace(CF_range(1),CF_range(2),numCF);
an_sout=zilany2014(stimulus,Fs,cf,'fiberType',fiberType,'nrep',fiber_num);
psth=an_sout;
  % with bruce2018
% cf = audspace(CF_range(1),CF_range(2),numCF);
% numH=12; numM=4; numL=4; 
% kv = {'numH',numH,'numM',numM,'numL',numL,'nrep',10,'psthbinwidth_mr',0.0005}; 
% tmp = bruce2018(stimulus,Fs,cf,kv{:});
% an_sout=tmp.psth_ft/(numH+numM+numL)*100000;
% psth=an_sout;
%% run the model            
BMF = 100;                        
[ic_sout_BE,ic_sout_BS,cn_sout_contra] = carney2015(an_sout,BMF,Fs);

%% plot the output
%account for zero-padding within the model
[sizeN,sizeM] = size(an_sout);
ic_sout_BE    = ic_sout_BE(1:sizeN,1:sizeM);
ic_sout_BS = ic_sout_BS(1:sizeN,1:sizeM);
cn_sout_contra    = cn_sout_contra(1:sizeN,1:sizeM);

t = 1/Fs:1/Fs:length(psth)/Fs;

figure
subplot(3, 1, 1)
plot(t, an_sout)
title('Generated AN neurogram')
xlim([0 length(psth)/Fs])
xlabel('Time index')
grid on
subplot(3, 1, 2)
plot(t, cn_sout_contra)
title('Cochlear nucleus contralateral response')
xlim([0 length(psth)/Fs])
xlabel('Time index')
grid on
subplot(3, 1, 3)
plot(t, ic_sout_BE)
hold on
plot(t, ic_sout_BS)
xlim([0 length(psth)/Fs])
xlabel('Time index')
title('Inferior colliculus responses')
legend('Inferior colliculus BE', 'Inferior colliculus BS')
grid on


