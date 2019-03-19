% SOFiA example  : AE9 CAD Reflection Analysis
% SOFiA Version  : R13-0306

clear all
clc
  
skipMatrixGeneration = true;

if ~skipMatrixGeneration
% Step1: Read Array Data from miro object

load AE_Resources/SBS_VSA_110RS_SBC
timeData = SBS_VSA_110RS_SBC.miroToSOFiA;

% Step2: Generate an equal distance matrix: This may take a while...

    ac           = 2;
    frequency    = 2000;
    timeLimit    = 0.300;
    
    eqDistMTX = sofia_makeEqDistMTX(timeData, ac, frequency, timeLimit);
    
    save AE_Resources/eqDistMTX eqDistMTX
end

% Step3: Detect Peaks...

load AE_Resources/eqDistMTX eqDistMTX

dBsensitivity  = 12.5;
surfaceRange   = 20;
timeRange      = 1;
dBthreshold    = 40;

detPeaks = sofia_pkd(eqDistMTX,dBsensitivity,surfaceRange,timeRange,dBthreshold);

% Step4: Visualize Peaks in CAD Model

cadFile = 'AE_Resources/SBS_CAD.mat';

sofia_cadra(detPeaks, cadFile);
disp(' ')
disp('The figure shows the Small Broadcast Studio at the WDR Funkhaus in Cologne, Germany.')
disp('For more information visit: http://www.audiogroup.web.th-koeln.de and go to the section:')
disp('"Spatial audio impulse response compilation captured at the WDR broadcast studios".')

