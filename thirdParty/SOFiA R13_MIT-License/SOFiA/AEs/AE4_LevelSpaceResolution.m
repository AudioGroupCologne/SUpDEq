% SOFiA example 4: Level/Space Resolution
% SOFiA Version  : R13-0306

clear all
clc
  
% Read miro/VariSphear dataset
addpath('AE_Resources');

load SOFiA_A2
timeData = SOFiA_A2.miroToSOFiA(); 

% Transform time domain data to frequency domain and generate kr-vector   

[fftData, kr, f] = sofia_fdt(timeData);
    
% Spatial Fourier Transform

Nsft = 5;
Pnm  = sofia_stc(Nsft, fftData, timeData.quadratureGrid);

% Radial Filters for a rigid sphere array 

Nrf      = Nsft;    % radial filter order              
maxAmp   = 10;      % Maximum modal amplification in [dB]
ac       = 2;       % Array configuration: 2 = Rigid Sphere 

dn       = sofia_mf(Nrf , kr, ac, maxAmp); % radial filters 
dn       = sofia_rfi(dn);                  % radial filter improvement

% Make MTX  

Nmtx = Nsft;
krIndex = 128;   % Here we select the kr-bin (Frequency) to display. 
mtxData = sofia_makeMTX(Nmtx, Pnm, dn, krIndex);
      

% Plot the response

figure(1)
clf();

sofia_visual3D(mtxData, 0);
view(90,0)    
        
disp(' ');
disp(['The plot shows the response at a frequency of ',num2str(round(10*f(krIndex))/10),'Hz']);