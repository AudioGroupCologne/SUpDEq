% SOFiA example 8: Superdirective impulse responses   
% SOFiA Version  : R13-0306

  clear all
  clc 

% Read miro/VariSphear dataset
addpath('AE_Resources');

load SOFiA_A5

SOFiA_A5 = SOFiA_A5.setResampling(22050); %Downsample in order to avoid spatial aliasing
timeData = SOFiA_A5.miroToSOFiA();

% Transform time domain data to frequency domain and generate kr-vector
  
  [fftData, kr, f] = sofia_fdt(timeData);

% Spatial Fourier Transform

  Nsft = 5;
  Pnm  = sofia_stc(Nsft, fftData, timeData.quadratureGrid);

% Radial Filters for a rigid sphere array

  Nrf      = Nsft;   % radial filter order               
  maxAmp   = 10;     % Maximum modal amplification in [dB]
  ac       = 2;      % Array configuration: 2 = Rigid Sphere 
  
  dn                         = sofia_mf(Nrf , kr, ac, maxAmp);
  [dn, kernelSize, latency]  = sofia_rfi(dn); % Radial Filter Improvement

% Plane wave decomposition for directions given by OmegaL:

  Npdc     = Nsft;   %Plane wave decomposition order
  OmegaL   = [0 pi/2; pi pi/2];
  Y        = sofia_pdc(Npdc, OmegaL, Pnm, dn);

% Reconstruct directional impulse responses
  impulseResponses = sofia_tdt(Y, 1/8, timeData.downSample);

  figure(1)
  
  tscale = linspace(0,size(impulseResponses,2)/(timeData.FS*timeData.downSample), size(impulseResponses,2));
  plot(tscale, impulseResponses')
  axis([-0.05 0.8 -0.32 0.32])
  title('Directional Impulse Responses')
  xlabel('Time in s')
  ylabel('Signal Amplitude')
  
  disp(' ');
  disp('You like to get more measured array data sets? Visit: http://www.audiogroup.web.th-koeln.de') 
