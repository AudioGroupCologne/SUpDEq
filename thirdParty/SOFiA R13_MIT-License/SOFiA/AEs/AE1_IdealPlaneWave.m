% SOFiA example 1: Ideal unity plane wave simulation   
% SOFiA Version  : R13-0306

  clear all
  clc
  
% Generate an ideal plane wave using W/G/C (Wave Generator Core)
  
  N    = 9;     % Order
  r    = 0.5;   % Array radius
  ac   = 2;     % Array configuration, 2: Rigid sphere array
  FS   = 48000; % Sampling Frequency
  NFFT = 128;   % FFT-Bins
  AZ   = pi/3;  % Azimuth angle
  EL   = pi/3;  % Elevation angle

  [Pnm, kr] = sofia_wgc(N, r, ac, FS, NFFT, AZ, EL); 
 
% Make radial filters for the rigid sphere array

  Nrf      = 9;   % radial filter order                
  dn       = sofia_mf(Nrf, kr, ac);
    
% Make visualization MTX  

  Nmtx = 9;
  krIndex = 30;   % Here we select the kr-bin (Frequency) to display. 
  mtxData = sofia_makeMTX(Nmtx, Pnm, dn, krIndex);
  
% Visualization with Style 0 

  figure(1)
  clf();
  
  sofia_visual3D(mtxData, 0); 
        