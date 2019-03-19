% SOFiA example 7: Impulse response reconstruction on a simulated 
%                  sampled unity plane wave.
% SOFiA Version  : R13-0306

  clear all
  clc
  
% Generate a full audio spectrum plane wave using S/W/G
   
  r    = 0.2;             % Array radius
  ac   = 2;               % Rigid Sphere Array
  FS   = 24000;           % Sampling Frequency
  NFFT = 1024;            % FFT-Bins
  AZ   = 0;               % Azimuth angle
  EL   = pi/2;            % Elevation angle
  
  quadrature_grid = sofia_lebedev(110, 0);  %EXAMPLE GRID LEB110, No Plot
 
  [fftData, kr] = sofia_swg(r, quadrature_grid, ac, FS, NFFT, AZ, EL);

% Spatial Fourier Transform  

  Nsft = 5;                     % Transform order  
  Pnm = sofia_stc(Nsft, fftData, quadrature_grid);

% Make radial filters 

  Nrf  = Nsft;                  % radial filter order
  limit = 150;                  % Amplification Limit (Keep the result
                                % numerical stable at low frequencies)
  
  dn   = sofia_mf(Nrf, kr, ac, limit);
  
% Plane wave decomposition for different look directions
  
  Npdc   = Nsft;                 % Decomposition order
  OmegaL = [0 pi/2; pi/2 pi/2];  % Looking towards the wave and to one side
  
  Y = sofia_pdc(Npdc, OmegaL, Pnm, dn);
  
% Reconstruct impulse responses

  impulseResponses = sofia_tdt(Y,0);
    
% Make IR causal:

  impulseResponses = [impulseResponses(:,end/2+1:end), impulseResponses(:,1:end/2)];
  
% Plot results (Impulse Responses)

  figure(1)  
  subplot(1,2,1)
  plot(impulseResponses')
  title('Impulse response');
  xlabel('Samples');
  ylabel('Amplitude');
  
  axis([-100 NFFT -0.2 1.2])
  grid on

% Plot results (Spectra)

  spectrum = abs(fft(impulseResponses'));
  fscale = FS/2*linspace(0,1,NFFT/2+1);
  
  subplot(1,2,2)
  semilogx(fscale,20*log10(spectrum(1:end/2+1,:)))
  title('Spectrum');
  xlabel('Frequency in Hz')
  ylabel('Magnitude');  
  axis([50 FS/2 -60 30])
  grid on
