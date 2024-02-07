% SOFiA example 6: Impulse response reconstruction on a simulated 
%                  ideal unity plane wave.
% SOFiA Version  : R13-0306

  clear all
  clc
  
% Generate a full audio spectrum plane wave using I/W/G
   
  Nwave = 5;       % Wave order
  r     = 1;       % Array radius
  ac    = 2;       % Array configuration: 2-Rigid Sphere
  FS    = 48000;   % Sampling Frequency
  NFFT  = 1024;    % FFT-Bins
  AZ    = 0;       % Azimuth angle
  EL    = pi/2;    % Elevation Angle
    
  [Pnm, kr] = sofia_wgc(Nwave, r, ac, FS, NFFT, AZ, EL);

% Make radial filters 

  Nrf  = Nwave;                  % radial filter order
  dn   = sofia_mf(Nrf, kr, ac);
   
% Running a plane wave decomposition for different look directions
  
  Npdc   = Nwave;                % Decomposition order
  OmegaL = [0 pi/2; pi/2 pi/2];  % Looking towards the wave and to one side
  
  Y = sofia_pdc(Npdc, OmegaL, Pnm, dn);

% Reconstruct impulse responses

  impulseResponses = sofia_tdt(Y);
  
% Plot results (Impulse Responses)

  figure(1)  
  subplot(1,2,1)
  plot(impulseResponses')
  title('Impulse response');
  xlabel('Samples');
  ylabel('Amplitude');
  
  axis([-200 NFFT -0.2 1.2])
  grid on

% Plot results (Spectra)

  spectrum = abs(fft(impulseResponses'));
  fscale = FS/2*linspace(0,1,NFFT/2+1);
  
  subplot(1,2,2)
  semilogx(fscale,20*log10(spectrum(1:end/2+1,:)))
  title('Spectrum');
  xlabel('Frequency in Hz')
  ylabel('Magnitude');  
  axis([50 20000 -60 30])
  grid on
