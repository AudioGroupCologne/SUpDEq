% /// ASAR/MARA Research Group
%  
% Technology Arts Sciences TH Köln
% Technical University of Berlin
% Deutsche Telekom Laboratories
% University of Rostock
% WDR Westdeutscher Rundfunk
% IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis
% 
% F/D/T frequency domain transform R12-0603
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% [fftData, kr, f, ctSig] = sofia_fdt(timeData, FFToversize,
%                                               firstSample, lastSample)
% ------------------------------------------------------------------------
% fftData           Frequency domain data ready for the Spatial
%                   Fourier Transform (SOFiA S/T/C Core)
%                   ! To save processing power the SOFiA chain always
%                     works on the half-sided FFT spectrum only (NFFT/2+1)
% kr                kr-Values of the delivered data (required e.g. for the
%                   modal radial filters SOFiA M/F)
% f                 Absolute frequency scale. Not really required but good
%                   for control purposes or scaling a spectrum plot.
% ctSig             Center signal. If timeData carries a center transducer
%                   time domain signal, the frequency domain signal is
%                   returned.
% ------------------------------------------------------------------------
% timeData          Struct with minimum fields:
%
%                   * .impulseResponses     [Channels X Samples]
%                   * .FS
%                   * .radius               Array radius
%                   * .averageAirTemp       Temperature in [°C]
%
%                   (* .centerIR            [1 x Samples] )
%
% FFToversize       FFToversize rises the FFT Blocksize.   [default = 2]
%                   A FFT of the blocksize (FFToversize*NFFT) is applied
%                   to the time domain data,  where  NFFT is determinded
%                   as the next power of two of the signalSize  which is
%                   signalSize = (lastSample-firstSample).
%
%                   The function will pick a window of
%                   (lastSample-firstSample) for the FFT:
%
% firstSample       First time domain sample to be included
% lastSample        Last time domain sample to be included
%                   If firstSample and lastSample are not defined
%                   the full IR will be transformed:
%                   [ default: firstSample = 1;
%                     lastSample = size(timeData.impulseResponses,2);]
%
% Call this function with a running window (firstSample+td->lastSample+td)
% iteration increasing td to obtain time slices. This way you resolve the
% temporal information within the captured sound field.
%

function [fftData, kr, f, ctSig] = sofia_fdt(timeData, FFToversize, firstSample, lastSample)

disp('SOFiA F/D/T - Frequency Domain Transform R13-0306');

try
    impulseResponses = timeData.impulseResponses;
    FS               = timeData.FS;
    temperature      = timeData.averageAirTemp;
    radius           = timeData.radius;
catch
    fftData=[];
    kr=[];
    f=[];
    error('Input data not valid');
    return
end

if nargin < 2
    FFToversize = 1;
end

FFToversize = round(FFToversize);
if FFToversize < 1
    error('FFToversize must be greater or equal 1')
end

if nargin < 4
    firstSample = 1;
    lastSample  = size(impulseResponses,2);
end

if lastSample < firstSample
    error('lastSample must be greater than firstSample')
end

if firstSample < 1
    error('firstSample must be greater than zero')
end

if lastSample > size(impulseResponses,2)
    error('lastSample out of range')
end


totalSamples     = lastSample - firstSample + 1;
impulseResponses = impulseResponses(:, firstSample:lastSample);
NFFT             = 2^(nextpow2(totalSamples));
fftData          = fft(impulseResponses, NFFT * FFToversize, 2); 
fftData          = fftData(:,1:end/2+1);
fftData          = cast(fftData ,'double');

if isfield(timeData,'centerIR')
    centerIR         = timeData.centerIR(:,firstSample:lastSample);
    ctSig            = fft(centerIR, NFFT * FFToversize, 2);
    ctSig            = ctSig(:,1:end/2+1);
    ctSig            = cast(ctSig ,'double');
else
    ctSig            = [];
end

f  = linspace(0, FS/2, (NFFT*FFToversize)/2+1);
c  = 331.5+0.6*temperature;
k  = 2*pi*f/c;
kr = k*radius;



