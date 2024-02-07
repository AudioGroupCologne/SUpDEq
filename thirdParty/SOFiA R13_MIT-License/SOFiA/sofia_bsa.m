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
% B/S/A BEMA Spatial Anti Aliasing R13-0603
%
%       BEMA Bandwidth Extension for Microphone Arrays
%       -> AES Convention 2012, San Francisco USA
%          Convention Paper 8751
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
%
% Pnm = sofia_bsa(Pnm, ctSig, dn, transition, avgBandwidth, fade)
% ------------------------------------------------------------------------
% Pnm               Alias-Free Spatial Fourier Coefficients
%
% ------------------------------------------------------------------------
% Pnm               Spatial Fourier Coefficients
% ctSig             Signal of the center microphone
% dn                Radial filters for the current array configuration
% transition        Highest stable bin
%                   Approx: transition = (NFFT/FS+1) * (N*c)/(2*pi*r)
% avgBandwidth      Averaging Bandwidth in oct
% fade              Fade over {true} or hard cut {false}  
%

function Pnm = sofia_bsa(Pnm, ctSig, dn, transition, avgBandwidth, fade)

disp('SOFiA B/S/A - BEMA Spatial Anti Aliasing R13-0306');

startaverage = floor(transition/(2^(avgBandwidth))); %Averaging BW
dn_inverted  = 1./dn;                                %Inverted Radial Filters

SpatialImage = zeros(size(Pnm,1),1);
AveragePower          = 0;
avgBins               = length(startaverage:transition);
N                     = sqrt(size(Pnm,1))-1;

cnt=0;

%Extraction of the Spatial Image
for n=0:N
    for m=-n:n
        cnt=cnt+1;
        for bincounter = startaverage:transition
            compensatedPnmBin  = Pnm(cnt,bincounter) / dn_inverted(n+1,bincounter);
            SpatialImage(cnt)  = SpatialImage(cnt) + compensatedPnmBin;
            AveragePower       = AveragePower      + abs(compensatedPnmBin).^2;
        end
    end
end

%Normalization of the Image
SumEnergyOfAllModesInFP  = sum(abs(SpatialImage).^2);
powerNormalizationFactor = AveragePower/avgBins;
SpatialImage    = SpatialImage*sqrt(powerNormalizationFactor/SumEnergyOfAllModesInFP);

%Normalization of the Center Signal
ctSigAverage = mean(abs(ctSig(startaverage:transition)).^2);
ctSig        = ctSig/sqrt(ctSigAverage);

%Add spectral information from Center Signal to Spatial Image
Pnm_bema = zeros(size(Pnm));
cnt = 0;

for n=0:N
    for m=-n:n
        cnt=cnt+1;
        for f=startaverage:size(Pnm,2)
            Pnm_bema(cnt, f) = SpatialImage(cnt) * dn_inverted(n+1,f)* ctSig(f);
        end
    end
end

%Phase correction
phaseoffset = angle(Pnm(1, transition)/Pnm_bema(1,transition));
Pnm_bema = Pnm_bema.*exp(1i*phaseoffset);

%Replace high bins with BEMA-coefficients
Pnm = [Pnm(:, 1:transition-1) Pnm_bema(:, transition:end)];

%Fade over original coefficients and BEMA-coefficients
if fade
    PnmFade_O = Pnm(:, startaverage:transition-1);
    PnmFade_S = Pnm_bema(:, startaverage:transition-1);
    
    fader    = linspace(0,1,size(PnmFade_O,2));
    fadeup   = repmat(fader,size(PnmFade_O,1),1); %Make Matrix UP
    fadedown = fliplr(fadeup);                    %Make Matrix DOWN
    
    PnmFade_O = PnmFade_O.*fadedown;
    PnmFade_S = PnmFade_S.*fadeup;
    
    PnmFade = PnmFade_O + PnmFade_S;
    Pnm(:, startaverage:transition-1) = PnmFade;
end

