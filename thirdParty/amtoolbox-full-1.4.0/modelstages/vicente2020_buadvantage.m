function [weighted_BUAdv,BUAdv_perFB] = vicente2020_buadvantage(fs, fc, MaskerSig, TargetLeftSpectrum, TargetRightSpectrum, TargetIPD, InternalNoise_L, InternalNoise_R, weightings)
%VICENTE2020_BUADVANTAGE calculates the bmld advantage
%   Usage: [weighted_BUAdv,BUAdv_perFB] = vicente2020_buadvantage(fs, fc, MaskerSig, TargetLeftSpectrum, TargetRightSpectrum, TargetIPD, InternalNoise_L, InternalNoise_R, weightings)
%
%   VICENTE2020_BUADVANTAGE computes the SII-weighted BU advantage in all frequency bands centered at fc, 
%   considering the audiogram.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/vicente2020_buadvantage.php


%   #StatusDoc: Perfect
%   #StatusCode: Good
%   #Verification: Qualified
%   #Requirements: MATLAB
%   #Author: Matthieu Lavandier
%   #Author: Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

BUAdv_perFB = zeros(1,length(fc));

for n = 1:length(fc)
    LeftMasker = auditoryfilterbank(MaskerSig(:,1),fs,fc(n), 'lavandier2022');        
    RightMasker = auditoryfilterbank(MaskerSig(:,2),fs,fc(n), 'lavandier2022');
    if 20*log10(rms(LeftMasker)) > InternalNoise_L(n) &&  20*log10(rms(RightMasker)) > InternalNoise_R(n) &&  TargetLeftSpectrum(n) > InternalNoise_L(n) &&  TargetRightSpectrum(n) > InternalNoise_R(n)
        [MaskerIPD, MaskerCoherence] = local_do_xcorr(LeftMasker, RightMasker, fs, fc(n)); 
        BUAdv_perFB(n) = bmld(MaskerCoherence, TargetIPD(n), MaskerIPD, fc(n));                                       
    end
end
%SII weightings+integreation across frequency
weighted_BUAdv = sum(BUAdv_perFB .* weightings');

end

function [phase, coherence] = local_do_xcorr(left, right, fs, fc)
    [iacc, lags] = xcorr(left,right,round(fs/(fc*2)),'coeff'); %round(fs/(fc*2)) is for conformity with Durlach's 1972 formulation which allows time delays up to 
                                                               %+/- half the period of the channel centre frequency.
    [coherence, delay_samp] = max(iacc);
    phase = fc*2*pi*lags(delay_samp)/fs;
end


