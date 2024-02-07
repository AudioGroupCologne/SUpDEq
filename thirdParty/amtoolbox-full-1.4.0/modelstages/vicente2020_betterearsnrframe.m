function [weighted_BE_Aud_SNR, BESNR_perFB] = ...
    vicente2020_betterearsnrframe(fs, fc, MaskerSig, targ_left,targ_right, InternalNoise_L, InternalNoise_R, ceiling, weightings)
%VICENTE2020_BETTEREARSNRFRAME calculates the better ear SNR
%   Usage: [weighted_BE_Aud_SNR, BESNR_perFB] = vicente2020_betterearsnrframe(fs, fc, MaskerSig, targ_left,targ_right, InternalNoise_L, InternalNoise_R, ceiling, weightings)
%
%   VICENTE2020_BETTEREARSNRFRAME computes the better ear SNR
%   in all frequency bands centered at fc, considering the audiogram.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/vicente2020_betterearsnrframe.php


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


BESNR_perFB = zeros(1,length(fc));

for n = 1:length(fc)
    % Get the filtered signals 
    LeftMasker = auditoryfilterbank(MaskerSig(:,1),fs,fc(n), 'lavandier2022');        
    RightMasker = auditoryfilterbank(MaskerSig(:,2),fs,fc(n), 'lavandier2022');
    % Compute the effective masker
    EffectiveLeftMasker = max(20*log10(rms(LeftMasker)), InternalNoise_L(n));
    EffectiveRightMasker = max(20*log10(rms(RightMasker)), InternalNoise_R(n));
    % Compute the left and right SNR
    LeftSNR = targ_left(n) - EffectiveLeftMasker;
    RightSNR = targ_right(n) - EffectiveRightMasker;
    % Compute better-ear SNR
    BESNR_perFB(n) = min(ceiling,max(LeftSNR,RightSNR));
end
%SII weightings+integreation across frequency
weighted_BE_Aud_SNR = sum(BESNR_perFB.*weightings');

end


