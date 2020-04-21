%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function ild = supdeq_calcILD( hrir )
%
% This function calculates the ILD (Interaural Level Difference) of a HRIR 
%
% Output:
% ild               - ILD of the HRIR in dB
%
% Input:
% hrir              - [N x 2] HRIR with N samples and 2 channels    
%
% Dependencies: -
%
%   
% (C) 2020 by JMA, Johannes M. Arend  
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function ild = supdeq_calcILD( hrir )

    %Check type of array
    if size(hrir,2) > size(hrir,1)
        hrir = hrir.';
    end

    %Calculate ILD
    ild = 10*log10(sum(abs(hrir(:,1).^2))/sum(abs(hrir(:,2).^2)));

end
