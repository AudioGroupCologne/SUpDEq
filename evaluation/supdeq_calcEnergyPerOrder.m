%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function results = supdeq_calcEnergyPerOrder( HRTFdataset, plotResults )
%
% This function calculates the energy of the SH-coefficients with respect
% to the transform order N. The function uses AKshEnergy, calculating  
% the energy using Parseval's theorem and normalizing by 1/(4*pi).
%
% Output:
% results               - Results struct with average energy per order N
%                         for SH-coefficients of left/right ear
%
% Input:        
% testHRTFdataset       - Struct with (SH-coefficients) of the test 
%                         HRTF dataset (for example the de-equalized HRTF set
%                         or any other HRTF set like the referenceHRTFdataset) 
%                         for the left (Hl_nm) and right (Hr_nm), channel/ear, 
%                         absolute frequency scale f, transform order N, and FFToversize
% plotResults          -  plotResults true/false. Simple plot...
%                         Default: false
%
% Dependencies: AKtools
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function results = supdeq_calcEnergyPerOrder( HRTFdataset, plotResults )

if nargin < 2
    plotResults = false;
end

%% Calculate energy for SH-coefficients of letft and right ear

[~,~,~,avgEperN_Hl_nm] = AKshEnergy(HRTFdataset.Hl_nm);
[~,~,~,avgEperN_Hr_nm] = AKshEnergy(HRTFdataset.Hr_nm);

%Write results struct
results.avgEperN_Hl_nm = avgEperN_Hl_nm;
results.avgEperN_Hr_nm = avgEperN_Hr_nm;
results.N              = HRTFdataset.N;

%% Plot if intended

if plotResults
   figure;
   plot(0:results.N,10*log10(results.avgEperN_Hl_nm),'k','LineWidth',1.5);
   hold on;
   plot(0:results.N,10*log10(results.avgEperN_Hr_nm),'r','LineWidth',1.5);
   legend('Left','Right');
   grid on;
   title('Total Energy per Order N')
   xlabel('Order N')
   ylabel('Energy [dB]')
end

end

