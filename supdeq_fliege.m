%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [gridData, Npoints, Nmax] = supdeq_fliege(N)
%
% This function returns a Fliege sampling grid with azimuth and elevation
% in degree plus the respective sampling weights.
%
% Output:
% gridData      - Q x 3 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
% Npoints       - Total number of sampling points / nodes [(N+1)^2]
% Nmax          - Highest stable grid order  
%
% Input:
% N             - Spatial order of the Fliege grid. 
%
% Dependencies: -
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [gridData, Npoints, Nmax] = supdeq_fliege(N)

if nargin == 0
    error('Please specify the desired spatial order N!');
elseif N <= 0
    error('Sorry, but N >= 1 is mandatory');
elseif N > 29
    error('Sorry, but we only have Fliege grids up to N = 29');
else
    % Load nodes from resources
    % Nodes copied from "Spherical Harmonic Transform Library" by Archontis
    % Politis, http://research.spa.aalto.fi/projects/sht-lib/sht.html
    load('fliegeMaierNodes_1_30.mat');
    
    cartCoord = fliegeNodes{N+1}(:,1:3);
    [gridData(:,1), gridData(:,2)] = cart2sph(cartCoord(:,1), cartCoord(:,2), cartCoord(:,3));
    %Weights
    w = fliegeNodes{N+1}(:,4);
    w = w/sum(w);
    gridData(:,3) = w;
    
    %Convert from rad 2 deg
    gridData(:,1:2) = gridData(:,1:2)*180/pi;
    
    %Transform to SUpDEq coordinate system
    az = gridData(:,1);
    az(az<0)=az(az<0)+360;
    gridData(:,1) = az;
    gridData(:,2) = 90-gridData(:,2);
    
    %Write other variables
    Npoints = size(gridData,1);
    Nmax = N;
    
end

end

