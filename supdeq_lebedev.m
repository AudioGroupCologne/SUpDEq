%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [gridData, Npoints, Nmax] = supdeq_lebedev(numOfSP)
%
% This function returns a Lebedev sampling grid with azimuth and elevation
% in degree plus the respective sampling weights.
%
% Output:
% gridData      - Q x 3 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
% Npoints       - Total number of sampling points / nodes
% Nmax          - Highest stable grid order  
%
% Input:
% numOfSP       - Number of sampling points / nodes. 
%                 Call supdeq_lebedev() to get a list of valid numbers.
%
% Dependencies: SOFiA toolbox
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [gridData, Npoints, Nmax] = supdeq_lebedev(numOfSP)

if nargin == 0
    sofia_lebedev();
else
    [gridData, Npoints, Nmax] = sofia_lebedev(numOfSP,0);
    %Convert from rad 2 deg
    gridData(:,1:2) = gridData(:,1:2)*180/pi;
end

end

