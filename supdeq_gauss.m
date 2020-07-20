%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [gridData, Npoints, Nmax] = supdeq_gauss(N)
%
% This function returns a Gaussian sampling grid with azimuth and elevation
% in degree plus the respective sampling weights.
%
% Output:
% gridData      - Q x 3 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
% Npoints       - Total number of sampling points / nodes [2*(N+1)^2]
% Nmax          - Highest stable grid order  
%
% Input:
% N             - Spatial order of the gauss grid. Twice as many samples 
%                 are always distributed along the azimuth compared to 
%                 the elevation.
%
% Dependencies: SOFiA toolbox
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [gridData, Npoints, Nmax] = supdeq_gauss(N)

if nargin == 0
    error('Please specify the desired spatial order N!');
else
    AZnodes = 2*(N+1);
    ELnodes = AZnodes/2;
    [gridData, Npoints, Nmax] = sofia_gauss(AZnodes,ELnodes,0);
    %Convert from rad 2 deg
    gridData(:,1:2) = gridData(:,1:2)*180/pi;
end

end

