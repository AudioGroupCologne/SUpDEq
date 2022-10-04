%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [gridData, Npoints, Nmax] = supdeq_extremal(N)
%
% This function returns an extremal sampling grid with azimuth and 
% elevation in degree plus the respective sampling weights.
% The nodes for 1 <= N <= 100 are from 
% https://web.maths.unsw.edu.au/~rsw/Sphere/Extremal/New/index.html
%
% Output:
% gridData      - Q x 3 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
% Npoints       - Total number of sampling points / nodes [(N+1)^2]
% Nmax          - Highest stable grid order  
%
% Input:
% N             - Spatial order of the extremal grid. 
%
% Dependencies: -
%
% (C) 2020 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [gridData, Npoints, Nmax] = supdeq_extremal(N)

if nargin == 0
    error('Please specify the desired spatial order N!');
elseif N <= 0
    error('Sorry, but N >= 1 is mandatory');
elseif N > 100
    error('Sorry, but we only have extremal grids up to N = 100');
else
    % Load nodes from materials
    load('extremalNodes_1_100.mat');
    
    %Get respective grid data already in correct format
    gridData = extremalNodes{N};
    
    %Write other variables
    Npoints = size(gridData,1);
    Nmax = N;
    
end

end

