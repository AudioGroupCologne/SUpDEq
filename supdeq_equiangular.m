%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [gridData, Npoints, Nmax] = supdeq_equiangular(N)
%
% This function returns a equiangular sampling grid with azimuth and 
% elevation in degree plus the respective sampling weights.
%
% Output:
% gridData      - Q x 3 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
% Npoints       - Total number of sampling points / nodes [4*(N+1)^2]
% Nmax          - Highest stable grid order  
%
% Input:
% N             - Spatial order of the equiangular grid. 
%
% Dependencies: -
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [gridData, Npoints, Nmax] = supdeq_equiangular(N)

if nargin == 0
    error('Please specify the desired spatial order N!');
elseif N <= 0
    error('Sorry, but N >= 1 is mandatory');
elseif N > 44
    error('Sorry, but we only have equiangular grids up to N = 44');
else
    % Load nodes from resources
    % Nodes calculated with ITA toolbox (ita_sph_sampling_equiangular)
    load('equiangularNodes_1_44.mat');
    
    %Get respective grid data already in correct format
    gridData = equiangularNodes{N};
    
    %Write other variables
    Npoints = size(gridData,1);
    Nmax = N;
    
end

end

