%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [gridData, Npoints, Nmax] = supdeq_lebedev(numOfSP, N)
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
% N             - Instead of numOfSP, the spatial order N of the Lebedev
%                 grid can be passed. If N is passed, numOfSP will be
%                 ignored.
%
% Dependencies: SOFiA toolbox
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [gridData, Npoints, Nmax] = supdeq_lebedev(numOfSP, N)

if nargin == 0
    sofia_lebedev();
end
    
if nargin == 2
    %Get Lebedev sampling grid according to N
    lebN = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65];
    lebNumPoints = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810];
    
    if isempty(find(N == lebN))
        error('No Lebedev grid available with given spatial order N.');
    end
 
    %Get numOfSP according to N
    numOfSP = lebNumPoints(find(N == lebN));
end

[gridData, Npoints, Nmax] = sofia_lebedev(numOfSP,0);
%Convert from rad 2 deg
gridData(:,1:2) = gridData(:,1:2)*180/pi;

end

