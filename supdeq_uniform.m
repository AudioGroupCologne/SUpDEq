%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function gridData = supdeq_uniform(numOfSP)
%
% This function returns a uniform sampling grid based on a random normal
% distribution of 'numOfSP' directions
%
% Output:
% gridData      - Q x 3 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
%
% Input:
% numOfSP       - Number of sampling points / nodes. 
%                 Call supdeq_lebedev() to get a list of valid numbers.
%
% Dependencies: -
%
% (C) 2020 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function gridData = supdeq_uniform(numOfSP)

sp = randn(numOfSP,3);
sp = bsxfun(@rdivide,sp,sqrt(sum(sp.^2,2)));

[az,el,~] = cart2sph(sp(:,1),sp(:,2),sp(:,3));
az = az/pi*180; az = mod(az,360);
el = el/pi*180; el = 90-el;
w = ones(numOfSP,1)*(1/numOfSP);
gridData = [az,el,w];

end

