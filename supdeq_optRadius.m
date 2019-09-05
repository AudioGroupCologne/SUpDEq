%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function r_opt = supdeq_optRadius(headWidth, headHeight, headLength)
%
% This function calculates the optimal radius r (in m) as a linear 
% combination of head width, head height, and head depth according to 
% Algazi et al. [1].
%
% Output:
% r             - Optimal radius r in m
%
% Input:
% headWith      - Head with in m
% headHeight    - Head height in m
% headDepth     - Head length/depth in m
%
% Dependencies: -
%
% References:
% [1] - Algazi VR., Avendano C., Duda RO. - Estimation of a spherical-head 
% model from anthropometry. J. Audio Eng. Soc. 2001; 49(6):472-479.
%   
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function r_opt = supdeq_optRadius(headWidth, headHeight, headLength)

%r_opt = 0.51*X1 + 0.019*X2 + 0.18*X3 + 3.2 cm
%with X1 = head half width, X2 = head half height, X3 = head half length

r_opt = 0.51*(headWidth/2) + 0.019*(headHeight/2) + 0.18*(headLength/2) + 0.032;

end

