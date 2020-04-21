%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function r_opt = supdeq_optRadius(headWidth, headHeight, headLength)
%
% This function calculates the optimal radius r (in m) as a linear 
% combination of head width, head height, and head depth according to 
% Algazi et al. [1] or Bahu & Romblom [2]
%
% Output:
% r             - Optimal radius r in m
%
% Input:
% headWith      - Head with in m
% headHeight    - Head height in m
% headDepth     - Head length/depth in m
% method        - 'Algazi' (ILD optimized) or 'Bahu' (ITD optimized)  
%                 'Algazi' is default
%
% Dependencies: -
%
% References:
% [1] - Algazi VR., Avendano C., Duda RO. - Estimation of a spherical-head 
% model from anthropometry. J. Audio Eng. Soc. 2001; 49(6):472-479.
%
%?[2] H. Bahu and R. David, ?Optimization and prediction of the spherical 
% and ellipsoidal ITD model parameters using offset ears,? 
% in Proceedings of the AES International Conference 
% on Spatial Reproduction - Aesthetics and Science, 2018, pp. 1?11.
%   
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function r_opt = supdeq_optRadius(headWidth, headHeight, headLength, method)

if nargin < 4 || isempty(method)
    method = 'Algazi';
end


if strcmp(method,'Algazi')
    %r_opt = 0.51*X1 + 0.019*X2 + 0.18*X3 + 3.2 cm
    %with X1 = head half width, X2 = head half height, X3 = head half length
    r_opt = 0.51*(headWidth/2) + 0.019*(headHeight/2) + 0.18*(headLength/2) + 0.032;
end

if strcmp(method,'Bahu')
    %r_opt = ?0.44*X1 + 0.23*X3 + 3.2 cm
    %with X1 = head half width, X2 = head half height, X3 = head half length
    r_opt = 0.44*(headWidth/2) + 0.23*(headLength/2) + 0.032;
end

end

