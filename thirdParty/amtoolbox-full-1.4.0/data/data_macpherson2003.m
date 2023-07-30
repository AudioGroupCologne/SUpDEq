function data = data_macpherson2003
%DATA_MACPHERSON2003  Listener averages of polar error rates
%   Usage: data = data_macpherson2003
%
%   Output parameters:
%     data     : struct containing the data
%
%   DATA_MACPHERSON2003 returns listener averages of polar error rates
%   (PERs in %) from Macpherson & Middlebrooks (2003).
%
%   The data struct contains the following fields:
%
%     'density'  probed ripple densities in ripples/oct (Exp. I)
%     'depth'    probed ripple depths (peak-to-trough) in dB (Exp. II)
%     'phase'    probed ripple phases in radians (Exp. III)
%     'pe_flat'  PER for flat spectrum
%     'pe_exp1'  increase in PER as a function of ripple density at 
%                a ripple depth of 40dB (col. 1: 0-phase, col. 2: pi-phase)
%     'pe_exp2'  increase in PER as a function of ripple depth at
%                a ripple density of 1 ripple/octave
%                (col. 1: 0-phase, col. 2: pi-phase)
%     'pe_exp3'  increase in PER as a function of ripple phase at a
%                ripple density of 1 ripple/octave and a ripple depth
%                of 40dB
%
%
%   References:
%     E. A. Macpherson and J. C. Middlebrooks. Vertical-plane sound
%     localization probed with ripple-spectrum noise. J. Acoust. Soc. Am.,
%     114:430--445, 2003.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_macpherson2003.php


%   #Author: Robert Baumgartner (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Fig. 6 Ripple-spectrum error rate as a function of ripple density and phase
data.density = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8]; % ripples/oct
z40db = zeros(6,10);  % 0-phase, 40dB
pi40db = zeros(6,10); % pi-phase, 40dB
zflat = zeros(6,10);  % 0-phase, flat pectrum
piflat = zeros(6,10); % pi-phase, flat pectrum
% S63
z40db(1,:) =   [26,14,35,30,14,12,2,7,0,0];  % filled squares
pi40db(1,:) =  [0,27,35,16,30,27,17,2,0,0];  % filled triangles
zflat(1,:) =   zeros(1,10);                  % open squares
piflat(1,:) =  [zeros(1,9),5];               % open triangles
% S64
z40db(2,:) =   [14,10,42,80,36,20,10,16,12,4];  % filled squares
pi40db(2,:) =  [10,42,37,26,22,20,7,2,0,0];  % filled triangles
zflat(2,:) =   [0,1,0,0,3,0,0,0,0,0];        % open squares
piflat(2,:) =  [0,0,7,0,0,0,0,0,0,0];        % open triangles
% S65
z40db(3,:) =   [36,16,25,32,24,16,24,7,8,4]; % filled squares
pi40db(3,:) =  [0,28,28,30,35,28,8,8,0,0];   % filled triangles
zflat(3,:) =   [0,0,2,5,4,2,0,5,5,2];        % open squares
piflat(3,:) =  [0,0,2,0,0,5,2,0,0,0];        % open triangles
% S66
z40db(4,:) =   [22,26,32,24,22,5,12,6,2,0];  % filled squares
pi40db(4,:) =  [0,28,26,35,20,0,0,4,2,0];    % filled triangles
zflat(4,:) =   [0,2,0,2,0,2,0,0,0,0];        % open squares
piflat(4,:) =  [2,0,0,0,0,0,0,0,2,0];        % open triangles
% S67
z40db(5,:) =   [24,30,44,60,26,64,28,34,26,16];   % filled squares
pi40db(5,:) =  [12,26,38,28,50,30,38,34,24,18];   % filled triangles
zflat(5,:) =   [6,4,8,6,12,4,6,8,18,10];          % open squares
piflat(5,:) =  [6,6,14,4,6,4,4,6,6,6];            % open triangles
% S77
z40db(6,:) =   [28,34,54,40,54,42,28,24,30,20];   % filled squares
pi40db(6,:) =  [30,42,44,50,40,34,32,42,34,18];   % filled triangles
zflat(6,:) =   [14,14,18,12,22,18,16,16,18,20];   % open squares
piflat(6,:) =  [4,18,14,14,16,14,18,12,16,16];  	% open triangles

data.pe_exp1(:,:,2) = pi40db;
data.pe_exp1(:,:,1) = z40db;


%% Fig. 9: Ripple-spectrum error rate as a function of ripple depth and phase
data.depth = 10:10:40; % ripple depth in dB
flatm = zeros(6,1);  % flat mean
z1rip = zeros(6,4);  % 0-phase, 1 ripple/oct
pi1rip = zeros(6,4); % pi-phase, 1 ripple/oct
z0rip = zeros(6,4);  % 0-phase, flat spectrum
pi0rip = zeros(6,4); % pi-phase, flat spectrum
% S63
flatm(1) =     0;             % open diamond
z1rip(1,:) =   [0,0,20,30];   % filled squares
pi1rip(1,:) =  [0,2,6,16];    % filled triangles
z0rip(1,:) =   [0,0,0,0];     % open squares
pi0rip(1,:) =  [0,0,4,0];     % open triangles
% S64
flatm(2) =     1;             % open diamond
z1rip(2,:) =   [16,44,62,80]; % filled squares
pi1rip(2,:) =  [6,4,16,26];   % filled triangles
z0rip(2,:) =   [0,0,0,0];     % open squares
pi0rip(2,:) =  [0,2,4,0];     % open triangles
% S65
flatm(3) =     2;             % open diamond
z1rip(3,:) =   [0,2,16,32];   % filled squares
pi1rip(3,:) =  [0,2,25,30];   % filled triangles
z0rip(3,:) =   [0,0,4,5];     % open squares
pi0rip(3,:) =  [5,0,2,0];     % open triangles
% S66
flatm(4) =     1;             % open diamond
z1rip(4,:) =   [4,8,18,24];   % filled squares
pi1rip(4,:) =  [6,4,14,35];   % filled triangles
z0rip(4,:) =   [0,2,0,2];     % open squares
pi0rip(4,:) =  [2,0,2,0];     % open triangles
% S67
flatm(5) =     10;            % open diamond
z1rip(5,:) =   [24,28,48,60]; % filled squares
pi1rip(5,:) =  [12,22,14,28]; % filled triangles
z0rip(5,:) =   [12,12,12,6];  % open squares
pi0rip(5,:) =  [18,12,4,4];   % open triangles
% S77
flatm(6) =     14;            % open diamond
z1rip(6,:) =   [22,40,48,40]; % filled squares
pi1rip(6,:) =  [20,38,43,50]; % filled triangles
z0rip(6,:) =   [14,16,14,12]; % open squares
pi0rip(6,:) =  [10,14,22,14];     % open triangles

data.pe_flat = flatm;  % mean error rate for flat spectrum
data.pe_exp2(:,:,2) = pi1rip;
data.pe_exp2(:,:,1) = z1rip;

%% Fig. 11: Ripple-spectrum error rate as a function of ripple phase, density 1 ripple/oct, depth 40dB
data.phase = -3/4*pi:pi/4:pi; % ripple phase in radians
exp3rip  = zeros(5,8);  % 1 ripple/oct with a depth of 40dB
exp3flat = zeros(5,8);   % flat spectrum
% S64
exp3rip(1,:) =  [16,32,60,80,48,18,12,27];   % filled symbols
exp3flat(1,:) = [4,0,0,0,0,2,0,0];           % open symbols
% S65
exp3rip(2,:) =  [30,28,32,32,34,44,48,30];   % filled symbols
exp3flat(2,:) = [0,2,2,4,0,2,0,0];           % open symbols
% S66
exp3rip(3,:) =  [14,16,18,24,14,22,30,34];   % filled symbols
exp3flat(3,:) = [0,0,2,2,0,0,0,0];           % open symbols
% S67
exp3rip(4,:) =  [32,48,42,60,52,50,38,28];   % filled symbols
exp3flat(4,:) = [8,12,8,6,10,16,10,4];       % open symbols
% S77
exp3rip(5,:) =  [34,40,24,40,18,28,24,50];   % filled symbols
exp3flat(5,:) = [24,18,22,12,24,26,24,14];   % open symbols

data.pe_exp3 = exp3rip;

end


