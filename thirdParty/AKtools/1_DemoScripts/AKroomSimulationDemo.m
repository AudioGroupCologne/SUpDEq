% This script shows an example for the AKroomSimulation. This is a
% simulation of a rectangular room that uses a combination of image sources
% for the early reflections and decaying noise for the late reflections.
% The model includes source and receiver directivities, frequency dependend
% wall absorption, and air absorption.
% The model is quite simple, and no big care was taken for an efficient
% implementation. It is intendend for demonstrating basic principles of
% room acoustic simulation. A separate demo for simulating late reflections
% can be found in 2_Tools/RoomSimulation/AKdiffuseReverbTailDemo.m
%
% The model uses AKism.m for calculating the image sources, and
% AKdiffuseReverbTail.m for the late reverberation. Please read the head of
% these functins for more information.
%
% [1] Allen, J. B. & Berkley, D. A.: "Image method for efficiently
%     simulating small-room acoustics." J. Acoust. Soc. Am., 65(4),
%     943-950 (1979).
% [2] Lehmann, E. A. & Johansson, A. M.: "Prediction of energy decay in
%     room impulse responses simulated with an image-source model."
%     J. Acoust. Soc. Am., 124(1), 269-277 (2008).
% [3] Borss, C. & Martin, R.: "An improved parametric Model for perception-
%     based design of virtual acoustics." AES 35th International Conference
%     London, UK, 2009.
% [4] Brinkmann, F. & Erbes, V. & Weinzierl, S.: "Extending the closed form
%     image source model for source directivity." Fortschritte der Akustik
%     DAGA 2018, Munich, Germany, March 2018.
%
% 2018/02 - fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
close all; clear; clc


%% ------------------------------------------------------ scene description
% (the scene geometry can be plotted below)

% Room dimensions in x/y/z direction [m]
% The room is located in the first octant (+++) of a right handed three-
% dimensional carthesian coordinate system. One corner of the room is
% located in the origin of coordinates.
rs.L = [25.66 11.95 7.64];

% wall absorption coefficients alpha:
% [x1,1 x2,1 y1,1 y2,1 z1,1 z2,1
%  x1,2 x2,2 y1,2 y2,2 z1,2 z2,2]
%    .    .    .    .    .    .
%  x1,K x2,K y1,K y2,K z1,K z2,K]
% where x, y, z denote walls with constant x, y, z coordinates, the first
% index denotes the wall position (1 beeing closer to the origin of
% coordinates) and the second index the frequency
rs.alpha = [.24 .11 .22 .22 .37 .17
            .17 .09 .16 .16 .34 .13
            .18 .07 .16 .16 .25 .13
            .25 .07 .21 .21 .28 .14
            .16 .08 .19 .19 .26 .20
            .20 .12 .18 .18 .36 .07
            .19 .12 .19 .19 .35 .12
            .20 .12 .20 .20 .35 .14
            .19 .12 .21 .21 .37 .15
            .20 .13 .22 .22 .40 .17
            .19 .13 .22 .22 .38 .18
            .19 .13 .22 .22 .36 .17
            .17 .12 .21 .21 .35 .17
            .17 .12 .20 .20 .34 .16
            .13 .10 .17 .17 .39 .13
            .14 .10 .19 .19 .41 .15
            .17 .13 .21 .21 .33 .19
            .17 .13 .19 .19 .32 .24
            .19 .15 .20 .20 .30 .26
            .18 .13 .21 .21 .39 .22
            .17 .14 .20 .20 .39 .21
            .17 .14 .21 .21 .40 .21];

% frequencies of absorption coefficients
rs.f = [.063 .08 .1 .125 .16 .2 .25 .315 .4	.5 .63	.8 1 1.25 1.6 2	2.5 3.15 4 5 6.3 8]' * 1e3;

% x/y/z position of the receiver [m]
rs.recPos  = [9.14 5.98 1.23];

% viewing direction of the receiver
% Azimuth:
%      90: rec. points to positive y direction
%       0: rec. points to positive x direction
% -90/270: rec. points to negative y direction
% Elevation:
%      90: rec. points to positive z direction
%     -90: rec. points to negative z direction
rs.recView = [0 0]; 

% N by 2 matrix that gives the rotation of the receiver in azimuth and
% elevation [az1 el1; az2 el2; ... ; azN elN]. This can for example be done
% to simulate a head rotation of a dummy head, and results in N impulse
% responses for the image source model. Positive azimuth values denote a
% counter clockwise rotation. Positive elevation values denote an upwards
% rotation. E.g., [(-30:30)' zeros(61,1)] rotates the receiver for +/-30
% degrees to the left and right.
% pass false to omit the rotation of the receiver
rs.recRot = false;

% specify whether the source position is given in spherical or carthesian
% coordinates
rs.srcCoordinates = 'carthesian';

% - source position(s) in azimuth [deg], elevation [deg], and radius [m]
%   relative to the receiver position rs.RecPos (coordinate convention
%   according to rs.recRot), or
% - absolute source position in x/y/z coordinates [m]
% Each row represents one source.
rs.srcPos = [19    7.48 2.03
             20.43 6.73 2.03
             20.43 5.23 2.03
             19    4.48 1.81];

% Azimuth and elevation orientation of the source(s) following the
% coordinate convention from rs.RecAz, and rs.RecEl
% []                      : sources are oriented towards the receiver
% [az1 el1; ..., azN elN] : sources are oriented according to specified
%                           values (one line per source)
rs.srcView = [];


%% ---------------------------------------- source and receiver description
%  So far source and receiver directivities are only supported in the form
%  of spherical harmonics coefficients. The source data are stored in third
%  octaves, and the receiver data in evenly spaced frequency bins.

% specify the source
% 'OMNI'   - omni-directional source with flat frequency response
% 'QSC-K8' - QSC-K8 two-way PA speaker from the AKtools DemoData folder
rs.src = 'QSC-K8';

% specify the receiver
% 'OMNI'   - omni-directional receiver with flat frequency response
% 'FABIAN' - uses HRTF from the FABIAN database
rs.rec = 'FABIAN';

% pass 'true' to compensate for the diffuse field transfer function of the
% HRTF dataset. This might be a good choice, when listening through
% headphones without applying a headphone compensation filter. In this
% case, the inverse diffuse field HRTF transfer function serves as an
% approximation for a headphone filter.
rs.recDTFcompensation = true;


%% ------------------------------------------------------- model parameters
% specify if you want to use the image source model (ISM, [1-2])
rs.ISM = true;

% specify the truncation of the image source model in a cell array:
% - to include reflections up to order three pass:
%   {'N' 3}
% - to include all refluctions up to 0.1 seconds pass:
%   {'t' 0.1}
% - to truncate after 1.5 times the estimated perceptual mixing time pass:
%   {'tmix' 1.5}
rs.ISMtruncation = {'tmix' 2};

% specify if you want to use stochastic reverb (SR, [3])
% if you are using the ISM and the SR, they will be combined at the point
% of the rs.ISMtruncation
rs.SR = true;

% specifies how the stochastic reverberation is calculated (see
% AKdiffuseReverbTailDemo for options)
rs.SRdecayMode = 'fft_ola';

% specify the dynamic of the SR in dB (a higher dynamic will result in
% longer duration of the SR).
rs.SRdynamic = 90;

% The level of the stochastic reverberation is automatically matched to the
% level of the image source model. A gain coorection in dB can be specified
% if the matching is not good enough. (0 dB = no gain correction)
rs.SRgain = -6;

% The stochastic reverberation reaches it's maximum value at the time
% specified by rs.ISMtruncation. An additional fade in can be applied, if
% desired.
% Specify the fade duration in seconds or pass false to start the
% fade at the position of the first reflection
rs.SRfadeDuration = false;

% specifiy the fade in for the stocahstic reverb
% false - do not fade in
% lin   - linear fade in
% sin   - sin fade in (you can specify the sine power by passing 'sin_2',
%         which for example will apply a squared sine fade in).
rs.SRfadeType = 'lin';


%% ----------------------------------------------------- general parameters

% relative humidity in percent
rs.h_r = 50;

% temperature in degree Celsius
rs.T = 20;

% atmospheric pressure in Pa
rs.p_a = 101325;

% set the speed of sound [m/s]
rs.c = AKspeedOfSound(rs.T, rs.h_r, rs.p_a);

% specify if you want to include air absorption
rs.airAbsorption = true;

% sampling rate in [Hz]
rs.fs = 44100;

% toggle verbosity
rs.verbose = true;

% show a plot of the room geometry and source-receiver configuration
rs.plotScene = true;

%% ------------------------------------------------ plot the scene geometry
if rs.plotScene
    AKroomSimulationPlot(rs)
end


%% ----------------------------------------------------- run the simualtion

[hISM, hSR, hHybrid, ISM, SR] = AKroomSimulation(rs);


%% -------------------------------- plot an example for an impulse response
%  This plot can be used for a quick verification. The level of the
%  stochastic reverb needs to be checked, because the room model estimates
%  it on a very rough basis.

% source number to plot
plotSrc = 1;
% receiver orientation to plot
recRot  = 1;

if rs.ISM && rs.SR

    N = size(hISM,3);

    AKf(40,20)
    for nn = 1:N
        subplot(N,1,nn)
            AKp(hSR(:,nn,plotSrc), 'et2d', 'c', .8)             % stochastic reverb
            AKp(hISM(:,recRot,nn,plotSrc), 'et2d', 'x', [0 1000])    % ISM
            dr = get(gca, 'YLim');
            AKp(hHybrid(:,recRot,nn,plotSrc), 'ed2d', 'x', [0 1000], 'norm_d', db( max(abs(hISM(:,1,nn,plotSrc))) ), 'c', 'b', 'dr', dr) % energy decay of hybrid IR
            hold on
            plot([SR.mixingTime SR.mixingTime]*1000, dr, '--r')
            
            title(['channel ' num2str(nn) ', source 1 (black=ISM, gray=SR, blue=EDC, red=t_{mix})'])
    end

end

clear N nn plotSrc dr recRot


%% ------------------------ plot some reflection paths for the first source

% type 'help AKroomSimulationPlot' for more information
AKroomSimulationPlot(rs, ISM, 1, 'N', 0:1)


%% ----------------------------------------------------------------- listen
%  uncomment the last line to listen to the audio - beware of the level!!!

% choose the IR for auralization ('ISM', 'SR', 'hybrid', or 'anechoic')
auralization = 'hybrid';
% choose audio ('string', or 'noise')
audio        = 'string';
% chose sources to listen to
sources      = 1:4;

% choose the receiver orientation
recID = 1;

% load audio content
if strcmpi(audio, 'string')
    x = audioread('StringQuartet_4channel.wav');
else
    x = [AKpulsedNoise(3, rs.fs, .5, .02, 1) ...
         AKpulsedNoise(3, rs.fs, .5, .02, 1) ...
         AKpulsedNoise(3, rs.fs, .5, .02, 1) ...
         AKpulsedNoise(3, rs.fs, .5, .02, 1)];
end
x = [x; zeros(size(hSR,1), size(x, 2))];

% filter audio with IRs
if any( strcmpi(auralization, {'ISM', 'hybrid', 'anechoic'}) )
    % select the IR
    if strcmpi(auralization, 'ISM')
        h = hISM;
    elseif strcmpi(auralization, 'hybrid')
        h = hHybrid;
    else
        if ~strcmp(rs.rec, 'FABIAN')
            error('Only works if receiver is FABIAN')
        end 
        
        % get source positions
        if strcmpi(rs.srcCoordinates, 'carthesian')
            xyz = rs.srcPos - rs.recPos;
            [az, el, d] = cart2sph(xyz(:,1), xyz(:,2), xyz(:,3));
            az = az / pi * 180;
            el = el / pi * 180;
        else
            az = rs.srcPos(:,1);
            el = rs.srcPos(:,2);
            d = rs.srcPos(:,3);
        end
        
        az = az(sources);
        el = el(sources);
        d  = d(sources);
        
        % get relative gains
        gain = d / max(d);
        
        % get HRIRs
        if rs.recDTFcompensation
            [l, r] = AKhrirInterpolation(az, el, 0, 'measured_sh_dir');
        else
            [l, r] = AKhrirInterpolation(az, el, 0, 'measured_sh_hrir');
        end
        
        % level HRIRs
        l = l .* gain';
        r = r .* gain';
        
        % format data
        h = zeros(size(r,1), 1, 2, numel(sources));
        for nn = 1:numel(sources)
            h(:,:,1,nn) = l(:,nn);
            h(:,:,2,nn) = r(:,nn);
        end
        
        clear xyz az el d gain l r
        
    end
    % convolve
    y = zeros(size(x,1), size(hHybrid,3));
    for nn = 1:numel(sources)
        y = y + fftfilt(squeeze( h(:,recID,:,sources(nn)) ), x(:,sources(nn)));
    end

else
    % convolve
    y = zeros(size(x,1), size(hHybrid,3));
    for nn = 1:numel(sources)
        y = y + fftfilt(squeeze( hSR(:,:,sources(nn)) ), x(:,sources(nn)));
    end
end

% soundsc(y, rs.fs)
