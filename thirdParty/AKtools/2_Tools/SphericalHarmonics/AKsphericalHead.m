% [h, offCenterParameter] = AKsphericalHead(sg, ear, offCenter, a, r_0, Nsh, Nsamples, fs, c)
%
% calculates head-realated impulse responses (HRIRs) of a spherical head
% model with offset ears using the formulation from according to [1]. HRIRs
% are calculated by dividing the pressure on the sphere by the free field
% pressure of a point source in the origin of coordinates.
%
% See AKsphericalHeadDemo.m for example use cases
% Additional information on the model itself can be found in
% 2_Tools/SphericalHarmonics/AKsphericalHead.pdf
%
%
% I N P U T:
% sg        - [N x 3] matrix with the spatial sampling grid, where the
%             first column specifies the azimuth (0 deg. = front,
%             90 deg. = left), the second column the elevation
%             (90 deg. = above, 0 deg. = front, -90 deg. = below), and the
%             third column the radius [m]. If only two columns are given
%             the radius is set to a*100 (see below)
%             (default: AKgreatCircleGrid(90:-10:-90, 10, 90) )
% ear       - four element vector that specifies position of left and right
%             ear: [azimuth_l elevation_l azimuth_r elevation_r]. If only
%             two values are passed, symmetrical ears are assumed.
%             (defualt = [85 -13], average values from [2], Tab V,
%             condition All, O-A) 
% a         - radius of the spherical head in m (default = 0.0875)
% r_0       - distance of the free-field point source in m that used as
%             reference (by default the radius from the sampling grid is
%             taken: r_0 = sg(1,3) )
% offCenter - false   : the spherical head is centered in the coordinate
%                       system (default)
%             true    : the interaural axis (i.e., the connection between
%                       the two ears) is centered. This is done be
%                       averaging the ear azimuth and elevation, and
%                       and a translation the sampling grid
%             [x y z] : x/y/z coordinates of the center of the sphere [m].
%                       E.g., [-4e3 1e3 2e3] moves the spherical head 4 mm
%                       to the back, 1 mm to the left (left ear away from
%                       the origin of coordinates) and 2 mm up.
% Nsh       - spherical harmonics order (default = 100)
% Nsamples  - length in samples (default = 1024)
% fs        - sampling rate in Hz (default = 44100)
% c         - speed of sound [m/s] (default = 343)
%
%
% O U T P U T
% h                  - spherical head impulse responses given in matrix of
%                      size [Nsamples x N x 2]: Left ear = h(:,:,1),
%                      right ear = h(:,:,2)
% offCenterParameter - spherical head model parameters after translation
%                      and changing the ear position (if applied)
%                      ear    : new ear position (see above)
%                      sg     : new sampling grid (see above)
%                      r      : radius for each point of sg
%                      azRot  : rotation above z-axis (azimuth) that was
%                               applied to get the new ear azimuth
%                      elRot  : rotation above x-axis (elevation) that was
%                               applied to get the new ear elevation
%                      xTrans : translation of the spherical head in x-
%                               direction, that was applied to center the
%                               interaural axis
%                      zTrans : translation of the spherical head in z-
%                               direction, that was applied to center the
%                               interaural axis
% 
%
% [1] R. O. Duda and W. L. Martens "Range dependence of the response of a
%     spherical head model." J. Acoust. Soc. Am., 104(5), 3048-3058 (1998).
% [2] H. Ziegelwanger and P. Majdak "Modeling the direction-continuous
%     time-of-arrival in head-related transfer functions" J. Acoust. Soc.
%     Am., 135(3), 1278-1293 (2014).
%
%
% 12/2017 - fabian.brinkmann@tu-berlin.de, kokabi@campus.tu-berlin.de,
%           Silke-Boegelein@web.de,

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [h, offCenterParameter] = AKsphericalHead(sg, ear, offCenter, a, r_0, Nsh, Nsamples, fs, c)

%% -------------------------------------------------- set default parameter
if ~exist('sg', 'var')
    sg = AKgreatCircleGrid(90:-10:-90, 10, 90);
end
if ~exist('c', 'var')
    c = 343;
end
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('Nsamples', 'var')
    Nsamples = 1024;
end
if ~exist('Nsh', 'var')
    Nsh = 100;
end
if ~exist('a', 'var')
    a = .0875;
end
if ~exist('r_0', 'var')
    if size(sg, 2) < 3
        r_0 = 100*a;
    else
        r_0 = sg(1,3);
    end
end
if ~exist('offCenter', 'var')
    offCenter = false;
end
if ~exist('ear', 'var')
    ear = [85 -13];
end

% check format of the sampling grid
if size(sg, 2) < 3
    sg = [sg r_0*ones(size(sg,1), 1)];
end

% check format of ear vector
if numel(ear) == 2
    ear = [ear 360-ear(1) ear(2)];
end

%% -------------------------------- rotate and translate the spherical head
% center the interaural axis
if isnumeric(offCenter)
    
    % sampling grid in carthesian coordinates
    [sgX, sgY, sgZ] = sph2cart(sg(:,1)/180*pi, sg(:,2)/180*pi, sg(:,3));
    
    % translate the sampling grid
    sgX = sgX - offCenter(1);
    sgY = sgY - offCenter(2);
    sgZ = sgZ - offCenter(3);
    
    % translted sampling grid in spherical coordinates
    [sgAz, sgEl, sgR] = cart2sph(sgX, sgY, sgZ);
     sg                = [sgAz/pi*180 sgEl/pi*180 sgR];
     sg                = round(sg*10000) / 10000;
     sg(:,1)           = mod(sg(:,1), 360);
    
    % save parameter
    offCenterParameter.sg = sg;
    
    clear sgX sgY sgZ sgAz sgEl sgR
    
elseif offCenter
    
    % check if the ear azimuths are symmetrical
    if ear(1) ~= 360-ear(3)
        earAz = mean([ear(1) 360-ear(3)]);
        offCenterParameter.azimuthRotation = earAz - ear(1);
        ear([1 3]) = [earAz 360-earAz];
        clear earAz
    else
        offCenterParameter.azimuthRotation = 0;
    end
    
    % check if the ear elevations are symmetrical
    if ear(2) ~= ear(4)
        earEl = mean(ear([2 4]));
        offCenterParameter.elevation         = earEl;
        offCenterParameter.elevationRotation = earEl - ear(2);
        ear([2 4]) = earEl;
        clear earEl
    else
        offCenterParameter.elevationRotation = 0;
    end
    
    % sampling grid in carthesian coordinates
    [sgX, sgY, sgZ] = sph2cart(sg(:,1)/180*pi, sg(:,2)/180*pi, sg(:,3));
    doTranslate     = false;
    
    % check for translation in x-direction (front/back)
    if ear(1) ~= 90
        offCenterParameter.xTranslation = sind(ear(1)-90)*a;
        doTranslate = true;
        sgX         = sgX - offCenterParameter.xTranslation;
    else
        offCenterParameter.xTranslation = 0;
    end
    
    % check for translation in z-direction (up/down)
    if ear(2) ~= 0
        offCenterParameter.zTranslation = sind(-ear(2))*a;
        doTranslate = true;
        sgZ         = sgZ - offCenterParameter.zTranslation;
    else
        offCenterParameter.zTranslation = 0;
    end
    
    % transform grid to spherical coordinates again
    if doTranslate
        [sgAz, sgEl, sgR] = cart2sph(sgX, sgY, sgZ);
        sg                = [sgAz/pi*180 sgEl/pi*180 sgR];
        sg                = round(sg*10000) / 10000;
        sg(:,1)           = mod(sg(:,1), 360);
    end
    
    offCenterParameter.sg  = sg;
    offCenterParameter.ear = ear;
    
    clear earAz earEl sgX sgY sgZ sgAz sgEl sgR doTranslate
    
else
    
    offCenterParameter = false;
    
end

% check parameter values
if any(sg(:,3)<a) || r_0 < a
    error('AKsphericalHead:Input', 'sg(:,3), and r_0 must be smaller than a')
end

%% --------------------------------------------------- spherical head model
% calculate great circle distances between the sampling grid and the ears
gcd = [acosd( sind(sg(:,2))*sind(ear(2)) + cosd(sg(:,2))*cosd(ear(2)) .* cosd(sg(:,1)-ear(1)) ); ...
       acosd( sind(sg(:,2))*sind(ear(4)) + cosd(sg(:,2))*cosd(ear(4)) .* cosd(sg(:,1)-ear(3)) )];

% get unique list of great circle distances and radii
[GCD, ~, gcdID] = unique([gcd repmat(sg(:,3), 2, 1)], 'rows');
% gcd = reshape(GCD(gcdID), size(gcd));
r   = GCD(:,2);
GCD = GCD(:,1);

% get list of frequencies to be calculated
f = 0:fs/Nsamples:fs/2;

% calculate complex the transfer function in the frequency domain
H = sphericalHead_Duda1998( a, r, r_0, GCD/180*pi, f, c, Nsh );

% set 0 Hz bin to 1 (0 dB)
H(1,:) = 1;

% make sure bin at fs/2 is real
if f(end) == 22050
    H(end,:) = abs(H(end,:));
end

% mirror the spectrum
H = AKsingle2bothSidedSpectrum(H, 1-mod(Nsamples, 2));

% get the impuse responses
hUnique = ifft(H, 'symmetric');

% add delay to shift the pulses away from the very start
hUnique = circshift(hUnique, [round(1.5e-3*fs) 0]);

% resort to match the desired sampling grid
h = zeros(Nsamples, size(sg,1), 2);
h(:,:,1) = hUnique(:, gcdID(1:size(sg,1) )    );
h(:,:,2) = hUnique(:, gcdID(size(sg,1)+1:end) );

end

function [ H ] = sphericalHead_Duda1998( a,r, r_0, theta,f,c,Nsh )
% spherical head model according to [1]. Input parameters as in
% AKsphericalHead, with the exception of:
% theta - degree of incidence [rad]

% allocate space for output (1st dimension: freq., 2nd dimension: angle) 
H = zeros(numel(f), numel(theta));

% get unique list of radii
[rUnique, ~, rID] = unique(r);

% normalized distance - Eq. (5) in [1]
rho_0 = r_0     ./ a;
rho   = rUnique ./ a;

% normalized frequency - Eq. (4) in [1]
mu = (2*pi*f*a) / c;


    
% Calculate H
for i = 1:length(theta)
    
    % argument for Legendre polynomial in Eq. (3) in [1]
    x = cos(theta(i));
    
    % initialize the calculation of the Hankel fraction.
    % Appendix A in [1]
    zr = 1./( 1i* mu * rho( rID(i) ) );
    za = 1./(1i * mu);
    Qr2 = zr;
    Qr1 = zr .* (1-zr);
    Qa2 = za;
    Qa1 = za .* (1-za);
    
    % initialize legendre Polynom for order m=0 (P2) and m=1 (P1)
    P2 = 1;
    P1 = x;
    
    % initialize the sum - Eq. (A10) in [1]
    sum = 0;
    
    % calculate the sum for m=0
    term = zr./(za.*(za-1));
    sum = sum + term;
    
    % calculate sum for m=1
    if Nsh > 0
        term = (3 * x * zr .* (zr-1)) ./ (za .* (2*za.^2 - 2*za+1));
        sum = sum + term;
    end
    
    % calculate the sum for 2 <= m <= Nsh
    for m = 2:Nsh
        
        % recursive calculation of the Legendre polynomial of order m
        % (see doc legendreP)
        P = ((2*m-1) * x * P1 - (m-1) * P2) / m;
        
        % recursive calculation of the Hankel fraction
        Qr = - (2*m-1) * zr .* Qr1 + Qr2;
        Qa = - (2*m-1) * za .* Qa1 + Qa2;
        
        % update the sum and recursive terms
        term    = ((2*m+1) * P * Qr) ./ ((m+1) * za .* Qa - Qa1);
        id      = ~isnan(term);         % this might become NaN for high SH orders and low frequencies. However, we usually don't need the high orders for low frequencies anyhow...
        sum(id) = sum(id) + term(id);
        
        Qr2 = Qr1;
        Qr1 = Qr;
        Qa2 = Qa1;
        Qa1 = Qa;
        P2  = P1;
        P1  = P;
    end
    
    % calculate the pressure - Eq. (A10) in [1]
    H(:,i) = (rho_0 * exp( 1j*( mu*rho(rID(i)) - mu*rho_0 - mu) ) .* sum) ./ (1i*mu);
    
end

% [1] uses the Fourier convention with the negative exponent for the
% inverse transform - cf. Eq. (13). Since Matlab uses the opposite
% convention H is conjugated
H = conj(H);

% eof sphericalHead_Duda1998
end