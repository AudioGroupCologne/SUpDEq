function [hrir_data, hrir_angles, errorsig] = enzner2008(x, y, P)
%ENZNER2008  Calculate HRIR set using the method of Enzner et al. (2008)
%   Usage: [hrir_data,hrir_angles,fs] = enzner2008(mu,delta_phi,...)
%
%   Input parameters:
%     x : reference signal (excitation signal used for the measurements)
%
%     y : recorded binaural signal
%
%     P : structure with parameters
%
%       .mu          : NLMS stepzize, e.g., 0.75 or 1
%
%       .delta_phi   : azimuthal resolution (delta_phi) in degree to store
%                      hrir,e.g., 0.1 or 1
%
%       .h_length : length of the impulse responses in samples, e.g., 256 (or 308 for single channel perfect sweeps)
%
%       .adapt : overhead at the end and the beginnig, depends on the recording. 
%          No. of symmetrically overlapping samples of the recorded ear signals
%          arround phi = 180 deg, used as adaptation buffer before HRIR data will
%          be stored, also used to shift the algorithm's input signals to enshure
%          causality, e.g, .adapt = 20000 for the given examples
%
%       .sys_latency : system latency, No. of samples to shift the input signals against
%          each other (ensure causality), e.g. 30 if using a reference
%          recording or -290 if using loudspeaker driving signals(playback signals)
%          whereupon the loudspeaker distance is approx. 2 m (fs = 44100)
%
%   Output parameters:
%     hrir_data : sampled HRIR data at the azimuth-resolution delta_phi with the structure: 
%
%       hrir_data(filter coefficients, left/right, no. of channels, azimuthal index)
%
%       filter coefficients: see h_length
%
%       left/right: 1 = left, 2 = right
%
%       no. of channels: allways 1 (single channel NLMS-algorithm)
%
%       azimuthal index: corresponding to an azimuth phi. The rotational direction during the
%         recording of the ear signals is counterclockwise! (1 = -180 deg, 2 = -180 deg + delta_phi, end = 180 deg - delta_phi )
%
%     hrir_angles : vector with azimuthal angles corresponding to the azimuthal index of hrir_data
%
%     errorsig: error signal (???)
%
%   ENZNER2008 calculates a set of HRIRs using the normalized LMS-algorithm.
%   A test signal in mono, e.g., white noise, perfect sweeps, or a reference
%   recording at the position in the middle of the listeners head is used as
%   the input of the algorithm, whereas the other input of the algorithm is
%   given by the corresponding spatially-continuous (i.e., dynamical) binaural
%   recording.
%
%   This recording contains the measured ear signals along the trajectory of
%   interest, e.g., the horizontal plane, plus some symmetric overhead.  The
%   overhead is used to ensure capturing of all data of interest, to give the
%   algorithm a scope to adapt and to be able to shift the signals against
%   each other to ensure causality (see sys_latency). The binaural recording
%   of the ear signals has to begin/end at the rear of the subject.  Thus the
%   first recorded sample number subsequent to the required overhead (see
%   adapt) corresponds to an azimuth of phi = 180° (rear).
%
%   From a given set of example files, the HRIR data will be calculated. Per
%   default the measured ear signals (stimulus: white noise) and the
%   corresponding reference recording will be used for the computation. If
%   you want to use the loudspeaker driving signals or a perfect sweep data
%   set, please uncomment only the case of interest in the section
%   "changeable parameters" in lines 74-104. In this section you can also
%   adjust the used filter length.
%
%   The computation is performed continuously for each sample, in compliance with
%   a continuous-azimuth HRIR representation. The storage of the HRIR data is
%   sampled with an arbitrary azimuth-spacing delta_phi. The HRIR data will be
%   written into the array hrir_data.
%
%   References:
%     
%     
%     
%
%   See also: exp_enzner2008
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/enzner2008.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%  Authors: Michael Weinert (Michael.Weinert@rub.de), Gerald Enzner (Gerald.Enzner@rub.de)
%  Date: 21-01-2013 (original file)
%  Date: 28-06-2013 (some modifications from the authors)
%  Date: 21-01-2015 (adapted to AMT, Piotr Majdak)


%%% sample numbers at which HRIR will be saved %%%
%%% Note: Computation is continuous corresponding to every sample %%%
circle = 360;                                               % range of measured trajectory in degree
delta_phi_samples = (length(y)-2*P.adapt)/circle * P.delta_phi; % length of ear signals contains full circle plus some overhead (2 * adapt)
if delta_phi_samples < 1
    amt_disp(['delta_phi can not be smaller than ', num2str(circle/(length(y)-2*P.adapt)),' deg = 1 Sample']);
    amt_disp(['delta_phi will be set to ', num2str(circle/(length(y)-2*P.adapt)),' deg = 1 Sample']);
    delta_phi_samples = 1;
end;
save_data_samples = zeros(round((length(y)-2*P.adapt)/delta_phi_samples),1);
for k = 1:length(save_data_samples)
    save_data_samples(k) = round((k-1)*delta_phi_samples)+1;
end;

%%% run NLMS-algorithm %%%
hrir_data = zeros(P.h_length,2,1,length(save_data_samples));
hrir_angles = zeros(length(save_data_samples),1);
h0 = zeros(P.h_length,2);
errorsig = zeros(length(x),2);
n = 1;                                                  % counter
k = P.adapt/2;                                            % counter, calculation begins at sample #adapt/2
                                                        % samples before adapt/2 are reserved for sys_latency
                                                        % samples from adapt/2 to adapt are used as adaptation time
while k <= length(x)-P.sys_latency-P.adapt/2                % note that adapt is also an overhead at the end of measurement
    x_buffer = x(k+P.sys_latency:-1:k+P.sys_latency-P.h_length+1,:);    % input vector backwards, considering sys_latency
    y_estimate(1) = h0(:,1).' * x_buffer(:);            % estimated ear signal left
    y_estimate(2) = h0(:,2).' * x_buffer(:);            % estimated ear signal right
    errorsig(k,1) = y(k,1) - y_estimate(1);                % error signal left
    errorsig(k,2) = y(k,2) - y_estimate(2);                % error signal right
    x_norm = x_buffer' * x_buffer;
    h0(:,1) = h0(:,1) + P.mu/(x_norm) .* errorsig(k,1) .* x_buffer(:); % estimated left ear hrir at sample #k
    h0(:,2) = h0(:,2) + P.mu/(x_norm) .* errorsig(k,2) .* x_buffer(:); % estimated right ear hrir at sample #k

    if n <= length(save_data_samples)
        if (k-P.adapt == save_data_samples(n))            % store hrir at specified azimuth
            hrir_data(:,1,n) = h0(:,1);                 % left
            hrir_data(:,2,n) = h0(:,2);                 % right
            hrir_angles(n) = ((k-1-P.adapt)/delta_phi_samples*P.delta_phi)-180; % corresponding azimuth to stored hrir_data
            n = n+1;
        end;
    end;

k = k+1;
end;

hrir_data = hrir_data./(max(max(max(max(abs(hrir_data))))));





