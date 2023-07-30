function [localization_error,perceived_direction,desired_direction,x,y,x0] = ...
        wierstorf2013(X,Y,phi,xs,src,L,varargin);
%WIERSTORF2013 Sound localization in wavefield synthesis
%   Usage: [...] = wierstorf2013(X,Y,phi,xs,src,L,method,...);
%
%   Input parameters:
%       X           : range of the x-axis [xmin,xmax] or a single point x / m
%       Y           : range of the y-axis [ymin ymax] or a single point y / m
%       phi         : orientation of the listener in rad (0 is in the direction
%                     of the x-axis)
%       xs          : position of the point source in m / direction of the 
%                     plane wave
%       src         : source type
%                     ps for a point source
%                     pw for a plane wave
%       L           : length/diameter of the loudspeaker array
%       method      : reproduction setup
%                     wfs for wave field synthesis
%                     setreo for stereophony
%
%   Output parameters:
%       localization_error  : deviation from the desired direction, defined as
%                             perceived_direction - desired_direction / degree
%       perceived_direction : the direction of arrival the binaural model has
%                             estimated for our given setup / degree
%       desired_direction   : the desired direction of arrival indicated by the
%                             source position xs / degree
%       x                   : corresponding x-axis
%       y                   : corresponding y-axis
%       x0                  : position and directions of the loudspeakers in the
%                             form n x 7, where n is the number of loudspeakers
%
%   WIERSTORF2013(X,Y,phi,xs,src,L,method) calculates the localization
%   error for the defined wave field synthesis or stereophony setup. The
%   localization error is defined here as the difference between the perceived
%   direction as predicted by the dietz2011 binaural model and the desired
%   direction given by xs. The loudspeaker setup for the desired reproduction
%   method is simulated via HRTFs which are than convolved with white noise
%   which is fed into the binaural model.
%
%   The following parameters may be passed at the end of the line of
%   input arguments:
%
%     'resolution',resolution
%                      Resolution of the points in the listening
%                      area. Number of points is resoluation x resolution. If
%                      only one point in the listening area is given via single
%                      values for X and Y, the resolution is automatically set
%                      to 1.
%
%     'nls',nls        Number of loudspeaker of your WFS setup.
%                      Default value is 2.
%
%     'array',array    Array type to use, could be 'linear' or 'circle'.
%                      Default value is 'linear'.
%
%     'hrtf',hrtf      HRTF database in SOFA format.
%                      Default HRTF set is the 3m one from TU-Berlin measured
%                      with the KEMAR.
%
%     'lookup',lookup  Lookup table to map ITD values to angles. This can be
%                      created by the itd2angle_lookuptable function. Default
%                      value is the lookup table itd2angle_lookuptable.mat that comes with AMT.
%
%
%   For the simulation of the wave field synthesis or stereophony setup this
%   functions depends on the Sound-Field-Synthesis Toolbox, which is available
%   here: http://github.com/sfstoolbox/sfs. It runs under Matlab and Octave. The
%   revision used to generate the figures in the corresponding paper was
%   a8914700a4. This one uses a newer version (>=2.4), but should return very
%   similar results.
%
%   See also: wierstorf2013_estimateazimuth, dietz2011, itd2angle
%   itd2angle_lookuptable data_wierstorf2013
%
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592--605, 2011. [1]http ]
%     
%     H. Wierstorf, M. Geier, A. Raake, and S. Spors. A free database of
%     head-related impulse response measurements in the horizontal plane with
%     multiple distances. In Proceedings of the 130th Convention of the Audio
%     Engineering Society, 2011.
%     
%     H. Wierstorf, A. Raake, and S. Spors. Binaural assessment of
%     multi-channel reproduction. In J. Blauert, editor, The technology of
%     binaural listening, chapter 10. Springer, Berlin--Heidelberg--New York
%     NY, 2013.
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S016763931000097X
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/wierstorf2013.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA SFS
%   #Author: Hagen Wierstorf

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.



%% ===== Checking of input parameters and dependencies ===================
nargmin = 7;
nargmax = 17;
narginchk(nargmin,nargmax);
if length(xs)==2, xs = [xs 0]; end

definput.flags.method = {'stereo','wfs'};
definput.keyvals.array = 'linear';
definput.keyvals.nls = 2;
definput.keyvals.resolution = 21;
definput.keyvals.hrtf = [];
definput.keyvals.lookup = [];
definput.keyvals.showprogress = 1;
[flags,kv] = ...
    ltfatarghelper({'resolution','nls','array','hrtf','lookup','showprogress'},definput,varargin);
array = kv.array;
resolution = kv.resolution;
number_of_speakers = kv.nls;
hrtf = kv.hrtf;
lookup = kv.lookup;
showprogress = kv.showprogress;


% Get SFS version
sfsver = strsplit(SFS_version,'.');
% Checking for the Sound-Field-Synthesis Toolbox
if ~exist('SFS_start') | sfsver{1}<2 | (sfsver{1}==2 && sfsver{2}<4)
    error(['%s: you need to install the Sound-Field-Synthesis Toolbox.\n', ...
        'You can download it at https://github.com/sfstoolbox/sfs.\n', ...
        'You need version 2.4.0 or higher of the Toolbox.'], ...
        upper(mfilename));
else
    SFS_start;
end


% Check if we have only one position or if we have a whole listening area
if length(X)==1 && length(Y)==1
    resolution = 1;
end


%% ===== Loading of additional data ======================================
% Load default 3m TU-Berlin KEMAR HRTF, see http://doi.org/10.5281/zenodo.55418
if isempty(hrtf)
    hrtf = amt_load('wierstorf2013', 'qu_kemar_anechoic_3m.sofa');
end
% Get sampling rate from the HRTFs
fs = hrtf.Data.SamplingRate;
% Get lookup table for mapping of ITDs to azimuth angles
if isempty(lookup)
    lookup = amt_load('wierstorf2013','itd2angle_lookuptable.mat');
end


%% ===== Configuration ===================================================
% The following settings are all for the Sound Field Synthesis Toolbox
conf = SFS_config;
conf.N = 4096;
conf.secondary_sources.geometry = array;
conf.secondary_sources.number = number_of_speakers;
conf.secondary_sources.size = L;


%% ===== Simulate the binaural ear signals ===============================
% Get loudspeaker positions
x0 = secondary_source_positions(conf);
% Calculate the stop frequency for the WFS pre-equalization filter
conf.wfs.hprefhigh = aliasing_frequency(x0,conf);

% Selection of loudspeakers for WFS
if flags.do_wfs && strcmpi('circle',array)
    x0 = secondary_source_selection(x0,xs,src);
end

% Get a grid of the listening positions
conf.resolution = resolution; % TODO: is this needed?
x = linspace(X(1),X(end),resolution);
y = linspace(Y(1),Y(end),resolution);

% 700 ms white noise burst
sig_noise = sig_whitenoiseburst(fs);
for ii=1:length(x)
    if showprogress
        amt_disp([num2str(ii) ' of ' num2str(length(x))]);
    end
    for jj=1:length(y)
        X = [x(ii) y(jj) 0];
        if strcmpi('circle',array) && norm(X)>L/2
            desired_direction(ii,jj) = NaN;
            perceived_direction(ii,jj) = NaN;
            localization_error(ii,jj) = NaN;
        else
            desired_direction(ii,jj) = source_direction(X,phi,xs,src);
            if flags.do_stereo
                % Summation of two loudspeakers
                ir1 = ir_point_source(X,phi,x0(1,1:3),hrtf,conf);
                ir2 = ir_point_source(X,phi,x0(2,1:3),hrtf,conf);
                ir = (ir1+ir2)/2;
            else % WFS
                conf.xref = X;
                ir = ir_wfs(X,pi/2,xs,src,hrtf,conf);
            end
            % Convolve impulse response with noise burst
            sig = auralize_ir(ir,sig_noise,1,conf);
            % Estimate the perceived direction.
            % This is done by calculating ITDs with the dietz2011 binaural model,
            % which are then mapped to azimuth values with a lookup table
            perceived_direction(ii,jj) = ...
                wierstorf2013_estimateazimuth(sig, ...
                                              lookup, ...
                                              'dietz2011', ...
                                              'no_spectral_weighting', ...
                                              'remove_outlier');
            localization_error(ii,jj) = ...
                perceived_direction(ii,jj) - desired_direction(ii,jj);
        end
    end
end

end % of main function


%% ----- Subfunctions ----------------------------------------------------
function direction = source_direction(X,phi,xs,src)
    if strcmp('pw',src)
        [direction,~,~] = cart2sph(xs(1),xs(2),xs(3));
        direction = (direction+phi)/pi*180;
    elseif strcmp('ps',src)
        x = xs-X;
        [direction,~,~] = cart2sph(x(1),x(2),x(3));
        direction = (direction-phi)/pi*180;
    end
end


