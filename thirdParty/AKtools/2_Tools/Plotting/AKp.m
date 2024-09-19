% dataOut = AKp(data, plot_type, varargin)
% AKp (short for AKplot) is intended for quick generation of plots that
% look quite ok. If you want you really good plots, use handels dataOut.h
% for manual adjustment.
% 
% See AKplotDemo.m for examples
%
% Only 'data' is required, i.e.
% AKp(data) show a plot of the time signal, energy time signal, magnitude
% response and phase
%
% 'plot_type' can be used to show specific plots. Other arguments are
% passed as attribute value pairs.
%
% Examples:
% AKp(data, 'm2d') .  .  .  .  .  .  .  .  .  . plots spectrum 2D
% AKp(data, 'm3d') .  .  .  .  .  .  .  .  .  . plots spectrum 3D
% AKp(data, 'm4', 'az', az, 'el', el)   .  .  . plots spherical spectrum
% AKp(data, 'itd1', , 'az', az, 'el', el)  .  . plots spherical ITD
%
%
%
% I N P U T:
% (default values are given in brackets. (M) is Matlab default value)
%
% Required input ----------------------------------------------------------
% data      : - time signals of size [S C], or itaAudio object
%               (S = number of samples; C = number of channels)
%             - if plot type is 'itd' or 'ild' data must be a struct
%               containing data.l and data.r with left and right ear hrir's
%
% plot_type : string specifiying the plot followed by '2D', '3D' or a
%             number specifying a spherical plot. (See below)
% 
%   't'   : time curve    .  .  .  .  .  .  .  .  .  .  {2D, 3D}
%   'et'  : energy time curve   .  .  .  .  .  .  .  .  {2D, 3D}
%   'ed'  : energy decay curve  .  .  .  .  .  .  .  .  {2D, 3D}
%   'sr'  : step response .  .  .  .  .  .  .  .  .  .  {2D, 3D}
%   'm'   : magnitude spectrum  .  .  .  .  .  .  .  .  {2D, 3D, spherical}
%   'p'   : phase spectrum   .  .  .  .  .  .  .  .  .  {2D, 3D, spherical}
%   'pu'  : phase spectrum  unwrapped .  .  .  .  .  .  {2D, 3D, spherical}
%   'gd'  : group delay spectrum   .  .  .  .  .  .  .  {2D, 3D}
%   'toa' : time of arrival  .  .  .  .  .  .  .  .  .  {2D, spherical}
%   'itd' : interaural time difference   .  .  .  .  .  {2D, spherical}
%   'ild' : interaural level difference  .  .  .  .  .  {2D, spherical}
%   'x'   : plot data as it is (no fft, no nothing)  .  {spherical}
%           NOTE: in this case data MUST only contain
%           one row but multiple cloumns.
%
% spherical plot types.
% (Use with plot_type, e.g. 'm1' for ballon plot of magnitude spectrum)
%       1 : Balloon, fixed radius, variable color
%       2 : Balloon, variable radius and color
%       3 : Ballon, radius is magnitude, color is phase
%       4 : Balloon, variable radius fixed color
%           (c can be used to pass a color, only [r g b] mode)
%       5 : Planar plot
%       6 : Polar plot (In this case only specifiy az, see below)
%       7 : Spherical Harmonics Plot. Uses radius to decode the
%           magnitude, and the color to decode the sign
%
% If 'plot_type' is a cell or starts with a number, AKpMulti is used
% (e.g. {'ir2d', 'm2d'} or '1a. See AKpMulti for help)
%
%
% --------------------------------------------------------- Universal input
%
% unit          : indicates the physical unit of the input data for plots
%                 and calcualtion the log. values for plots of the
%                 magnitude spectra
%                 - 'Pa'  : (Pascal). In this case the reference for taking
%                            the logarithm is 2e-5 Pascal.
%                 - false : plot data as it is (default)
% N             : truncate or zero-pad data to N samples
%                 (size(data,1) is the default)
% dr            : - two element vector specifiying minimum and maximum 
%                   values to be plotted , or                 
%                 - scalar value specifying the plot range
% du            : data unit
%                 - 'dB' or 'lin' for plots of the magnitude spectra, and
%                    spherical plots wit color coded magnitude (dB)
%                 - 'rad', 'deg' for plots of the phase spectra, and
%                   spherical plots with color coded phase (deg)
%                 - 'n' for samples, 's' for seconds, 'ms' for milli sec., 
%                   'us' for micro seconds when plotting TOA, ITD, or
%                    group delay (us)
% xu            : axis unit
%                 - time axis: 'n' for samples,
%                              's', 'ms', 'us' for (milli/micro) seconds,
%                              'm' for distance in meter according to a
%                                  speed of sound of 343 m/s
%                              ('us' for ITD, 'ms' otherwise)
%                 - frequency-axis : 'Hz', or 'kHz' (Hz)
% norm_d (off)  : normalize data to norm_d (unit according to du)
% fs (44100)    : sampling frequency in Hz
% overlay (on)  : switches hold 'on' or 'off' (see help hold)
% labeling (on) : for default axis labels and figure titles, 'on', 'off'
% f_size (M)    : font size used for labeling plots. Axis font size is (M).
% frac (off)    : fractional octave smoothing for magnitude spectra
%                 eg. frac=3 applies 1/3 octave smoothing
% hv            : pass 'off' to plot something that should not show up in
%                 the legend. Ignore if everything should be available to
%                 the legend. See doc legend.
%
%
% -------------------------------------------------- Arguments for 2D plots
%
% x             : two element vector specifying x-axis limits (i.e. time or
%                 frequency axis). Values according to 'xu'.
%                 NOTE: x differs for toa, itd, ild, plots, see below.
% c             : string or RGB vector specifiying the colors for plotting.
%                 for available color strings call AKcolors(true)
%                 - rgb matrix [C x 3] or string with C elemets
%                   (C = num. of channels) specifies color for each channel
%                 - rgb matrix [2 x 3] or string with 2 elements 
%                   interpolates between first and second entry to obtain
%                   C colors (C = num. of channels)
%                 - single value RGB vector denotes greyscale
%                   (.5 = [.5 .5 .5])
%                 - 'cyc' cycles through the Matlab default colors
%                 (default: 'k'   for single channel data, 'br'  for two 
%                  channel data, 'cyc' otherwise)
% lw (1)        : line width in points
% ls (-)        : line style as string (see doc linespec)
% dash (false)  : 4 element vector containing length of dash1, gap1, dash2,
%                 and gap2 [mm]. If dash is a scalar or 2 element vector,
%                 it will be repeated accordingly.
%                 Only works for 2D plots. Uses dashline.m (Edward Abraham)
% lm (none)     : line marker (see doc linespec)
%
%
% -------------------------------------------------- Arguments for 3D plots
%
% x             : see 2D plots
% hp_view(top_v): Sets view angle for 3D plots. 'side', 'top_h', 'top_v',
%                 or [az el] (see doc view for [az el] mode)
% y             : vector with C elements (C = num. of channels) specifying
%                 the data points (e.g. angles y=-90:90).
%                 Must be monotonic increasing/decreasing.
%                 (1:size(data,2) is the default)
% cm            : (a) string specifying the colormap (cf. doc colormap)
% (RdBu_flip)         append '_flip' to invert the colormap
%                     (e.g. 'gray_flip' wil have white assigned to the
%                      smallest and black assigned to the largest value)
%                 (b) string specifying AKcolormaps (see AKcolormaps.m)
%                     append '_flip' to invert the colormap
%                 (c) colormap as [N x 3] matrix with values between zero
%                     and one, where each row represents a color and N the
%                     number of steps. cr will be obsolete in this case
% cb            : location of colorbar, or false for not showing it.
%                 ('EastOutside', see doc colorbar)
% cr            : scalar value specifying the colormap resolution with unit
%                 according to 'du'. The minimum value defined by 'dr'
%                 will be adjusted if cr does not perfectly fit the range.
%                 By default, a colormap with 128 steps is created.
% ct            : title of colorbar, pass false to omitt the title
% cf (false)    : freezes the colorbar and colormap to allow different
%                 setttings in different subplots. Warning - experimental!
% plot_func     : 'surf' uses Matlabs surf function for plotting (default
%                  if hp_view is [az el])
%                 'img' uses imagesc (default if hp_view is a string)
%                 
%
% ------------------------------------------- Arguments for spherical plots
%
% - x, y and z axis are marked by a red, green and blue line
% - positive x,y and z are markerd by litlle circles
% - x,y and z axis are limited by default
%
% az            : specification of azimuth. This has to be either a
%                 - vector with one entry per channel of data, or
%                 - a matrix with same size as data. In this case data is
%                   plotted as it is
%                 If az contains values > 2pi, it is assumed that its unit
%                 is degree (if not rad are assumed). See 'coord' parameter
%                 below for the used coordinate convention.
% el            : specification of elevation. Same as 'az'. For polar 
%                 plots, only specify az.
% g             : alternative specification of azimuth and eleveation in a
%                 [N x 2] matrix, where the first column holds the azimuth
%                 and the second the elevation.
% sph_f (1000)  : Frequency in Hz to be plotted
% sph_proc      : Data display for balloon plots
%                 - 'interp' : for balloon and planar plots, data is 
%                (default)     interpoalted/extrapolted to a full sphere
%                              with 1 degree resolution.
%                              for polar plots, data is interpolated to
%                              full circle if range(az) > 180 degree.
%                              NOTE: This simple interpolation distortes
%                              the data close to the poles.
%                 - 'interpSplineN' : uses spherical splines of order N for
%                              interpolation. N can be 1, 2, or 3, e.g.
%                              'interpSpline1' uses first order splines
%                 - 'tri'    : data is displayed after triangularization
%                 - 'none'   : can be used to plot data from a meshgrid,
%                              i.e. data, az, and el are [m x n] matrices
% triPop (false): excludes oddly shaped triangles from the plot if 
%                 sph_proc = 'tri'. This is for example helpfull if your
%                 data is not full-spherical.
%                 triPop = 1 will work in most cases. Larger values will
%                 discard more triangles, smaller values less.
% coord (1)     : specify the used coordinate convention. This also defines
%                 the orientation of x,y and z axis. New coordinate systems
%                 can be implemented in AKpCoordinateTransform.m.
%                 Implemente conventions are described in the same file.
%                 1: Matlab default
%                    (default, see doc sph2cart.)
%                 2: Mathematical
%                    (az 0=front, 90=left; el 0=North Pole, 90=front)
%                 NOTE: units of ccordinates are automatically set
%                 (see parameter 'az')
% axis_l (1)    : line width of axis
% axis_m (20)   : marker size for marking posititve x,y and z in spherical
%                 plots
% axis_s        : style of the axis
%                 0: do not draw axis
%                 1: x, y, and z axis are marked by red, green, and blue
%                    lines. Positive x, y, and z direction is marked by dot
%                 2: like 1, but additonal lines mark intermediate angles
%                 (Default is 2 for fixed radius plots, and 1 otherwise)
% cb, cm, cr,   : see 3D plots
% cf
%
% -------------------------------------------- Plot type specific arguments
% toa, itd, ild
%
% x             : vector with C elements (C = num. of channels) specifying
%                 the data points,  e.g. angles [-90:90]
%                 (1:size(data,2))
% toa_us (10)   : upsampling factor for onset detection, (only toa, itd)
% toa_th (6)    : threshold for onset detection in dB (only toa, itd)
% toa_mode (rel): Onset detection mode. Type 'help AKonsetDetect' in the
%                 command window for information
%                 
%
%
% O U T P U T
% dataOut       : processed input data. This is a struct containing the
%                   processed data and
%                 - time or freq axis (2D plots)
%                 - time or freq and y atrix (3D plots)
%                 - azimuth and elevation matrix (balloon and planar plots)
%                 - azimuth axis (polar plots)
%                 - plot handle h
%
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin,
% DFG research unit 'SEACEN', 7/2012
% 
% Thanks for hints in form of code:
% Alexander Lindau, Frank Schultz, Vera Erbes,
% SOFiA Toolbox (Bennjamin Bernschuetz)

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

% PROCESSING
% 0:  Categorization of plot type into
%     - time, freq, other
%     - 2D, 3D, spherical
% 1:  Parsing of varargin
% 2:  Calculation of auxiliary variables
% 3:  data transformation (e.g. fft), data type (e.g. Magnitude in dB),
%     and plot title (e.g. Magnitude spectrum)
% 4:  Fractional octave smoothing
% 5:  Normalization
% 6:  Autorange
% 7:  Process data for spherical plots
% 8:  Define x-axis and x-label
% 8:  Plotting
% 10: Delete output argument if not requested
%
% NOTE: The code within the blocks sometimes is a bit redundant, but that
% makes it easier to understand and maintain.
%
%
% Overload variables    Difference attribute
% x                     2D, 3D, toa, itd, ild
% WARNING process overloaded varibales only within switch or if cases to
% avoid misbehaviour of function.
%
% Hardcoded variables
% min_dB = -300; to avoid clipping of -inf dB
% set(gcf, 'color', [1 1 1])
%
% toDo (ordered according to importance)
% - sph_proc 'none' makes problems...
% - spherical planar plot using world-like maps
% - output format of az and el in data struct is allways matlab default
% - spectogram, waterfall diagram
%   [sc]  : spectogram, single channel   .  .  .  . {3D}
%   [w]   : waterfall diagram , single channel .  . {3D}

function dataOut = AKp(data, plot_type, varargin)
%% ----------------------------------------- 0. categorization of plot type

% overview plot
if nargin == 1
    AKpMulti(data, '1a')
    return
end

% pass to AKpMulti
if ( ~isnan(str2double(plot_type(1))) && isreal(str2double(plot_type(1))) ) || iscell(plot_type)
    dataOut = AKpMulti(data, plot_type, varargin{:});
    
    if ~nargout
        clear dataOut
    end
    
    return
end

% ... return with AKp ...

% squeeze input data in case matrix with the dimensions [N,1,C] was given
if ndims(data) == 3 && size(data,2) == 1
    data = squeeze(data);
end

% separate plot_type and plot_look
if ~isnan(str2double(plot_type(end)))
    plot_look = str2double(plot_type(end));
    plot_type = plot_type(1:end-1);
elseif length(plot_type) >= 3
    if strcmpi(plot_type(end-1:end), '2D')
        plot_look = '2D';
        plot_type = plot_type(1:end-2);
    elseif strcmpi(plot_type(end-1:end), '3D')
        plot_look = '3D';
        plot_type = plot_type(1:end-2);
    end
else
    error('AKp:InputPlotType',[plot_type ' is not a valid plot_type. See doc hp for a list of valid plot_types'])
end

% catch deprecated plot looks
if any(strcmpi(plot_type, {'tc' 'ir'}))
    plot_type = 't';
end
if strcmpi(plot_type, 'etc')
    plot_type = 'et';
end
if strcmpi(plot_type, 'edc')
    plot_type = 'ed';
end
if any(strcmpi(plot_type, {'ms' 's'}))
    plot_type = 'm';
end
if any(strcmpi(plot_type, {'ps' 'ph'}))
    plot_type = 'p';
end
if strcmpi(plot_type, 'pus')
    plot_type = 'pu';
end
if any(strcmpi(plot_type, {'g' 'gds'}))
    plot_type = 'gd';
end

% get the plot domain
if any(strcmpi(plot_type, {'t' 'et' 'ed' 'sr'}))
    plot_domain = 'time';
elseif any(strcmpi(plot_type, {'m' 'p' 'pu' 'gd' 'sc' 'w'}))
    plot_domain = 'freq';
elseif any(strcmpi(plot_type, {'toa' 'itd' 'ild' 'x'}))
    plot_domain = plot_type;
else
    error('AKp:InputPlotType',[plot_type ' is not a valid plot_type. See doc hp for a list of valid plot_types'])
end

if any(strcmpi(plot_type, {'itd' 'ild'})) && ...
    (~isfield(data, 'l') || ~isfield(data, 'r'))
    error('AKp:InputLR', 'input variable ''data'' must be a struct with fields ''l'', and ''r''')
end

%% ------------------------------------------ 1a. define default parameters
inputnames = {...
'N' 'fs' 'overlay' 'norm_d' 'hv' 'unit' ...    % universal parameters
'labeling' 'f_size', ...                       % plot layout
'axis_m' 'axis_l' 'axis_s', ...
'xu' 'frac' 'du' ...                           % frequency and time parameters
'hp_view' 'x' 'y' 'dr'...                      % 2D and 3D plots
'c' 'lw' 'ls' 'dash' 'lm'...
'cb' 'cm' 'cr' 'ct' 'cf' 'plot_func'...
'az' 'el' 'g' 'sph_f' 'sph_proc' 'triPop' ...  % spherical plots parameters
'coord' ...
'toa_us' 'toa_th' 'toa_mode'...                % for plots toa, itd, g
'line_w' 'line_s' 'line_m' ...                 % deprecated
};

% Define default values
if ~strcmpi(class(data), 'itaAudio')
    def.N        = size(data,1);
else
    def.N        = size(data.time,1);
end
def.fs       = 44100;
def.overlay  = 'on';
def.hv       = false;
def.unit     = false;
def.norm_d   = false;
def.labeling = 'on';

def.f_size  = get(0,'DefaultAxesFontSize');
def.axis_l  = 1;
def.axis_m  = 10;
if isnumeric(plot_look)
    if plot_look == 1
        def.axis_s  = 2;
    else
        def.axis_s  = 1;
    end
else
    def.axis_s  = 1;
end

if strcmpi(plot_domain, 'time')
    def.xu  = 'ms';
else
    def.xu  = 'Hz';
end
def.frac    = false;

if strcmpi(plot_type, 'itd')
    def.du  = 'us';
elseif any(strcmpi(plot_type, {'p' 'pu'})) || ~ischar(plot_look)
    def.du  = 'deg';
elseif strcmpi(plot_type, 'm')
    def.du  = 'dB';
else
    def.du  = 'us';
end

if strcmpi(plot_look, '3D')
    def.hp_view = 'top_v';
elseif ~ischar(plot_look) && plot_look <= 4
    def.hp_view = [135 15];
else
    def.hp_view = false;
end
def.x       = false;
def.y       = false;
def.dr      = false;

if size(data,2)==1
    def.c = 'k';
elseif size(data,2)==2
    def.c = 'br';
else
    def.c = 'cyc';
end    
def.lw    = 1;
def.ls    = '-';
def.dash  = false;
def.lm    = 'none';

def.cb        = 'EastOutside';
def.cm        = 'RdBu_flip';
def.cr        = false;
def.ct        = [];
def.cf        = false;
def.plot_func = false;

def.az      = false;
def.el      = false;
def.g       = false;
def.sph_f   = 1000;
if plot_look ~= 6
    def.sph_proc= 'interp';
else
    def.sph_proc= 'interp';%'none';
end
def.triPop  = false;
def.coord   = 1;

def.toa_us   = 10;
def.toa_th   = 6;
def.toa_mode = 'rel';

% deprecated
def.line_w = NaN;
def.line_s = NaN;
def.line_m = NaN;



% --------------------------------------------------------- 1b. parse input
% Check the number of arguments to be even
if mod(size(varargin,2),2)
    error('AKp:InputCount','Odd number of input arguments. Check attribute-value pairs.')
end

% Stock the values
for n=1:2:size(varargin,2)-1
    % checks for all possible attribute/value pairs
    is_parameter = 0;
    for m=1:size(inputnames,2) 
        % detect submitted attribute/value pairs
        if strcmp(inputnames{m},varargin{n})
            % create input-variables from submitted attribute-value-pairs (supress output)
            [~] = evalc([inputnames{m},'=','varargin{n+1}']);
            is_parameter = 1;
        end
    end
    if ~is_parameter
        error('AKp:InputParameter',['No such parameter: ' varargin{n}])
    end
end

% Create missing input-variables with default-values
for m=1:size(inputnames,2)
    default = ['def.' inputnames{m}];
    % if input-variable m hasn't already been defined...
    if ~exist(inputnames{m},'var')
        [~] = evalc([inputnames{m},'=',default]); % (supress output)
    end
end

% check if old input parameters were passed and use them to overwrite the
% correspondingn default values of the new input parameters
if ~isnan(line_w); lw = line_w; end
if ~isnan(line_s); ls = line_s; end
if ~isnan(line_m); lm = line_m; end

clear line_w line_s line_m

% check size of dash
dash = reshape(dash, [1 numel(dash)]);
dash = repmat(dash, [1 4/numel(dash)]);

% delete unnecassary variables to have a clean workspace inside the
% function (better for debugging)
clearvars def default in_args inputnames m n ninput trash varargin


%% ------------------------------------------------ 2.  auxiliary variables
% check input data type
if strcmpi(class(data), 'itaAudio')
    fs   = data.samplingRate;
    data = data.timeData;
end

% Zero pad or truncate
if ~isstruct(data)
    if size(data,1)>N        %#ok<*NODEF>
        data = data(1:N,:);
    else
        data(end+1:N,:) = 0;
    end
else
    if size(data,1)>N
        data.l = data.l(1:N,:);
        data.r = data.r(1:N,:);
    else
        data.l(end+1:N,:) = 0;
        data.r(end+1:N,:) = 0;
    end
end

% copy original input
data_cp = data;

% number of samples and channels
if ~isstruct(data)
    [N,C] = size(data);
else
    [N,C] = size(data.l);
end

% frequency vector
f     = (0:N-1)' * fs/N;

% minimum to be displayed in 3D plots if data has inf entries
min_dB = -300;

% define plot function for 3D plots
if ~plot_func
    if strfind(hp_view, 'top') %#ok<*STRIFCND>
        plot_func = 'img';
    else
        plot_func = 'surf';
    end
else
    if strfind(hp_view, 'side')
        plot_func = 'surf';
    end
end

% check data format of az for spherical plots
if ~islogical(g)
    % check dimension
    if size(g,1) ~= 2 && size(g,2) ~= 2
        error('AKp:Input', '''g'' must be of size [N x 2]')
    end
    if size(g, 2) ~= 2 && size(g, 1) == 2
        g = g';
    end
    
    az = g(:,1);
    el = g(:,2);
    
end

if ~islogical(az)
    if size(az,1) == size(data,1) && size(az,2) == size(data,2)
        sph_proc = 'none';
        % disable processing of the data as done in next step
        plot_type_cp = plot_type;
        plot_type = 'spherical no transformation';
    else
        % transpose azimuth and elevation if necessary
        if size(az,1) < size(az,2)
            az = az';
        end
        if exist('el', 'var') && size(el,1) < size(el,2)
            el = el';
        end
    end
end

% get nice colors from Matlab defaults and add some more
colors = AKcolors;


%% ----------------------- 3. data transformation, data type and plot title

if strcmpi(unit, 'pa')
    dBref = 2e-5;
else
    dBref = 1;
end

% --------------------------------------------------------- t (time signal)
if strcmpi(plot_type, 't')
    d_type = 'Amplitude';
    t_str  = 'Time signal';
end
% -------------------------------------------------- et (energy time curve)
if strcmpi(plot_type, 'et')
    data   = 20*log10(abs(data)/dBref);
    data( data == -inf ) = db(eps);
    d_type = 'Amplitude in dB';
    t_str  = 'Energy Time Curve (ETC)';
end
% ------------------------------------------------- ed (energy decay curve)
if strcmpi(plot_type, 'ed')
    data   = AKedc(data);
    d_type = 'EDC in dB';
    t_str  = 'Energy Decay Curve (EDC)';
end
% ------------------------------------------------------ sr (step response)
if strcmpi(plot_type, 'sr')
    for k = 1:C
        data(:,k) = stepz(data(:,k), 1, N, fs);
    end
    clear k
    d_type = 'Amplitude';
    t_str  = 'Step response';
end
% -------------------------------------------------- m (magnitude spectrum)
if strcmpi(plot_type, 'm')
    if strcmpi(du, 'lin')
        data = abs(fft(data));
        d_type = 'Amplitude';
        t_str  = 'Magnitude spectrum';
    else
        data = 20*log10(abs(fft(data/dBref)));
        data( data == -inf ) = db(eps);
        d_type = 'Amplitude in dB';
        t_str  = 'Magnitude spectrum';   
    end
end
% -------------------------------------------------------- gd (group delay)
if strcmpi(plot_type, 'gd')
    for k = 1:C
        data(:,k)   = grpdelay(data(:,k), 1, f, fs);
    end
    clear k
    t_str  = 'Group delay';
    
    if strcmpi(du, 'n')
        d_type = 'group delay in samples';
    elseif strcmpi(du, 's')
        data   = data/fs;
        d_type = 'group delay in s';
    elseif strcmpi(du, 'ms')
        data   = data/fs*10^3;
        d_type = 'group delay in ms';
    elseif strcmpi(du, 'us')
        data   = data/fs*10^6;
        d_type = 'group delay in \mus';
    else
        error('AKp:InputDataUnit',[du ' is not a valid data unit. See doc hp for a list of valid d_units'])
    end
end
% -------------------------------------------- p (phase spectrum - wrapped)
if strcmpi(plot_type, 'p')
    data   = angle(fft(data));
    t_str  = 'Phase spectrum';
    
    if strcmpi(du, 'rad')
        d_type = 'phase in radians';
    elseif strcmpi(du, 'deg')
        data   = rad2deg(data);
        d_type = 'phase in degree';
    else
        error('AKp:InputDataUnit',[du ' is not a valid data unit. See doc hp for a list of valid d_units'])
    end
end
% ----------------------------------------- pu (phase spectrum - unwrapped)
if strcmpi(plot_type, 'pu')
    data   = unwrap(angle(fft(data)));
    t_str  = 'Unwrapped phase spectrum';
    
    if strcmpi(du, 'rad')
        d_type = 'phase in radians';
    elseif strcmpi(du, 'deg')
        data   = rad2deg(data);
        d_type = 'phase in degree';
    else
        error('AKp:InputDataUnit',[du ' is not a valid data unit. See doc hp for a list of valid d_units'])
    end
end
% --------------------------------------------------- toa (time of arrival)
if strcmpi(plot_type, 'toa')
    data   = AKonsetDetect(data, toa_us, toa_th, toa_mode);
    d_type = 'TOA in samples';
    t_str  = 'Time of Arrival (TOA)';
end
% ---------------------------------------- itd (interaural time difference)
if strcmpi(plot_type, 'itd')
    tmp.l  = AKonsetDetect(data.l, toa_us, toa_th, toa_mode);
    tmp.r  = AKonsetDetect(data.r, toa_us, toa_th, toa_mode);
    data   = tmp.l - tmp.r;
    d_type = 'ITD in samples';
    t_str  = 'Interaural Time Difference (ITD)';
    clear tmp
end

if max(strcmpi(plot_type, {'toa' 'itd'}))
    if strcmpi(du, 'n')
        d_type = [upper(plot_type) 'in samples'];
    elseif strcmpi(du, 's')
        data   = data/fs;
        d_type = [upper(plot_type) 'in s'];
    elseif strcmpi(du, 'ms')
        data   = data/fs*10^3;
        d_type = [upper(plot_type) 'in ms'];
    elseif strcmpi(du, 'us')
        data   = data/fs*10^6;
        d_type = [upper(plot_type) 'in \mus'];
    end
end
% --------------------------------------- ild (interaural level difference)
if strcmpi(plot_type, 'ild')
    % same calculation as in sound field synthesis toolbox, but sign is
    % switched
    data   = db(rms(data.l)) - db(rms(data.r));
    d_type = 'ILD in dB';
    t_str  = 'Interaural Level Difference (ILD)';
    clear tmp
end


if strcmpi(plot_type, 'spherical no transformation')
    plot_type = plot_type_cp;
    d_type = '';
    clear plot_type_cp
end

if strcmpi(plot_type, 'x')
    d_type = '';
    t_str  = '';
end

% add unit to data type
if strcmpi(d_type, 'amplitude')
    if strcmpi(unit, 'pa')
        d_type = [d_type ' in Pascal'];
    end
end

if strfind(d_type, 'in dB')
    if strcmpi(unit, 'pa')
        d_type = [d_type ' SPL'];
    end
end

%% ----------------------------------------- 4. fractional octave smoothing
% (discard frequency bin at fs/2 if N is even)
if ~islogical(frac) && (strcmpi(plot_type, 'm'))
    data = AKfractOctSmooth(data(1:ceil(N/2),:), 'welti', fs, frac);
    f    = f(1:ceil(N/2));
end


%% ------------------------------------------------------- 5. normalization
if ~islogical(norm_d)
    if strfind(upper(d_type), 'DB')
        data = data -max(max(data)) + norm_d;
    else
        data = data / max(max(abs(data))) * norm_d;
    end
end

%% ----------------------------------------------------------- 6. autorange
if islogical(dr)
    % find maximum and set visible plot range to next highest values that
    % is a divisor of 10
    tmp_max = max(max(max(data)));
    tmp_max = tmp_max - mod(tmp_max, 10)+10;
    if ceil(max(max(max(data)))) == tmp_max
        tmp_max = tmp_max + 10;
    end
    
    % set minimum value according to plot_type
    if strcmpi(plot_type, 'et') || strcmpi(plot_type, 'ed')
        dr = [tmp_max-100 tmp_max];
    elseif strcmpi(plot_type, 'm')
        dr = [tmp_max-60 tmp_max];
    end
elseif numel(dr) == 1
    % find maximum and set visible plot range to next highest values that
    % is a divisor of 10
    tmp_max = max(max(max(data)));
    tmp_max = tmp_max - mod(tmp_max, 10)+10;
    if ceil(max(max(max(data)))) == tmp_max
        tmp_max = tmp_max + 10;
    end
    
    if any(strcmpi(plot_type, {'et' 'ed' 'm'}))
        dr = [tmp_max-dr tmp_max];
    end
elseif numel(dr) == 2
    dr = sort(dr);
end

clear tmp_max

%% ------------------------------------ 7. process data for spherical plots
if ~ischar(plot_look) && any(strcmpi(plot_type, {'m' 'p' 'pu' 'toa' 'itd' 'ild' 'x'}))
    
    % index of frequency to be plotted
    fID    = round(sph_f / (fs/N)) + 1;
    
    if (strcmpi(sph_proc, 'none') && plot_look ~= 6) || ...
        strcmpi(plot_type, 'x')     
        data_c = data;
    else
    
        % get magnitude and phase spectrum or toa, itd, ild and assign to
        % new variables for further processing
        if strcmpi(plot_type, 'm')
            data_abs = data;
            data_ang = angle(fft(data_cp));
        elseif any(strcmpi(plot_type, {'p' 'pu'}))
            data_abs = 20*log10(abs(fft(data_cp/dBref)));
            data_abs( data_abs == -inf ) = db(eps);
            data_ang = data;
        elseif any(strcmpi(plot_type, {'toa', 'itd', 'ild'}))
            if plot_look ~= 1
                error('AKp:InputPlotType',['plot_type ' plot_type ' only supports ' plot_type '1'])
            end
            data_abs = data;
            data_ang = data;
        end
        
        % get data to be plotted (In ballon and polar plots this determines
        % the in planar plots it determines the color)
        if plot_look == 1
            data = ones(size(data));
        elseif strcmpi(plot_type, 'm') || plot_look == 3 || ...
               max(strcmpi(plot_type, {'toa', 'itd', 'ild'}))
            data   = data_abs;
            d_type = 'amplitude in dB';
        else
            data = data_ang;
            if strcmpi(du, 'rad')
                d_type = 'phase in radians';
                if strcmpi(plot_type, 'pu')
                    if ~isnumeric(dr)
                        dr = [min(min(data(fID,:))) max(max(data(fID,:)))];
                        if dr(1) == dr(2)
                            dr(1) = dr(1)-1;
                            dr(2) = dr(2)+1;
                        end
                    end
                else
                    dr = [-pi pi];
                end
            elseif strcmpi(du, 'deg')
                if strcmpi(plot_type, 'm')
                    data   = rad2deg(data);
                end
                d_type = 'phase in degree';
                if strcmpi(plot_type, 'pu')
                    if ~isnumeric(dr)
                        dr = [min(min(data(fID,:))) max(max(data(fID,:)))];
                        if dr(1) == dr(2)
                            dr(1) = dr(1)-1;
                            dr(2) = dr(2)+1;
                        end
                    end
                else
                    dr = [-180 180];
                end
            else
                error('AKp:InputDataUnit',[du ' is not a valid data unit. See doc hp for a list of valid d_units'])
            end
        end
        
        % get data that determines the color to be plotted.
        % auto set data range for color, if the phase determines the color.
        if (max(plot_look == [1 2 4 5 6]) && strcmpi(plot_type, 'm')) || ...
            max(strcmpi(plot_type, {'toa', 'itd', 'ild'}))
            data_c = data_abs;
        else
            data_c = data_ang;
            if strcmpi(du, 'rad')
                d_type = 'phase in radians';
                if strcmpi(plot_type, 'pu')
                    if ~isnumeric(dr)
                        dr = [min(min(data_c(fID,:))) max(max(data_c(fID,:)))];
                        if dr(1) == dr(2)
                            dr(1) = dr(1)-1;
                            dr(2) = dr(2)+1;
                        end
                    end
                else
                    dr = [-pi pi];
                end
            elseif strcmpi(du, 'deg')
                if strcmpi(plot_type, 'm')
                    data_c   = rad2deg(data_c);
                end
                d_type = 'phase in degree';
                if strcmpi(plot_type, 'pu')
                    if ~isnumeric(dr)
                        dr = [min(min(data_c(fID,:))) max(max(data_c(fID,:)))];
                        if dr(1) == dr(2)
                            dr(1) = dr(1)-1;
                            dr(2) = dr(2)+1;
                        end
                    end
                else
                    dr = [-180 180];
                end
            else
                error('AKp:InputDataUnit',[du ' is not a valid data unit. See doc hp for a list of valid d_units'])
            end 
        end
        
        % pick the frequency to be plotted - or use plot type for title
        if any(strcmpi(plot_type, {'m' 'p' 'pu'}))
            fID    = round(sph_f / (fs/N)) + 1;
            data   = data(fID,:);
            data_c = data_c(fID,:);
            sph_f  = f(fID);
        else
            sph_f = upper(plot_type);
        end
        
        clear data_abs data_ang fID
    end
    
    % coordinate conversion from user to matlab default format
    co = AKpCoordinateTransform(az, el, coord);
    
    % clip data to dr if its information is used as radius in balloon or
    % polar plots
    if ~islogical(dr) && max(plot_look == [2:4 6])
        data = max(data, dr(1));
        data = min(data, dr(2));
        if exist('data_c', 'var')
            data_c = max(data_c, dr(1));
            data_c = min(data_c, dr(2));
        end
    % clip data that is used for coloration to min_dB. Otherwise the
    % colorbar can not be constructed, if data_c contains -inf values.
    else
        data(isinf(data)) = min_dB;
        if exist('data_c', 'var')
            data_c(isinf(data_c)) = min_dB;
        end
    end
    % if data is as radius in balloon plots, negative values don't make any
    % sense. Everthing has to be shiftet to positive values
    if min(min(data)) <= 0 && max(plot_look == 2:4)
        data = data - min(min(data)) + 1;
    end

    % processing for polar plots
    if plot_look == 6
        % data has to be sorted
        [co.az, id] = sort(co.az);
        data = data(id);
        % data has to be clipped for polar plots
        if ~islogical(dr)
            data = max(data, dr(1));
            data = min(data, dr(2));
        else
            dr = [min(data) max(data)];
        end
    end
    
    % interpolation and triangularization
    if plot_look == 6
    % polar plots
        if strcmpi(sph_proc, 'interp')
            if range(co.az) > pi
            % interpolate to full circle
                % copy smalles value to 0 and 360 degree to avoid
                % extrapolation which would result in a discontinuity at 0
                [tmp_min, tmp_id] = min(co.az);
                if tmp_min ~= 0
                    co.az = [0; co.az];
                    data  = [data(tmp_id) data];
                elseif max(co.az) < 2*pi
                    co.az = [co.az; 2*pi];
                    data  = [data data(tmp_id)];
                end
                % interpolate
                data = interp1(co.az, data, (0:360)/180*pi);
                co.az = (0:360)'/180*pi;
                
                clear tmp_min tmp_id
            else
            % interpolate to half circle
                % construct new azimuth
                tmp = min(co.az):pi/180:max(co.az);
                % check if the last value of co.az is included in the new
                % azimuth vector
                if co.az(end) ~= tmp(end)
                    tmp = [tmp co.az(end)];
                end
                % interpolate
                data = interp1(co.az, data, tmp);
                co.az = tmp;
            end
        elseif strcmpi(sph_proc, 'none') 
        elseif strcmpi(sph_proc, 'tri') 
            error('AKp:InputShpericalProcessing','Triangularization is not implemeted for polar plots')
        elseif strfind('interpspline', lower(sph_proc))
            error('AKp:InputShpericalProcessing','Spherical spline interpolation is not implemeted for polar plots')
        else
            error('AKp:InputShpericalProcessing',[sph_proc ' is not a valid shperical processing. See doc hp.'])
        end
    else
    % planar and balloon plots
        if strcmpi(sph_proc, 'interp')
            % this warning does not affect us
            warning('off', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
            
            % get radius data
            if plot_look == 1
                % ballon with radius 1
                [az_tmp, el_tmp] = meshgrid((0:360)*pi/180, (90:-1:-90)*pi/180);
                [X, Y, Z] = sph2cart(az_tmp, el_tmp, ones(181,361));
                clear az_tmp el_tmp
            elseif plot_look == 7
                % absolute value (i.e. for spherical harmonics)
                [X, Y, Z, data] = AKpDataInterpolation(co, abs(data));
            else
                % radius of ballon
                [X, Y, Z, data] = AKpDataInterpolation(co, data);
            end
            
            % get color data
            [~,~,~, data_c] = AKpDataInterpolation(co, data_c);
            
            warning('on', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
            
        elseif strfind(lower(sph_proc), 'interpspline')
            
            % get spline order
            if numel(sph_proc) > 12
                mSpline = str2double(sph_proc(end));
            else
                mSpline  = 1;
                sph_proc = [sph_proc '1']; %#ok<NASGU>
            end
            
            % get 1x1 grid
            [az_tmp, el_tmp] = meshgrid((0:360)/180*pi, (90:-1:-90)/180*pi);
            AZ_tmp = reshape(az_tmp, numel(az_tmp), 1);
            EL_tmp = reshape(el_tmp, numel(el_tmp), 1);
            
            
            % get radius data
            if plot_look == 1
                % ballon with radius 1
                [az_tmp, el_tmp] = meshgrid((0:360)*pi/180, (90:-1:-90)*pi/180);
                [X, Y, Z] = sph2cart(az_tmp, el_tmp, ones(181,361));
            elseif plot_look == 7
                % absolute value (i.e. for spherical harmonics)
                data = AKsphSplineInterp(co.az, co.el, reshape(abs(data), numel(co.az), 1), AZ_tmp, EL_tmp, mSpline, 0, 'rad', false);
                data = reshape(data, 181, 361);
                [X, Y, Z] = sph2cart(az_tmp, el_tmp, data);
            else
                % radius of ballon
                data = AKsphSplineInterp(co.az, co.el, reshape(data, numel(co.az), 1), AZ_tmp, EL_tmp, mSpline, 0, 'rad', false);
                data = reshape(data, 181, 361);
                [X, Y, Z] = sph2cart(az_tmp, el_tmp, data);
            end
            
            % get color data
            data_c = AKsphSplineInterp(co.az, co.el, reshape(data_c, numel(co.az), 1), AZ_tmp, EL_tmp, mSpline, 0, 'rad', false);
            data_c = reshape(data_c, 181, 361);
            
            % treat remainder like interpolated data
            sph_proc = 'interp';
            
            clear az_tmp el_tmp AZ_tmp EL_tmp mSpline
            
        elseif strcmpi(sph_proc, 'tri')
            % Balloon plots
            if plot_look ~= 5
                % create trangularization with unity radius
                [X, Y, Z] = sph2cart(co.az, co.el, ones(size(co.az)));
                % tri = DelaunayTri([X Y Z]);
                tri = delaunayTriangulation([X Y Z]);
                ch = convexHull(tri);
                % obtain x, y, z coordinates of vertices with correct
                % radius for plotting
                if plot_look == 7
                    % absolute value (i.e. for spherical harmonics)
                    [X, Y, Z] = sph2cart(co.az, co.el, abs(data'));
                else
                    [X, Y, Z] = sph2cart(co.az, co.el, data');
                end
                
            % Planar plots
            else
                 % create 2D trangularization
                 % tri = DelaunayTri(az, el);
                 tri = delaunayTriangulation(az, el);
            end
        elseif strcmpi(sph_proc, 'none')
            data = reshape(data, [numel(co.az) 1]);
            if plot_look == 1
                [X, Y, Z] = sph2cart(co.az, co.el, ones(size(data)));
            elseif plot_look == 7
                % absolute value (i.e. for spherical harmonics)
                [X, Y, Z] = sph2cart(co.az, co.el, abs(data));
            else
                [X, Y, Z] = sph2cart(co.az, co.el, data);
            end
        else
            error('AKp:InputShpericalProcessing',[sph_proc ' is not a valid argument. See doc hp'])
        end
    end    
    
elseif ~ischar(plot_look) && ~max(strcmpi(plot_type, {'m' 'p' 'pu' 'toa' 'itd' 'ild', 'x'}))
    error('AKp:InputShpericalProcessing',['spherical plots are not possible for plot_type ' plot_type])
end


%% ------------------------------------------------- 8. x-axis and x-labels
% check for consitency
if strcmpi(plot_domain, 'time') && ~any( strcmpi(xu, {'s' 'ms' 'us' 'm' 'n'}) )
    xu = 'ms';
end
if strcmpi(plot_domain, 'freq') && ~any( strcmpi(xu, {'hz' 'kHz'}) )
    xu = 'Hz';
end

% convert (from samples to sec, etc.)
if strcmpi(xu, 's')
    n = (0:N-1) * 1/fs;
    x_label = 'Time in s';
elseif strcmpi(xu, 'ms')
    n =  (0:N-1) * 1/fs * 10^3;
    x_label = 'Time in ms';
elseif strcmpi(xu, 'us')
    n = (0:N-1) * 1/fs * 10^6;
    x_label = 'Time in \mus';
elseif strcmpi(xu, 'n')
    n = 1:N;
    x_label = 'Time in samples';
elseif strcmpi(xu, 'm')
    n = (0:N-1) * 343/fs;
    x_label = 'Distance of flight in m';
elseif strcmpi(xu, 'hz')
    n = f;
    x_label = 'Frequency in Hz';
    % default xtick and label for frequency plots
    if fs <= 48e3
        [f_tick, f_label]  = AKfrequencyTicks('Hz', false, 20e3);
    elseif fs <= 100e3
        [f_tick, f_label]  = AKfrequencyTicks('Hz', false, 40e3);
    else
        [f_tick, f_label]  = AKfrequencyTicks('Hz', false, 90e3);
    end
elseif strcmpi(xu, 'khz')
    n = f;
    x_label = 'Frequency in kHz';
    % frequency axis   
    % default xtick and label for frequency plots
    if fs <= 48e3
        [f_tick, f_label]  = AKfrequencyTicks('kHz', false, 20e3);
    elseif fs <= 100e3
        [f_tick, f_label]  = AKfrequencyTicks('kHz', false, 40e3);
    else
        [f_tick, f_label]  = AKfrequencyTicks('kHz', false, 90e3);
    end
else
    error('AKp:InputAxisUnit',[xu ' is not a valid x axis unit. See doc hp for a list of valid x_units'])
end

% define default range and xlabels for toa, itd and ild
if max(strcmpi(plot_domain, {'toa' 'itd' 'ild'}))
    % range
    if islogical(x)
        n = 1:C;
    else
        n = x;
        x = [n(1) n(end)];
    end
    % xlabel is somehow a guess
    if strcmpi(plot_domain, {'toa'})
        x_label = 'data points';
    else
        x_label = 'source/head orientation';
    end
end



%% ----------------------------------------- 9. plot preparations and plots

warning('off', 'MATLAB:warn_r14_stucture_assignment');

% check x limits
if islogical(x)
    x = [min(n) max(n)];
    if strcmpi(plot_domain, 'freq')
        x(1) = max(n(2), 20);
        x(2) = fs/2;
    end
end

% convert color strings to RGB values
if ~strcmpi(c, 'cyc')
    if ischar(c)
        c_str = c;
        c     = zeros(numel(c), 3);
        
        for nn = 1:numel(c_str)
            c(nn,:) = colors.rgb(colors.string == c_str(nn),:);
        end
        
        clear nn c_str
    end
    
    % parse number of color and interpolate if necessary
    if ~ischar(c) && size(c,1) == 2
        tmp = c;
        c      = linspace(tmp(1,1), tmp(2,1), C)';
        c(:,2) = linspace(tmp(1,2), tmp(2,2), C)';
        c(:,3) = linspace(tmp(1,3), tmp(2,3), C)';
    elseif ~ischar(c) && numel(c) == 1
        c = [c c c];
    end
    
    % correct number of colors if to many or to little colors were passed
    if any(strcmpi(plot_type, {'itd', 'ild', 'toa'}))
        c = c(1,:);
    else
        if ~ischar(c) && size(c,1)~=C
            if size(c,1) > C
                c = c(1:C,:);
            else
                c = repmat(c, [ceil(C/size(c,1)), 1]);
                c = c(1:C,:);
            end
        end
    end
    
    if ~ischar(c) && any(c(:)>1)
        c = c/255;
    end
end

% make figure if none exists
if isempty(findall(0,'Type','Figure'))
    AKf
end

% ------------------ 2D plots ------------------ %
if strcmpi(plot_look, '2D')
    % plot
    eval(['hold ' overlay])
    
    % logarithmic x axis for frequency plots
    if strcmpi(plot_domain, 'freq')
        set(gca, 'xscale', 'log')
    end
    
    % limits, grid and tick
    if length(n) > 1
        xlim([min(x) max(x)])
    end
    if ~islogical(dr)
        ylim([dr(1) dr(2)])
    end
    if strcmpi(plot_domain, 'freq')
        set(gca, 'xTick', f_tick, 'xtickLabel', f_label)
    end
    grid on
    box on
    
    % plot
    if any(dash)
        for m = 1:size(data,2)
            h(m) = dashline(n, data(:,m), dash(1), dash(2), dash(3), dash(4), 'LineWidth', lw, 'LineStyle', ls, 'Marker', lm, 'HandleVisibility', 'off'); %#ok<AGROW>
            h2(m) = plot([0 4],[NaN NaN], 'LineWidth', lw, 'LineStyle', '--', 'Marker', lm);    %#ok<AGROW> % plot to get legend entry...
        end
    else
        h = plot(n, data, 'LineWidth', lw, 'LineStyle', ls, 'Marker', lm);
        h2 = h;
    end
    hold off
    
    % color
    if ~strcmpi(c, 'cyc')
        if size(c,1) ~= C
            set(h, 'color', c)
            set(h2, 'color', c)
        else
            for k = 1:C
                set(h(k), 'color', c(k,:))
                set(h2(k), 'color', c(k,:))
            end
            clear k
        end
    end
    
    % labeling
    if strcmpi(labeling, 'on')
        xlabel(x_label, 'fontsize', f_size)
        ylabel(d_type, 'fontsize', f_size)
        title(t_str, 'fontsize', f_size)
    end
    
    % remove tick marks (the grid already marks everything)
    set(gca, 'TickLength', [0 0])
    
    % set output struct
    dataOut.data = data;
    dataOut.n    = n;
    dataOut.h    = h;
    
    clear h2
 
    
% ------------------ 3D plots ------------------ %
elseif strcmpi(plot_look, '3D')
    
    % default data range
    if islogical(y)
        y = 1:C;
    end
    
    % clip data for plotting (do not clip the data that is returned by the function)
    data_p = data;
    if ~islogical(dr)
        data_p = max(dr(1), data_p);
        data_p = min(dr(2), data_p);
    end
    
    % generate mesh
    if strcmpi(plot_func, 'surf')
        [Y, N] = meshgrid(y, n);
        
        % surf with some hardcoded EdgeColor
        eval(['hold ' overlay])
        h = surf(N, Y, data_p, 'EdgeColor', 'none');
        hold off
        
        % set view
        if ischar(hp_view)
            if strcmpi(hp_view, 'top_h')
                if y(1)<y(2)
                    view([0 -90])
                else
                    view([0 90])
                end
            elseif strcmpi(hp_view, 'top_v')
                if y(1)<y(2)
                    view([90 -90])
                else
                    view([-90 90])
                end
            elseif strcmpi(hp_view, 'side')
                view([0 0])
            else
                error('AKp:InputViewType',[hp_view 'is not a valid view type. See doc hp for valid view types'])
            end
        else
            view([hp_view(1) hp_view(2)])
        end
        
        % logarithmic x axis for frequency plots
        if strcmpi(plot_domain, 'freq')
            set(gca, 'xscale', 'log')
        end
        
        % limits, grid and tick
        grid off
        set(gca, 'Layer', 'top');
        box on
        if strcmpi(plot_domain, 'freq')
            set(gca, 'xTick', f_tick, 'xtickLabel', f_label)
        end
        xlim([x(1) x(2)])
        ylim([min(y) max(y)])
        % clipping visualization of z axis (using 'dr') and defining colormap
        if islogical(dr)
            % if not specified, find min and max within given x-range (time, or
            % freqeuncy axis) and use that for clipping
            [~,a] = min(abs(x(1)-n));
            [~,b] = min(abs(x(2)-n));
            dr = [min(min(data(a:b,:))) max(max(data(a:b,:)))];
            % avoid zmin = inf
            dr(isinf(dr)) = min_dB;
            clear a b
        end
        
        % colorbar and colormap
        [dr,c_obj] = AKpSetColormapAndColorbar(cm, cb, cr, dr);
        set(gca, 'zlim', [dr(1) dr(2)]);
        
        % label colorbar
        if any(cb) && strcmpi(labeling, 'on')
            if ~isempty(ct) && ischar(ct)
                ylabel(c_obj, ct, 'fontsize', f_size);
            elseif ~ct
            else
                ylabel(c_obj, d_type, 'fontsize', f_size);
            end
        end
        
        if cf
            if cb ~= 0
                try %#ok<*TRYNC>
                    cbfreeze
                end
            end
            freezeColors
        end
        
        % labeling
        if strcmpi(labeling, 'on')
            xlabel(x_label, 'fontsize', f_size)
            ylabel('data points', 'fontsize', f_size)
            title(t_str, 'fontsize', f_size)
        end
        
        % set output struct
        dataOut.data = data;
        dataOut.N    = N;
        dataOut.Y    = Y;
        dataOut.h    = h;
        
    elseif strcmpi(plot_func, 'img')
        
        % first sample to plot
        plot_start = 1;
        
        % logarithmic x axis for frequency plots
        if strcmpi(plot_domain, 'freq')
            if n(1)==0
                n_img = logspace(log10(n(2)),log10(n(end)),numel(n)-1)';
                plot_start = 2;
            else
                n_img = logspace(log10(n(1)),log10(n(end)),numel(n))';
            end
        else
            n_img = n;
        end
        
        % interpolate data to log grid
        if strcmpi(plot_domain, 'freq')
            data_p = interp1(n ,data_p, n_img, 'spline');
        end
        
        % imagesc
        eval(['hold ' overlay])
        if strcmpi(hp_view, 'top_h')
            if(y(1)<y(end))
                if strcmpi(plot_domain, 'freq')
                    h = imagesc(plot_start:numel(n), -y, data_p');
                else
                    h = imagesc(n(plot_start:end), -y, data_p');
                end
                % yLimits
                ylim([min(-y) max(-y)])
                % set corect yTick
                y_t = -get(gca, 'yTick');
                set(gca, 'yTickLabel', y_t)
            else
                if strcmpi(plot_domain, 'freq')
                    h = imagesc(plot_start:numel(n), y, data_p');
                else
                    h = imagesc(n(plot_start:end), y, data_p');
                end
                % yLimits
                ylim([min(y) max(y)])
            end
        elseif strcmpi(hp_view, 'top_v')
            if(y(1)<y(end))
                if strcmpi(plot_domain, 'freq')
                    h = imagesc(y, plot_start:numel(n), data_p);
                else
                    h = imagesc(y, n(plot_start:end), data_p);
                end
                % yLimits
                xlim([min(y) max(y)])
            else
                if strcmpi(plot_domain, 'freq')
                    h = imagesc(-y, plot_start:numel(n), data_p);
                else
                    h = imagesc(-y, n(plot_start:end), data_p);
                end
                % yLimits
                xlim([min(-y) max(-y)])
                % set corect yTick
                y_t = -get(gca, 'xTick');
                set(gca, 'xTickLabel', y_t)
            end
        end
        hold off
        
        grid off
        set(gca, 'Layer', 'top');
        box on
        % tick for frequency axis
        if strcmpi(plot_domain, 'freq')
            x_t = f_tick(f_tick >= min(x) & f_tick <= max(x));
            x_l = f_label(f_tick >= min(x) & f_tick <= max(x));
            for m = 1:numel(x_t)
                [~, x_t(m)] = min(abs(x_t(m)-n_img));
            end
            
            if strcmpi(hp_view, 'top_h')
                set(gca, 'xTick', x_t, 'xtickLabel', x_l)
            elseif strcmpi(hp_view, 'top_v')
                set(gca, 'yTick', x_t, 'ytickLabel', x_l)
            end
            
            clear x_t x_l
        end
        
        % xlimits
        if strcmpi(hp_view, 'top_h')
            if strcmpi(plot_domain, 'freq')
                xlim([find(n_img>x(1),1,'first') find(n_img<x(2),1,'last')])
            else
                xlim([x(1) x(2)])
            end
        elseif strcmpi(hp_view, 'top_v')
            if strcmpi(plot_domain, 'freq')
                ylim([find(n_img>x(1),1,'first') find(n_img<x(2),1,'last')])
            else
                ylim([x(1) x(2)])
            end
        end
        
        % clipping visualization of z axis (using 'dr') and defining colormap
        if islogical(dr)
            % if not specified, find min and max within given x-range (time, or
            % freqeuncy axis) and use that for clipping
            [~,a] = min(abs(x(1)-n));
            [~,b] = min(abs(x(2)-n));
            dr = [min(min(data(a:b,:))) max(max(data(a:b,:)))];
            % avoid zmin = inf
            dr(isinf(dr)) = min_dB;
            clear a b
        end
        
        % colorbar and colormap
        [~,c_obj] = AKpSetColormapAndColorbar(cm, cb, cr, dr);
        
        % label colorbar
        if any(cb) && strcmpi(labeling, 'on')
            if ~isempty(ct) && ischar(ct)
                ylabel(c_obj, ct, 'fontsize', f_size);
            elseif ~ct
            else
                ylabel(c_obj, d_type, 'fontsize', f_size);
            end
        end
        
        if cf
            if cb ~= 0
                try %#ok<*TRYNC>
                    cbfreeze
                end
            end
            freezeColors
        end
        
        % labeling
        if strcmpi(labeling, 'on')
            if strcmpi(hp_view, 'top_h')
                xlabel(x_label, 'fontsize', f_size)
                ylabel('data points', 'fontsize', f_size)
            elseif strcmpi(hp_view, 'top_v')
                ylabel(x_label, 'fontsize', f_size)
                xlabel('data points', 'fontsize', f_size)
            end
            title(t_str, 'fontsize', f_size)
        end
        
        % set output struct
        dataOut.data = data;
        dataOut.N    = n;
        dataOut.Y    = y;
        dataOut.h    = h;
        
    else
        error('hp:input', ['hp: ' plot_func ' is not a valid argument for ''plot_func'''])
    end

    
% ------------------ spherical plots ------------------ %
elseif ~ischar(plot_look)
    
    % ------------------ balloon plots ------------------ %
    if plot_look <= 4 || plot_look == 7 
        
        % discard oddly shaped triangles
        if triPop > 0 && strcmpi(sph_proc, 'tri')
            % triangle focus points
            triFocus        = mean(X(ch),2);
            triFocus(:,:,2) = mean(Y(ch),2);
            triFocus(:,:,3) = mean(Z(ch),2);
            
            % max distance from vertices of each triangle to focus point 
            triDeformation = max([sum(AKm(triFocus(:,:,1), X(ch), '-').^2,2)  ...
                                  sum(AKm(triFocus(:,:,2), Y(ch), '-').^2,2)  ...
                                  sum(AKm(triFocus(:,:,3), Z(ch), '-').^2,2)  ...
                                 ], [], 2);
            
            %normalize by mean and reference value (chosen by try and error)
            triDeformation = triDeformation ./ mean(triDeformation) / 4;
            
            % discard oddly shaped triangles
            ch = ch(triDeformation <= 1/triPop,:);
            clear triFocus triDeformation
        end
        
        % surf with some hardcoded EdgeColor
        eval(['hold ' overlay])
        
        % plot with or without color information
        if plot_look == 4
            if strcmpi(sph_proc, 'tri')
                if ~ischar(c)
                    h = trisurf(ch, X, Y, Z, data_c, 'EdgeColor', 'none', 'FaceColor', c);
                else
                    h = trisurf(ch, X, Y, Z, data_c, 'EdgeColor', 'none', 'FaceColor', [1 121/255 4/255]);
                end
            else
                if ~ischar(c)
                    if strcmpi(plot_type,'x') % Isaac 4/10/2020: added this to prevent error with 'x4' plots
                        c = c(1,:);
                    end
                    h = surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', c);
                else
                    h = surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', [1 121/255 4/255]);
                end
            end
        elseif plot_look == 7
            if strcmpi(sph_proc, 'tri')
                h = trisurf(ch, X, Y, Z, sign(data_c), 'EdgeColor', 'none');
            else
                h = surf(X, Y, Z, sign(data_c), 'EdgeColor', 'none');
            end
        else
            if strcmpi(sph_proc, 'tri')
                h = trisurf(ch, X, Y, Z, data_c, 'EdgeColor', 'none');
            else
                h = surf(X, Y, Z, data_c, 'EdgeColor', 'none');
            end
        end
        
        axis equal
        axis off
        box off
        hold off
        rotate3d on
                
        % colorbar and colormap or fixed color
        if plot_look == 4
            light;
            lighting phong;
        elseif plot_look == 7
            light;
            lighting phong;
            [~, c_obj] = AKpSetColormapAndColorbar(cm, cb, cr, [-1.3 1.3]);
            
            % labeling color bar
            if any(cb) && strcmpi(labeling, 'on')
                if ~isempty(ct) && ischar(ct)
                    ylabel(c_obj, ct, 'fontsize', f_size);
                elseif ~ct
                else
                    ylabel(c_obj, d_type, 'fontsize', f_size);
                end
            end
            
        else
            if islogical(dr)
                clear dr
                dr(1) = min(min(data_c));
                dr(2) = max(max(data_c));
            end
            if dr(1) == dr(2)
                dr(1) = dr(1)-1;
                dr(2) = dr(2)+1;
            end
            [~, c_obj] = AKpSetColormapAndColorbar(cm, cb, cr, dr);
            
            % labeling colorbar
            if any(cb) && strcmpi(labeling, 'on')
                if ~isempty(ct) && ischar(ct)
                    ylabel(c_obj, ct, 'fontsize', f_size);
                elseif ~ct
                else
                    ylabel(c_obj, d_type, 'fontsize', f_size);
                end
            end
            
        end
        
        if cf
            if cb ~= 0
                try %#ok<*TRYNC>
                    cbfreeze
                end
            end
            freezeColors
        end
        
        % draw coordinate system and limit x,y,z axis
        co.L = 1.2 * max(max(sqrt(X.^2 + Y.^2 + Z.^2)));
        % in case your data is zero dB everywhere, the axis limits are set
        % to -1 and 1.
        if co.L == 0
            co.L = 1;
        end
        if axis_s
            AKpDrawCoordinates(co, axis_s, axis_m, axis_l)
        end
        axis([-co.L co.L -co.L co.L -co.L co.L])
        
        % set view and use default if not specified
        if ischar(hp_view)
            if strcmpi(hp_view, 'top_h') || strcmpi(hp_view, 'top_v')
                view([0 90])
            elseif strcmpi(hp_view, 'side')
                view([0 0])
            else
                error('AKp:InputViewType',[hp_view 'is not a valid view type. See doc hp for valid view types'])
            end
        elseif ~islogical(hp_view)
            view([hp_view(1) hp_view(2)])
        end
        
        % plot title (axis are not labeled because they are off)
        if ischar(sph_f)
            title_str = [sph_f ' '];
        elseif sph_f < 10000 && ~any(strcmpi(plot_type, {'toa' 'itd' 'ild'}))
            title_str = ['f = ' num2str(round(sph_f)) ' Hz '];
        elseif sph_f >= 10000 && ~any(strcmpi(plot_type, {'toa' 'itd' 'ild'}))
            title_str = ['f = ' num2str(round(sph_f/10)*10/1000) ' kHz '];
        else
            title_str = '';
        end
        
        if strcmpi(plot_type, 'm')
            tmp = 'magnitude';
        elseif strcmpi(plot_type, 'p' )
            tmp = 'phase';
        elseif strcmpi(plot_type, 'pu')
            tmp = 'phase';
        else
            tmp = upper(plot_type);
        end
        
        if strcmpi(labeling, 'on')
            if ~strcmpi(plot_domain, 'x')
                switch plot_look
                    case 1; title([title_str '(radius is fixed)']);
                    case 2; title([title_str '(radius is ' tmp ')']);
                    case 3; title([title_str '(radius is magnitude)']);
                    case 4; title([title_str '(radius is ' tmp ')']);
                    case 7; title([title_str '(radius is ' tmp ', color : red=+ blue=-)']);
                end
            else
                switch plot_look
                    case 1; title(['color is ' tmp]);
                    case 2; title(['color is ' tmp]);
                    case 3; title('color is phase');
                    case 4; title('color is fixed');
                    case 7; title('color : red=+ blue=-');
                end
            end
            clear title_str tmp
        end
        
        [az_tmp, el_tmp] = meshgrid((0:360)*pi/180, (90:-1:-90)*pi/180);
        
        % set output struct
        if strcmpi(sph_proc, 'interp')
            dataOut.data = data_c;
            dataOut.az   = az_tmp;
            dataOut.el   = el_tmp;
            dataOut.h    = h;
        else
            dataOut.data = data_c;
            dataOut.az   = co.az;
            dataOut.el   = co.el;
            dataOut.h    = h;
        end
        
        
    % ------------------ planar plots ------------------ %
    elseif plot_look == 5
        
        % plot data
        if strcmpi(sph_proc, 'tri')
            % get azimuth and elevation of triangles defined by tri
            az_tmp = az(tri(:,:))';
            el_tmp = el(tri(:,:))';
            
            % get corresponding color data (the mean of all three points
            % that constitute a triangle is used)
            c = zeros(3, size(az_tmp,2));
            for n = 1:size(az_tmp,2)
                for m = 1:3
                    c(m, n) = data(az == az_tmp(m,n) & el == el_tmp(m,n));
                end
            end
            
            % plot using patch
            h = patch(az_tmp, el_tmp, mean(c), 'EdgeColor', 'none');
            
            view([0 90])
            xlim([min(az) max(az)])
            ylim([min(el) max(el)])
        else
            h = imagesc(data);
        end
        
        % colormap and colorbar
        if islogical(dr)
            clear dr
            dr(1) = min(min(data));
            dr(2) = max(max(data));
        end
        if dr(1) == dr(2)
            dr(1) = dr(1)-1;
            dr(2) = dr(2)+1;
        end
        [~, c_obj] = AKpSetColormapAndColorbar(cm, cb, cr, dr);
        
        % labeling colorbar
        if any(cb) && strcmpi(labeling, 'on')
            if ~isempty(ct) && ischar(ct)
                ylabel(c_obj, ct, 'fontsize', f_size);
            elseif ~ct
            else
                ylabel(c_obj, d_type, 'fontsize', f_size);
            end
        end
        
        if cf
            if cb ~= 0
                try %#ok<*TRYNC>
                    cbfreeze
                end
            end
            freezeColors
        end
        
        % x and y axis labeling and ticks
        if strcmpi(sph_proc, 'interp')
            set(gca, 'xTick', co.xtick, 'xTickLabel', co.xticklabel,...
                     'yTick', co.ytick, 'yTickLabel', co.yticklabel)
        end
        
        if strcmpi(labeling, 'on')
            xlabel 'Azimuth'
            ylabel 'Elevation'
        end
        box on
        
        % set output struct
        dataOut.data = data;
        dataOut.az   = co.az;
        dataOut.el   = co.el;
        dataOut.h    = h;
        
        
        
    % ------------------ polar plots ------------------ %
    elseif plot_look == 6
        
        eval(['hold ' overlay])
        
        % try to plot using matlabs polarplot, and use mmpolar if it does
        % not exist
        if exist('polarplot', 'file')
            
            % turn off axis if it is not a polar axis
            if ~strcmpi(class(gca), 'matlab.graphics.axis.PolarAxes')
                axis off
                pax = polaraxes;
            else
                pax = gca;
            end
            
            %plot and set axis properties
            h = polarplot(co.az, data);
            
            pax.ThetaZeroLocation = 'top';
            
            ThetaTicks = pax.ThetaTick;
            ThetaTicksCell = cell(size(ThetaTicks));
            for nn = 1:numel(ThetaTicks)
                ThetaTicksCell{nn} = [num2str(ThetaTicks(nn)) '\circ'];
            end
            pax.ThetaTickLabel = ThetaTicksCell;
            
            if ~islogical(dr)
                rlim(dr)
            end
            
        elseif exist('mmpolar.m', 'file')
            
            if islogical(dr)
                h = mmpolar(co.az, data, 'Style', 'compass');
            else
                h = mmpolar(co.az, data, 'Style', 'compass', 'Rlimit', dr);
            end
            
        else
            error('AKp:dependencies', 'This plot looks needs Matlabs ''polarplot'' or ''mmpolar'' from the Mathworks file exchange: https://www.mathworks.com/matlabcentral/fileexchange/38855-comprehensive-polar-plots/content/mmpolar.m')
        end
        hold off
        
        % set line color and width
        if strcmpi(c, 'cyc')
            c = 'k';
        end
        set(h, 'color', c(1,:), 'lineWidth', lw, 'LineStyle', ls, 'marker', lm)
    
        % set output struct
        dataOut.data      = data;
        dataOut.az        = co.az;
        dataOut.h         = h;
        if exist('pax', 'var')
            dataOut.polarAxes = pax;
        end
        
    % ------------------ wrong type ------------------ %
    else
        error('AKp:InputViewType',[plot_look 'is not a valid spherical plot. See doc hp for valid plots'])
    end
end

set(gcf, 'color', [1 1 1])
warning('on', 'MATLAB:warn_r14_stucture_assignment');


%% ---------- 10. delete output argument, and plot handles if not requested
if ~islogical(hv)
    if isfield(dataOut, 'h')
        set(dataOut.h, 'handleVisibility', 'off')
    end
end

if nargout == 0
    clear dataOut
end

end





%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% downloaded from http://de.mathworks.com/matlabcentral/fileexchange/1892-dashline
% no licence provided :(
function [h1, h2]=dashline(xdata, ydata, dash1, gap1, dash2, gap2, varargin)
%DASHLINE  Function to produce accurate dotted and dashed lines
%   DASHLINE(XDATA, YDATA, DASH1, GAP1, DASH2, GAP2, ...) plots 
%   the data XDATA and YDATA usig a two dash linestyle. The 
%   function allows the lengths of the dashes to be specified. 
%   DASH1 is the length of the first dash; GAP1 is the length of 
%   the first gap; DASH2 is the length of the second  dash and 
%   GAP2 is the length of the second gap, with all lengths being 
%   in millimeters. The two dash pattern is repeated along the 
%   length of the line. The routine is designed to allow greater
%   control over the linestyle than is available using PLOT.
%
%   Additional inputs may be given to specify line properties. 
%   These are passed through to the plot function. 
%
%   DASH1 or DASH2 may be strings specifying plot symbols. A 
%   dashlength of 0 will cause a dot to be plotted. If GAP1 and 
%   GAP2 are 0 then the arguments are simply passed through to a 
%   plot command
%
%   For example: 
%   clf;
%   dashline([1:10],rand(1,10),1,1,1,1) 
%      %produces a dashed line with 1mm dashes and 1mm gaps
%   hold on;
%   dashline([1:10],rand(1,10),4,2,'o',2,'markerfacecolor','y',...
%      'color','k') 
%      %produces a line with black dashes and yellow centred circles
%
%   [H1, H2]=DASHLINE(...) Outputs are handles to the plotted 
%      lines and markers
%
%   DASHLINE works by calculating the distance along the line to be 
%   plotted, using the existing axes position. If HOLD is OFF when 
%   DASHLINE is called, then the axis limits are set automatically; 
%   otherwise the existing limits are used. DASHLINE then 
%   interpolates XDATA and YDATA at the positions of the start and 
%   end of the dashes, using NaN's to lift the pen. Rescaling the 
%   axes after DASHLINE has been called, or changing the limits, or 
%   changing from log to linear plots, will mean that the dashes no 
%   longer have the specified lengths. This routine does not calculate 
%   distances correctly if the axes DataAspectRatio or 
%   PlotBoxAspectRatio have been manually specified. The lines will 
%   not display properly in a call to LEGEND.
%
%   Kluged together by Edward Abraham (e.abraham@niwa.cri.nz) 
%   27-Sept-2002 Minor change made so that it works with Matlab 6.5
%   27-June-2002 Original version




%Check shape of input data, and put into column form ..
[si, sj]=size(xdata);
[siy, sjy]=size(ydata);
if (si*sj ~= siy*sjy || ((si~=1 && sj~=1) && (si~=siy || sj~=sjy)))
    error('AKp:dashline','Vectors must be the same lengths.')
end
if ( si==1  || siy==1 )
    xdata=xdata(:);
    ydata=ydata(:);
    si=sj;
    sj=1;
end

%Check inputs for the length of the dashes and gaps are sensible ...
if ischar(gap1) || ischar(gap2) || ~isreal(gap1) || ~isreal(gap2) || ~isfinite(gap1) || ~isfinite(gap2) || gap1<0 || gap2<0
    error('AKp:dashline','Gaps must be positive lengths')    
end
if (~ischar(dash1) && ( ~isreal(dash1) || ~isfinite(dash1) || dash1<0)) || (~ischar(dash2) && ( ~isreal(dash2) || ~isfinite(dash2) || dash1<0))
    error('AKp:dashline','Dashes must either be positive lengths or plot strings.')
end   
% If there are no gaps between the dashes then pass data straight through to a plot comand
if (gap1==0 && gap2==0)
    p=plot(xdata,ydata,'-',varargin{:});
    if (nargout>0)
        h1=p;
        h2=[];
    end
    return
end

%Check inputs for string specification of the dashes (dashs of zero length are plotted as dots) ...
if ischar(dash1)
    Marker1=dash1;
    dash1=0;
elseif dash1==0
    Marker1='.';
end
if ischar(dash2)
    Marker2=dash2;
    dash2=0;
elseif dash2==0
    Marker2='.';
end


%Get Axes properties ...
AxesUnits=get(gca,'units');
IsXLog=strcmp(get(gca,'XScale'),'log');
IsYLog=strcmp(get(gca,'YScale'),'log');
IsHold=ishold;
set(gca,'units','centimeters');
Position=get(gca,'Position');

try
    % If hold is off then determine axes limits by initially plotting the data ..
    if ~IsHold
        cla
        p=plot(xdata, ydata);
        if IsXLog
            set(gca,'Xscale','log');
        end
        if IsYLog
            set(gca,'Yscale','log');
        end
    end
    hold on;
    XLim=get(gca,'Xlim');
    YLim=get(gca,'Ylim');
    if ~IsHold
        delete(p)
    end

    % Try to correct for the annoying fact that in Log mode the axes limits are not always the axes limits!
    if IsXLog && XLim(1)==0
        XTick=get(gca,'xtick');
        XLim(1)=XTick(1);
    end
    if IsYLog && YLim(1)==0
        YTick=get(gca,'ytick');
        YLim(1)=YTick(1);
    end


    %Work out position of datapoints...
    if ~IsXLog
        xpos=(xdata-XLim(1))/(XLim(2)-XLim(1))*Position(3);
    else
        xpos=(log10(xdata)-log10(XLim(1)))/(log10(XLim(2))-log10(XLim(1)))*Position(3);
    end

    if ~IsYLog
        ypos=(ydata-YLim(1))/(YLim(2)-YLim(1))*Position(4);
    else
        ypos=(log10(ydata)-log10(YLim(1)))/(log10(YLim(2))-log10(YLim(1)))*Position(4);
    end

    handles=NaN*ones(sj,3);
    %Process each column of data ...
    for i=1:sj
        xposi=xpos(:,i);
        yposi=ypos(:,i);
        xdatai=xdata(:,i);
        ydatai=ydata(:,i);
        f=find(~isreal(xposi) | isinf(xposi)  | isnan(xposi) | ~isreal(yposi) | isinf(yposi)  | isnan(yposi));
        if ~isempty(f)
            xposi(f)=[];
            yposi(f)=[];
            xdatai(f)=[];
            ydatai(f)=[];
        end
        %Calculate distance from the start of the line (in mm) ...
        dist=[0;cumsum(sqrt(diff(xposi).^2 + diff(yposi).^2))*10];

        start1=0:dash1+gap1+dash2+gap2:dist(end);
        dashes=zeros(6*length(start1),1);
        dashes(1:6:end)=start1;
        dashes(2:6:end)=start1+dash1;
        dashes(3:6:end)=NaN;
        dashes(4:6:end)=start1+dash1+gap1;
        dashes(5:6:end)=start1+dash1+gap1+dash2;
        dashes(6:6:end)=NaN;
    
        xdash=NaN*zeros(length(dashes)+length(xdata),1);
        ydash=NaN*zeros(length(dashes)+length(xdata),1);
        %Straight dashes ...
        if ~IsXLog
            xdash(1:length(dashes))=interp1(dist, xdatai, dashes);
        else
            xdash(1:length(dashes))=10.^interp1(dist, log10(xdatai), dashes);
        end 
        if ~IsYLog
            ydash(1:length(dashes))=interp1(dist, ydatai, dashes);
        else
            ydash(1:length(dashes))=10.^interp1(dist, log10(ydatai), dashes);
        end 
        % Get data for markers ...
        if (dash1==0)
            xdot1=xdash(1:6:end);
            ydot1=ydash(1:6:end);
            dashes(1:6:end)=NaN;
            dashes(2:6:end)=NaN;
            xdash(1:6:end)=NaN;
            xdash(2:6:end)=NaN;
            ydash(1:6:end)=NaN;
            ydash(2:6:end)=NaN;
        end
        if (dash2==0)
            xdot2=xdash(4:6:end);
            ydot2=ydash(4:6:end);
            dashes(4:6:end)=NaN;
            dashes(5:6:end)=NaN;
            xdash(4:6:end)=NaN;
            xdash(5:6:end)=NaN;
            ydash(4:6:end)=NaN;
            ydash(5:6:end)=NaN;
        end
    
        %Insert data points that fall within dashes (allows dashes to curve...)
        count=0;
        xlen=length(xdash);
        dashstart=zeros(2*length(start1),1);
        dashstart(1:2:end)=1:6:length(dashes);
        dashstart(2:2:end)=4:6:length(dashes);
        for j=1:length(dashstart)
            f=find(dist>dashes(dashstart(j)) & dist<dashes(dashstart(j)+1));
            if ~isempty(f)
                xdash(dashstart(j)+count+length(f)+1:xlen+length(f))=xdash(dashstart(j)+count+1:xlen);
                xdash(dashstart(j)+count+1:dashstart(j)+count+length(f))=xdatai(f);
                ydash(dashstart(j)+count+length(f)+1:xlen+length(f))=ydash(dashstart(j)+count+1:xlen);
                ydash(dashstart(j)+count+1:dashstart(j)+count+length(f))=ydatai(f);
                xlen=xlen+length(f);
                count=count+length(f);
            end
        end
    
    
        %Plot line and markers ...
        if (~(dash1==0 && dash2==0))
            handles(i,1)=plot(xdash, ydash, varargin{:});
        end
        if (dash1==0)
            handles(i,2)=plot(xdot1, ydot1,Marker1,varargin{:});
        end
        if (dash2==0)
            handles(i,3)=plot(xdot2, ydot2,Marker2,varargin{:});
        end
    end

    % Only return two handles (h2 is empty if there was only one call to plot for each column of data) ...
    if (nargout>0)
        handles(isnan(handles))=[];
        h1=handles(:,1);
        [~, hj]=size(handles);
        if (hj>1)
            h2=handles(:,2);
        else
            h2=[];
        end
    end

    %Restore Axes properties ...
    set(gca,'Units',AxesUnits,'XLim',XLim,'YLim',YLim);
    if ~IsHold
        hold off
    end
catch ME
    %Restore Axes properties ...
    set(gca,'Units',AxesUnits,'XLim',XLim,'YLim',YLim);
    if ~IsHold
        hold off
    end
    rethrow(ME);
end
end
