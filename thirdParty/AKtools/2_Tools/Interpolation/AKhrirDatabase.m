% [data, sg] = AKhrirDatabase(type, sg, DFC, dataFormat)
% get HRIRs or an auralization of moving sources from the HUTUBS HRIR
% database available. For information on the database see [1, 2]
%
% e.g.
% hrir = AKhrirDatabase('pp1_measured_ir', [0 90]);
% returns the measured HRIR to the left of particpant 1
% see AKhrirDemo.m for more examples
%
%
% I N P U T
% type       - is a string that specifies the participant and data that is
%              returned. It has three subsrings that are separated by
%              undersocres, i.e. 'type1_type2_type3'
%
%               type1: 'pp#' to select a participant from the database,
%                            e.g. 'pp1', or 'pp21'
%               type2: 'measured' returns acoustically measured HRIRs
%                      'simulated' returns numerically simulated HRIRs
%               type3: 'ir' loads the HRIRs from the SOFA files. In this
%                           case the HRIR that is closest to the desired
%                           source position is returend (see 'sg' below)
%                      'sh' loads the HRIRs from the spherical harmonics
%                           coefficients. This is a spatially continuous
%                           representation, and HRIRs at the desired source
%                           positions are returned (see 'sg' below)
%                      'auralization' returns an auralization of a source
%                                     moving in the horiztonal and median
%                                     plane.
% sg         - [Q x 2] vector that specifies the source positions that are
%              returned. Each row contains one source positions given by
%              azimuth and elevation in degree. The azimuth goes from 0 to
%              360 deg., where 0 is to the front and 90 to the left. The
%              elevation goes from 90 to -90 deg., where 90 is above, 0 to
%              the front and -90 below.
%              Note that sg is ignored if type3 = 'auralization'
% DFC        - flag to apply a diffuse field compensation to the HRIRs or
%              auralization (true, false=default)
% dataFormat - 'SOFA' will return the HRIRs in a SOFA file [3]
%            - 'mat' will return the HRIRs in a matrix of size [256, Q, 2]
%                    where 256 is the length of each HRIR in samples, Q are
%                    the number of source positions, and 2 are the HRIRs of
%                    the left and right ear.
%              An auralization will always be of size [464366 x 2].
%
% O U T P U T
% data       - HRIRs or auralization according to 'dataFormat' (see above)
% sg         - source positions of the HRIRs. If type3='ir' the actual
%              source positions can differ from the requested positions.
%              In this case the HRIRs closest to the request source
%              positions are returned in 'data' and 'sg'
% DTF        - Diffuse field transfer function obtained by energetic average
%              across source positions and ears. The DTF is 3rd octave
%              smoothed and minimum-phase.
% DTFinverse - Inverted DTF (3rd octave smoothed and minimum-phase). The
%              gain above 10 kHz is limited to 2 dB to avoid high gains due
%              to the high frequency roll-off in measured HRIRs above
%              approx. 20 kHz.
%
%
% [1] tba
% [2] tba
% [3] AES Standards Comittee (2015), AES69-2015: AES standard for file
%     exchange - Spatial acoustic data file format. Audio Engineering
%     Society .
%
% 09/2018 - fabian.brinkmann@tu-berlin.de

% ToDo: check DemoScript with final database
%       insert URL to database here, in AKdpendencies, in AKtoolsStart
%       add DemoFile
%       add references to this file

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
function [data, sg, DTF, DTFinverse] = AKhrirDatabase(type, sg, DFC, dataFormat)

% set defaults
if nargin < 4
    dataFormat = 'mat';
end
if nargin <3
    DFC = false;
end

% check type
if isstruct(type)
    % data was passed directly
    
    if isfield(type, 'Data')
        % SOFA object was passed
        H = type;
        clear type
        
        type{1} = H.GLOBAL_ListenerShortName;
        type{2} = '?';
        type{3} = 'ir';
        
    else
        % else SH coefficients were passed
        sh = type;
        clear type
        
        type{1} = sh.GLOBAL_ListenerShortName;
        if isfield(sh, 'TOA')
            type{2} = 'measured';
        else
            type{2} = 'simulated';
        end
        type{3} = 'sh';
        
    end
    
else
    % data is loaded based on specification in 'type'
    
    % parse input
    type  = strsplit(type, '_');
    
    if ~contains(type{1}, 'pp')
        error('AKhrirDatabase:Input', 'type must start with pp#_')
    end
    
    if ~any( strcmpi( type{2}, {'measured', 'simulated'} ) )
        error('AKhrirDatabase:Input', 'The second entry of type must be ''measured'' or ''simulated''')
    end
    
    AKdependencies('HUTUBS')
    
    % load data
    switch lower(type{3})
        case 'ir'
            
            H = SOFAload( [type{1} '_HRIRs_' type{2} '.sofa'] );
            
        case {'sh', 'auralization'}
            
            sh = load( [type{1} '_SHcoefficients_' type{2} '.mat'] );
            
        otherwise
            error('AKhrirDatabase:Input', 'The third entry of type must be ''ir'', ''sh'', or ''auralization''')
    end
end

% get HRIRs or auralization
switch lower(type{3})
    case 'ir'
        
        % find closest measured HRIR to requested points
        SG = H.SourcePosition;
        id = nan(size(sg,1),1);
        
        for nn = 1:size(sg,1)
            d = 2 * acosd( sind(sg(nn,2)).*sind(SG(:,2)) +  cosd(sg(nn,2)).*cosd(SG(:,2)).*cosd( abs(sg(nn,1)-SG(:,1)) ));
            [~, ID] = min(d);
            sg(nn,:) = SG(ID,1:2);
            id(nn)   = ID;
        end
        
        % get points from SOFA file
        H.SourcePosition = H.SourcePosition(id,:);
        H.Data.IR        = H.Data.IR(id,:,:);
        H                = SOFAupdateDimensions(H);
        
        % convert to desired output format
        switch lower(dataFormat)
            case 'sofa'
                data = H;
            case 'mat'
                data = AKsofa2ak(H);
            otherwise
                error('AKhrirDatabase:Input', 'dataFormat be ''SOFA'' or ''mat''')
        end
        
    case 'sh'
        
        % get HRIRs
        l = AKisht(sh.HRIR.CoeffLeft, true, [sg(:,1) 90-sg(:,2)], 'complex');
        r = AKisht(sh.HRIR.CoeffRight, true, [sg(:,1) 90-sg(:,2)], 'complex');
        
        % re-insert TOA in case of measured data
        if strcmpi(type{2}, 'measured')
            toa = AKisht(sh.TOA.CoeffLeft, false, [sg(:,1) 90-sg(:,2)], 'complex', false, false, 'real');
            l   = AKfractionalDelayCyclic(l, toa);
            
            toa = AKisht(sh.TOA.CoeffRight, false, [sg(:,1) 90-sg(:,2)], 'complex', false, false, 'real');
            r   = AKfractionalDelayCyclic(r, toa);
        end
        
        switch lower(dataFormat)
            case 'sofa'
                % load and fill empty sofa container
                data                   = SOFAload( [type{1} '_HRIRs_' type{2} '.sofa'], 'nodata');
                data.SourcePosition    = [sg repmat( data.SourcePosition(1,3), [size(sg,1) 1] )];
                data.Data.IR           = shiftdim( cat(3, l, r), 1 );
                data.Data.SamplingRate = 44100;
                data.Data.Delay        = [0 0];
                data                   = SOFAupdateDimensions(data);
                
            case 'mat'
                data = cat(3, l, r);
            otherwise
                error('AKhrirDatabase:Input', 'dataFormat be ''SOFA'' or ''mat''')
        end
        
    case 'auralization'
        
        % auralization of horizontal plane
        y1 = aura(sh, [0 360]', [0 0]');
        
        % auralization of median plane
        y2 = aura(sh, [0 0 NaN 180 180 NaN 0 0]', [0 90 NaN 90 -90 NaN -90 0]');
        
        data = [AKfade(y1, [], 441, 441); zeros(22000, 2); AKfade(y2, [], 441, 441)];
        
    otherwise
end


% apply diffuse field compensation
if DFC || nargout >= 3
    
    % load SH coefficients
    sh = load( [type{1} '_SHcoefficients_' type{2} '.mat'] );
    
    % calculate diffuse field filters
    DTF_l = sqrt( AKshEnergy(sh.HRIR.CoeffLeft) )';
    DTF_r = sqrt( AKshEnergy(sh.HRIR.CoeffLeft) )';
    
    DTF = sqrt( ( DTF_l.^2 + DTF_r.^2 ) / 2 );
    DTF = AKfractOctSmooth(DTF, 'welti', 44100, 3);
    
    DTFinverse = AKsoftLimit(1./DTF, 2, 2, [10e3 22050], 44100);
    
    DTF = AKsingle2bothSidedSpectrum(DTF);
    DTF = ifft(DTF, 'symmetric');
    DTF = AKphaseManipulation(DTF, 44100, 'min', 2, false);
    
    
    DTFinverse = AKsingle2bothSidedSpectrum(DTFinverse);
    DTFinverse = ifft(DTFinverse, 'symmetric');
    DTFinverse = AKphaseManipulation(DTFinverse, 44100, 'min', 2, false);
    
    if DFC
        % get data from SOFA file
        if isstruct(data)
            tmp  = data;
            data = AKsofa2ak(data);
        end
        
        % apply filter
        for cc = 1:size(data,3)
            data(:,:,cc) = fftfilt(DTFinverse, data(:,:,cc));
        end
        
        % get data to SOFA file
        if exist('tmp', 'var')
            tmp.Data.IR = shiftdim(data, 1);
            data        = tmp;
        end
    end
    
end

end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% code for auralization taken from AKshAura and added re-insertion of the
% TOA
function y = aura(sh, az, el)

% block size
N = 512;

% input signal
x = AKnoise(5*44100);
x = [x x];


% find segments in azimuth and elevation
Sid = find(isnan(az)==1); numel(az);
S   = numel(Sid) + 1;
if any(~isnan(el(Sid)))
    error('NaN must be at the same positions in az and el.')
end
Sid = [0; Sid; numel(az)+1];

% get number of Blocks for convolution
L = ceil(size(x,1)/N);
if mod(L, S)
    L = L + S-mod(L, S);
end

% zero Pad input signal to integer divisor of block size N
x(end+1:N*L,:) = 0;

% Interpolate trajectory
azI = zeros(L,1);
elI = azI;
for ss = 1:S
    azI((ss-1)*L/S+1:ss*L/S) = interp1(az(Sid(ss)+1:Sid(ss+1)-1), linspace(1, numel(Sid(ss)+1:Sid(ss+1)-1), L/S));
    elI((ss-1)*L/S+1:ss*L/S) = interp1(el(Sid(ss)+1:Sid(ss+1)-1), linspace(1, numel(Sid(ss)+1:Sid(ss+1)-1), L/S));
end

% get IRs
h1 = AKisht(sh.HRIR.CoeffLeft,  true, [azI 90-elI], 'complex');
h2 = AKisht(sh.HRIR.CoeffRight, true, [azI 90-elI], 'complex');

% re-insert TOA
if isfield(sh, 'TOA')
    toa = AKisht(sh.TOA.CoeffLeft,  false, [azI 90-elI], 'complex', true, false, 'real');
    h1  = AKfractionalDelayCyclic(h1, toa);
    
    toa = AKisht(sh.TOA.CoeffRight, false, [azI 90-elI], 'complex', true, false, 'real');
    h2  = AKfractionalDelayCyclic(h2, toa);
end

% length of IRs
Nh = size(h1,1);

% allocate output signal
y = zeros(N*L+Nh-1, size(x,2));

% fft filt
if size(x,2) == 1
    for ll = 1:L
        y((ll-1)*N+1:ll*N+Nh-1)   = y((ll-1)*N+1:ll*N+Nh-1)   + fftfilt(x((ll-1)*N+1:ll*N),    [h1(:,ll); zeros(N-1, 1)]);
    end
else
    for ll = 1:L
        y((ll-1)*N+1:ll*N+Nh-1,:) = y((ll-1)*N+1:ll*N+Nh-1,:) + fftfilt(x((ll-1)*N+1:ll*N,:), [[h1(:,ll) h2(:,ll)]; zeros(N-1, 2)]);
    end
end

end