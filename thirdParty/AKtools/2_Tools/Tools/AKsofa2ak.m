% [h, f, srcPos, recPos, unit] = AKsofa2ak(Obj, type)
% returns the data (impulse responses or complex spectra) from a SOFA 
% object in a format used by AKtools along with the most important meta
% data entries.
%
% For information on SOFA see
% [1] https://www.sofaconventions.org
% [2] AES Standards Comittee (2015): "AES69-2015: AES standard for file 
%     exchange - Spatial acoustic data file format."  Audio Engineering
%     Society.
%
%
% I N P U T:
% Obj    - SOFAobject of one of the following conventions
%          'GeneralFIR'
%          'SimpleFreeFieldHRIR'
%          'SimpleHeadphoneIR'
%          'SingleRoomDRIR'
%          'GeneralFIRE'
%          'MultiSpeakerBRIR'
%          'GeneralTF'
%          'SimpleFreeFieldTF'
% type   - desired type of srcPos and recPos (see below).
%          'cartesian' or 'spherical'
%
% O U T P U T:
% h      - impulse responses or complex spectra of size [N M R E], with
%          N: number of samples or frequency bins
%          M: number of measurements
%          R: number of receivers (e.g. microphone positions)
%          E: number of emitter (e.g. loudspeaker, only for GeneralFIRE and
%             MultiSpeakerBRIR)
% f      - the sampling rate in case Obj hold impulse responses or a vector
%          of frequencies in case Obj holds transfer functions
% srcPos - source positions of size [1 3] or [M 3]. The format of srcPos
%          can either be x/y/z or azimuth/elevation/distance
% recPos - receiver positions of size [R 3 1] or [R 3 M] The format of
%          recPos can either be x/y/z or azimuth/elevation/distance
% unit   - cell string that gives the type (spherical or cartesian) and
%          unit of srcPos, and recPos
%
%
% 3/2015 - fabian.brinkmann@tu-berlin.de

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
function [h, f, srcPos, recPos, unit] = AKsofa2ak(Obj, type)

% Check which SOFA convention we have and pass the corresponding data
if any(strcmpi(Obj.GLOBAL_SOFAConventions, {'GeneralFIR' 'SimpleFreeFieldHRIR' 'SimpleHeadphoneIR' 'SingleRoomDRIR'}))
    h = shiftdim(Obj.Data.IR, 2);
    f = Obj.Data.SamplingRate;
elseif any(strcmpi(Obj.GLOBAL_SOFAConventions, {'GeneralTF' 'SimpleFreeFieldTF'}))
    h = shiftdim(Obj.Data.Real, 2) + 1j*shiftdim(Obj.Data.Imag, 2);
    f = Obj.N;
elseif any(strcmpi(Obj.GLOBAL_SOFAConventions, {'GeneralFIRE' 'MultiSpeakerBRIR'}))
    h = shiftdim(Obj.Data.IR, 3);
    f = Obj.Data.SamplingRate;
else
    error(['SofaConvention ' Obj.GLOBAL_SOFAConventions ' is not supported'])
end

% get source and receiver positions
srcPos = Obj.SourcePosition;
recPos = Obj.ReceiverPosition;

% check if unit is as desired
if exist('type', 'var')
    if ~strcmpi(Obj.SourcePosition_Type, type)
        if strcmpi(type, 'spherical')
            srcPos = ak_spherical(srcPos);
        elseif strcmpi(type, 'cartesian')
            srcPos = ak_cartesian(srcPos);
        else
            error('AKsofa2ak:Input', 'type must be ''spherical'' or ''cartesian''')
        end
    end
    
    if ~strcmpi(Obj.ReceiverPosition_Type, type)
        if strcmpi(type, 'spherical')
            recPos = ak_spherical(recPos);
        elseif strcmpi(type, 'cartesian')
            recPos = ak_cartesian(recPos);
        else
            error('AKsofa2ak:Input', 'type must be ''spherical'' or ''cartesian''')
        end
    end
end

unit = {'srcPos' Obj.SourcePosition_Type   Obj.SourcePosition_Units; ...
        'recPos' Obj.ReceiverPosition_Type Obj.ReceiverPosition_Units};

end

function pos = ak_spherical(pos)
[pos(:,1), pos(:,2), pos(:,3)] = cart2sph(pos(:,1), pos(:,2), pos(:,3));

pos(:, 1:2) = pos(:,1:2) / pi * 180;
end

function pos = ak_cartesian(pos)
[pos(:,1), pos(:,2), pos(:,3)] = sph2cart(pos(:,1)/180*pi, pos(:,2)/180*pi, pos(:,3));
end