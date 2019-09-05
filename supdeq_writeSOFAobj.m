%% SUpDEq - Spatial Upsampling by Directional Equalization
% 
% function SOFAobj = supdeq_writeSOFAobj(HRIR_L, HRIR_R , samplingGrid, fs, earDistance, sourceDistance)
%
% This function can be used to generate a SOFA object based on a
% pre-defined spherical grid and the corresponding HRIRs.
%
% Output:
% SOFAobj according to the SOFA convention "SimpleFreeFieldHRIR".
% Specific attributes and metadata can be added manually afterwards.
% Use SOFAsave to save SOFA object with specific filename.
%
% Input:
% HRIR_L/R      - 2D Array [N X M] with
%                 N = Sample point on spherical grid, M = Samples per IR, 
%                 Example: HRIR_L/R with dimensions [2702x256]--> 2702
%                 HRIRs (according to samplingGrid with 2702 points) with
%                 256 samples per HRIR
% samplingGrid  - 2D Array [AZ x EL] with
%                 AZ = Azimuth, EL = Elevation in degree [DEG]!
% fs            - Sampling rate 
%                 Default: 48000
% earDistance   - earDistance of the receiver, only needed for HRIRs/BRIRs
%                 Two receivers (=ears) on a head with radius H (in meter): ReceiverPosition=[0 -H 0; 0 +H 0].
%                 Default: 0.165
% sourceDistance- Distance between receiver and source. Only important for
%                 measurement-based SOFA files (HRIRs/BRIRs/Array-IRs).
%                 Default: 3.00 m - In the far field...
%
% Dependencies: SOFA API
%
% SOFA conventions - https://www.sofaconventions.org/mediawiki/index.php/SOFA_conventions
%
% Standardized SOFA conventions
% Standardized SOFA conventions are those which have been standardized by the AES. As with AES69-2015, we have:
% GeneralFIR: General convention with FIR as DataType (no restrictions but DataType)
% GeneralTF: General convention with TF as DataType (no restrictions but DataType)
% SimpleFreeFieldHRIR: Free-field HRTFs stored as impulse responses, measured with an omnidirectional source for a single listener.
% Note: Any modification in one of these conventions changes its status to "stable", i.e., not standardized anymore, unless the modification will be approved the AES.
%
% Stable SOFA conventions
% Stable SOFA conventions are those for which SOFA files are publicly available and can be read/modified by at least one publicly available software package.
% GeneralFIRE: General convention with FIRE as DataType (no restrictions but DataType).
% SimpleHeadphoneIR: Conventions to store headphone IRs recorded for each emitter and each ear, single listener and no directionality of emitter/receiver considered.
% SimpleFreeFieldTF: as SimpleFreeFieldHRIR, but uses TF as DataType covering special needs coming from HRTF simulations.
% SimpleFreeFieldSOS: as SimpleFreeFieldHRIR, but uses SOS as DataType (second-order sections) covering special needs coming from HRTF rendering.
% SingleRoomDRIR: directional room impulse responses (DRIRs) measured with an arbitrary number of receivers (such as a microphone array) and an omnidirectional source in a single room.
% MultiSpeakerBRIR: binaural room impulse responses (BRIRs) measured with an arbitrary number of emitters (such as a loudspeaker array).
% These conventions, when used often and widely, can be proposed for standardization to the AES.
%
% Proposed SOFA conventions
% SimpleBRIR: Binaural room impulse responses measured with an omnidirectional source in a single reverberant space. Somebody wanted to have this, but the work stopped at the moment.
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function SOFAobj = supdeq_writeSOFAobj(HRIR_L, HRIR_R, samplingGrid, fs, earDistance, sourceDistance)

if nargin < 4 || isempty(fs)
    fs = 48000;
end

if nargin < 5 || isempty(earDistance)
    earDistance = 0.165;
end

if nargin < 6 || isempty(sourceDistance)
    sourceDistance = 3.00;
end

%Write IRs 3D Array
IRs = cat(3,HRIR_L',HRIR_R');

if size(IRs,2) ~= size(samplingGrid,1)
    error('Number of IRs must match number to number of sampling grid points!');
end


%%
%Write SOFA file
disp('Writing "SimpleFreeFieldHRIR" SOFA object');

%Get empty SOFA object
SOFAobj = SOFAgetConventions('SimpleFreeFieldHRIR');

%Get number of grid points
nPoints = size(samplingGrid,1);

%Transform grid to SOFA format
samplingGrid(:,2) = 90-samplingGrid(:,2);

%Fill obj with data
SOFAobj.Data.IR = IRs(:,:,1); %Channel1 - Left 
SOFAobj.Data.IR(:,:,2) = IRs(:,:,2); %Channel2 - Right
SOFAobj.Data.IR = shiftdim(SOFAobj.Data.IR,1); % convert from [N M R] to [M R N]
SOFAobj.Data.SamplingRate = fs;
%Fill the mandatory variables (0 -radius 0; 0 radius 0)
SOFAobj.ReceiverPosition = [0 -earDistance/2 0; 0 earDistance/2 0];  
SOFAobj.ListenerPosition = [0 0 0];
SOFAobj.ListenerView = [1 0 0];
SOFAobj.ListenerUp = [0 0 1];
SOFAobj.SourcePosition = [...
            samplingGrid(:,1) ... % AZ in DEG
            samplingGrid(:,2) ... % EL in DEG
            sourceDistance*ones(nPoints,1)];   % nPoints source distance in m
%Update dimensions
SOFAobj=SOFAupdateDimensions(SOFAobj);

disp('Done...');
end

