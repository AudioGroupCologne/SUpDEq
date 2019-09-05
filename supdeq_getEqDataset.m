%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function eqDataset = supdeq_getEqDataset(N, earDistance, NFFT, fs, waveType, sourceDistance, earPosition, sphereOffset, transformCore, swSettings)
%
% This function returns the sound pressure distribution of an incident 
% plane wave for different sound incidence directions on the left ear 
% position and on the right ear position of a sphere which models the 
% human head.
%
% Output:
% eqDataset     - Struct with SH-coefficients describing the sound 
%                 incidence at the letf/right ear position on the sphere 
%                 (Hl_nm/Hr_nm). Length of the SH-coefficients is NFFT/2+1 
%                 (single sided spectrum). Important input variables are
%                 stored in the struct.
%
% Input:        
% N             - Transform order N
%                 Default: 35
% earDistance   - Distance between both ear positions in m (radius/2)
%                 Default: 0.165
% NFFT          - FFT size of the SH-coefficients to be returned
%                 Default: 512
% fs            - Sampling rate
%                 Default: 48000
% waveType      - Boolean for plane wave (0) or spherical wave (1)
%                 Default: 0
% sourceDistance- Distance (in meter) between source and sphere. 
%                 Only considered if waveType = 1 (spherical wave)
%                 Default: 1
% earPosition   - 4 x 1 row vector describing the position of the ears in
%                 spherical coordinates in degree [azL, elL, azR, elR]
%                 Default: [90, 90, 270, 90] (left-right symmetrical)
% sphereOffset  - 3 x 1 row vector describing the offset of the sphere from
%                 the origin in m [xd, yd, zd], with x shift on "depth-axis" 
%                 of the sphere, y shift on "width-axis" of the sphere, 
%                 and z shift on "height-axis" of the sphere. See matlab 
%                 function "cart2sph" for describtion of the coordinate.
%                 Positive values 
%                 system.
%                 Default: [] - No offset
%                 Combinations of earPosition and sphereOffset are stored
%                 in the field "earPosition" in the output struct
% transformCore - String to define method to be used for the rigid spher
%                 transfer function synthesis
%                 'sofia - sofia_wgc from SOFiA toolbox
%                 'ak'   - AKsphericalHead from AKtools 
%                 The results are almost exactly the same
%                 Default: 'sofia'
% swSettings    - Settings only for sourceType = 1 (spherical wave), 
%                 indicated by integer 0,1,2,or 3.
%                 Reference distance of spherical wave is set to 1.00m for
%                 both transform cores!
%                 0 - Compensate time shift according to distance to
%                 reference distance (1.00m), but maintain level changes.
%                 1 - Compensate time shift and compensate level changes
%                 according to 1/r law (to level at r0 = 1m). The 1/r
%                 compensation is not necessarily accurate enough
%                 in the near-field!
%                 2 - Maintain time shift but compensate level changes 
%                 3 - Maintain time shift and maintain level changes (no
%                 compensation)
%                 Default: 1
%
% IMPORTANT NOTE: ALWAYS USE THE SAME TRANSFORM-CORE FOR EQUALIZATION AND
% DE-EQUALIZATION. AS THE LEVEL AND PHASE PROPERTIES OF THE MODELS ARE 
% DIFFERENT, MIXING THEM UP RESULTS IN INCORRECT RESULTS!
%
% Dependencies: SOFiA toolbox
%
% References:
% Benjamin Bernschütz: Microphone Arrays and Sound Field Decomposition 
% for Dynamic Binaural Recording. Ph.D. dissertation, Technical University
% Berlin (2016).
%   
% R. O. Duda and W. L. Martens, ?Range dependence of the response of a 
% spherical head model,? J. Acoust. Soc. Am., vol. 104, no. 5, 
% pp. 3048?3058, 1998.
%
% (C) 2018 by JMA, Johannes M. Arend
%             CP,  Christoph Pörschmann
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function eqDataset = supdeq_getEqDataset(N, earDistance, NFFT, fs, waveType, sourceDistance, earPosition, sphereOffset, transformCore, swSettings)

if nargin < 1 || isempty(N)
    N = 35;
end

if nargin < 2 || isempty(earDistance)
    earDistance = 0.165;
end

if nargin < 3 || isempty(NFFT)
    NFFT = 512;
end

if nargin < 4 || isempty(fs)
    fs = 48000;
end

if nargin < 5 || isempty(waveType)
    waveType = 0;
    sourceDistance = 1;
end

if nargin < 6 || isempty(sourceDistance)
    sourceDistance = 1;
end

if nargin < 7 || isempty(earPosition)
    earPosition = [90, 90, 270, 90];
end

if nargin < 8 || isempty(sphereOffset)
    sphereOffset = [];
end

if nargin < 9 || isempty(transformCore)
    transformCore = 'sofia';
end

if nargin < 10 || isempty(swSettings)
    swSettings = 1;
end
     
%% Get SH-Coefficients with sofia_wgc

%Define required parameters
radius  = earDistance/2;
if sourceDistance < radius
    error('Invalid source distance. Source must be outside the radius (earDistance/2)');
end
ac      = 2; %Array configuration, 2 - Rigid sphere with pressure transducers
c       = 343; %Speed of sound in m/s - Set to 343 here
%Time delay in s - Set to 0 for plane waves or, to compensate for the time shift induced by spherical wave generation,
%to negative time delay according to swSettings
if waveType == 0  
    delay = 0; 
else %waveType == 1
    
    %This way, sofia and ak IRs are at the same position when reference in ak is set to r0 = 1m. Applied to all sofia data as a global delay
    %0.5141 was determined empirically with several tests...
    globalDelay = -1*((1-0.5141)/c);
    globalDelay2 = -1*((0.5141)/c);
    
    if swSettings == 0 %Compensate time shift but maintain level changes
        delay = (-1*(sourceDistance/c)) - globalDelay2;
    elseif swSettings == 1 %Compensate time shift and compensate level changes (Default)
        delay = (-1*(sourceDistance/c)) - globalDelay2;
        levelComp = sourceDistance/1; %Reference at r0 = 1m
    elseif swSettings == 2 %Maintain time shift but compensate level changes
        delay = globalDelay; %Just the global delay to align sofia and ak
        levelComp = sourceDistance/1; %Reference at r0 = 1m
    elseif swSettings == 3 %Maintain time shift and maintain level changes (no compensation)
        delay = globalDelay; %Just the global delay to align sofia and ak
    end    
    
end

%Convert earPosition to radiant
appliedEarPosition = earPosition*pi/180;

%Shift earPosition if sphereOffset applied
if ~isempty(sphereOffset)
    azL = appliedEarPosition(1);
    elL = appliedEarPosition(2);
    azR = appliedEarPosition(3);
    elR = appliedEarPosition(4);
 
    [xL,yL,zL] = sph2cart(azL,pi/2-elL,radius);
    [xR,yR,zR] = sph2cart(azR,pi/2-elR,radius);
    
    xL = xL - sphereOffset(1); %Use minus because sphere offset is inverse to ear offset
    xR = xR - sphereOffset(1);
    yL = yL - sphereOffset(2);
    yR = yR - sphereOffset(2);
    zL = zL - sphereOffset(3);
    zR = zR - sphereOffset(3);
    
    [azL, elL, ~] = cart2sph(xL,yL,zL);
    azL = mod(azL,2*pi);
    elL = pi/2 - elL;
    [azR, elR, ~] = cart2sph(xR,yR,zR);
    azR = mod(azR,2*pi);
    elR = pi/2 - elR;
    
    appliedEarPosition = [azL,elL,azR,elR];
end

%Calculate SH-Coefficients
if strcmp(transformCore,'sofia')
    eqDataset.Hl_nm = sofia_wgc(N, radius, ac, fs, NFFT, appliedEarPosition(1), appliedEarPosition(2), delay,c,waveType,sourceDistance);
    eqDataset.Hr_nm = sofia_wgc(N, radius, ac, fs, NFFT, appliedEarPosition(3), appliedEarPosition(4), delay,c,waveType,sourceDistance);
    
    %Apply level compensation if chosen with respective swSettings
    if waveType == 1 && (swSettings == 1 || swSettings == 2)
        eqDataset.Hl_nm = eqDataset.Hl_nm*levelComp;
        eqDataset.Hr_nm = eqDataset.Hr_nm*levelComp;
    end
    
    
elseif strcmp(transformCore,'ak')
    
    %Transform coordinates first
    appliedEarPositionAK = appliedEarPosition*180/pi;
    appliedEarPositionAK(2) = 90 - appliedEarPositionAK(2);
    appliedEarPositionAK(4) = 90 - appliedEarPositionAK(4);
    
    %Get Lebedev sampling grid according to N
    lebN = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65];
    lebNumPoints = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810];
    
    if N > lebN(end)
        error('N = 65 is the highest spatial order possible with transformcore ak. Use sofia instead');
    end
 
    %Get closest N and add + 1 to be safe
    [~,cId] = min(abs(N-lebN));
    if lebN(cId) < N
        cId = cId + 1;
    end

    %Get sampling grid
    samplingGridAK = supdeq_lebedev(lebNumPoints(cId));
    samplingGridAKSH = samplingGridAK;
    
    %Transform coordinates again
    samplingGridAK(:,2) = 90-samplingGridAK(:,2);
    samplingGridAK = samplingGridAK(:,1:2);
    
    %Define SH order for AKsphericalHead
    shOrderAK = 100; %Maybe set to lower order?
    
    %Define distance - if wavetype = 0 (plane wave), set distance to 100 m
    if waveType == 0
        distanceAK = ones(size(samplingGridAK,1),1) * 100; 
        samplingGridAK(:,3) = distanceAK;
        %%Define referenceDistance
        refDistance = distanceAK(1);
    end   
        
    if waveType == 1 % wavetype = 1 (spherical wave)
        distanceAK = ones(size(samplingGridAK,1),1) * sourceDistance;
        samplingGridAK(:,3) = distanceAK;
        
        %Define parameters according to swSettings
        if swSettings == 0 %Compensate time shift but maintain level changes
            refDistance = distanceAK(1); %Setting reference distance and source distance to the same value leads to time- and level-compensation
            levelComp = 1/sourceDistance; %Reference at r0 = 1m
        elseif swSettings == 1 %Compensate time shift and compensate level changes (Default)
            refDistance = distanceAK(1); %Setting reference distance and source distance to the same value leads to time- and level-compensation
        elseif swSettings == 2 %Maintain time shift but compensate level changes
            refDistance = 1;
            levelComp = sourceDistance/1; %Reference at r0 = 1m
        elseif swSettings == 3 %Maintain time shift and maintain level changes (no compensation)
            refDistance = 1; %Set reference to 1 m and apply to compensation/shift
        end  
       
    end
    
    %Get IRs
    irsAK = AKsphericalHead(samplingGridAK,appliedEarPositionAK,false,radius,refDistance,shOrderAK,NFFT,fs,c);
    irsAK_L = squeeze(irsAK(:,:,1));
    irsAK_R = squeeze(irsAK(:,:,2));
    
    %Apply level compensation if chosen with respective swSettings
    if waveType == 1 && (swSettings == 0 || swSettings == 2)
        irsAK_L = irsAK_L*levelComp;
        irsAK_R = irsAK_R*levelComp;
    end
    
    %Transform IRs to SH domain with given order N
    eqDataset.Hl_nm = AKsht(irsAK_L,true,samplingGridAKSH,N,'complex',fs);
    eqDataset.Hr_nm = AKsht(irsAK_R,true,samplingGridAKSH,N,'complex',fs);

end
eqDataset.f = linspace(0,fs/2,NFFT/2+1);
eqDataset.N = N;
eqDataset.earDistance = earDistance;
eqDataset.radius = radius;
eqDataset.waveType = waveType;
eqDataset.sourceDistance = sourceDistance;
eqDataset.c = c;
if isempty(sphereOffset)
    eqDataset.earPosition = earPosition;
end
if ~isempty(sphereOffset)
    eqDataset.inputEarPosition = earPosition;
    eqDataset.sphereOffset = sphereOffset;
    eqDataset.appliedEarPosition = appliedEarPosition*180/pi;
end
%if waveType == 1
%    eqDataset.tsComp = tsComp;
%end
eqDataset.transformCore = transformCore;

end

