%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [denseHRTFdataset, denseHRIRdataset, denseHRTFdataset_sfd] = supdeq_deq(eqHRTFdataset, deqDataset, N, samplingGrid, windowLength, normalization, parallax)
%
% This function performs the de-equalization of a equalized sparse HRTF 
% dataset, here called eqHRTFdataset, (in SH-domain / stored as SH-coefficients) 
% with the pre-defined deqDataset (dataset for de-equalization in SH-domain / 
% stored as SH-coefficients). The result is a dense hrtf dataset (according
% to the (dense) spatial sampling grid), saved as in structs as HRTFs
% (denseHRTFdataset), HRIRs (denseHRIRdataset) or as SH-coefficients
% (denseHRTFdataset_sfd)
%
% Output:
% denseHRTFdataset      - Dense HRTF dataset as a struct 
% denseHRIRdataset      - Dense HRIR dataset as a struct
% denseHRTFdataset_sfd  - Dense HRIR dataset as SH-coefficients,
%                         transformed with given order N
%
% Input:        
% eqHRTFdataset         - Struct with (SH-Coefficients) of the equalized 
%                         HRTF dataset for the left (Hl_nm) and right (Hr_nm) 
%                         channel/ear, absolute frequency scale f, 
%                         transform order N, and FFToversize
% deqDataset            - Struct with de-equalization dataset in SH-domain / as 
%                         SH-coefficients. Can be the output of supdeq_getEqDataset.
% N                     - Transform order N
% samplingGrid          - Spatial sampling grid (Q x 2 or Q x 3 matrix), 
%                         where the first column holds the azimuth, 
%                         the second the elevation (both in degree), and
%                         optionally the third the sampling weights.
%                         Azimuth in degree
%                         (0=front, 90=left, 180=back, 270=right)
%                         (0 points to positive x-axis, 90 to positive y-axis)
%                         Elevations in degree
%                         (0=North Pole, 90=front, 180=South Pole)
%                         (0 points to positive z-axis, 180 to negative z-axis)
% windowLength          - Values of head and tail window in samples [a,b].
%                         Head and tail window will be a half sided hann
%                         window with the defined length. Windowing will be
%                         applied to the HRIRs. The HRTF dataset in frequency- 
%                         and SH-domain will be calculated based on these 
%                         windowed HRIRs. If the matrix is empty ([]), no
%                         windowing will be applied!
%                         Default: [], no windowing
% normalization         - Linear value (for example 0.99) which can be
%                         applied in order to normalize the de-equalized
%                         HRIRs. Can be useful if de-equalization with
%                         waveType == 1 (spherical wave) is applied.
%                         Default: [], no normalization.
% parallax              - Boolean with true or false. Only applies if 
%                         waveType in deq is 1 (spherical wave) and a
%                         distance shift is applied (DVF).
%                         false: acoustic parallax not considered
%                         true (Default: acoustic parallax considered
%
% Dependencies: SOFiA toolbox
%
% (C) 2018 by JMA, Johannes M. Arend
%             CP,  Christoph Pörschmann
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function [denseHRTFdataset, denseHRIRdataset, denseHRTFdataset_sfd] = supdeq_deq(eqHRTFdataset, deqDataset, N, samplingGrid, windowLength, normalization, parallax)

if nargin < 5 || isempty(windowLength)
    windowLength = [];
end

if nargin <6 || isempty(normalization)
    normalization = [];
end

if nargin < 7 || isempty(parallax)
    parallax = true;
end

%Check if eqDataset was with limiting but deqDataset is without limiting
if isfield(eqHRTFdataset,'limited')
    if eqHRTFdataset.limited
        if isfield(deqDataset,'limited')
            if ~deqDataset.limited
                error('Limited eq transfer functions applied for equalization, but not applied for de-equalization!');
            end
        else
            error('Limited eq transfer functions applied for equalization, but not applied for de-equalization!');
        end
    end
end

%Check if applied radius and fA is the same
if isfield(eqHRTFdataset,'limited') && isfield(deqDataset,'limited')
    if eqHRTFdataset.limited && deqDataset.limited
        if eqHRTFdataset.limitInfo.appliedRadius ~= deqDataset.limitInfo.appliedRadius
            error('Applied radius of limited equalization and de-equalization dataset must be the same!');
        end
        if eqHRTFdataset.limitInfo.fA ~= deqDataset.limitInfo.fA
            error('fA of limited equalization and de-equalization dataset must be the same!');
        end
    end
end

%Get fs
fs = eqHRTFdataset.f(end)*2;

%% PROCESSING
%Processing for de-equalization with plane wave deqDataset
if deqDataset.waveType == 0

    %Transform deqDataset to Fourier domain at dense sampling grid points
    %(inverse spherical Fourier transform)
    fprintf('Extracting %d deq transfer functions. This may take some time...\n',size(samplingGrid,1));
    [deqTF_L,deqTF_R] = supdeq_getArbHRTF(deqDataset,samplingGrid,[],[],'ak'); %Use AKtools...
    fprintf('%d deq transfer functions extracted...\n',size(samplingGrid,1))     
    
    %Get HRTF from equalized HRTF dataset by inverse spherical Fourier transform at
    %sampling points of dense sampling grid
    fprintf('Extracting %d HRTFs. This may take some time...\n',size(samplingGrid,1));
    [eqHRTF_L,eqHRTF_R] = supdeq_getArbHRTF(eqHRTFdataset,samplingGrid,[],[],'ak'); %Use AKtools...
    fprintf('%d HRTFs extracted...\n',size(samplingGrid,1))

    %Perform de-equalization (multiplication in frequency domain)
    HRTF_deequalized_L = zeros(size(samplingGrid,1),length(deqDataset.f));
    HRTF_deequalized_R = zeros(size(samplingGrid,1),length(deqDataset.f));
    for kk = 1:size(samplingGrid,1)
        HRTF_deequalized_L(kk,:)=eqHRTF_L(kk,:).*deqTF_L(kk,:);
        HRTF_deequalized_R(kk,:)=eqHRTF_R(kk,:).*deqTF_R(kk,:); 

        if ~mod(kk,100)
            fprintf('%d of %d HRTFs de-equalized...\n',kk,size(samplingGrid,1))
        end

    end

end


%Processing for de-equalization with spherical wave deqDataset (DVF - Distance Variation Function)
if deqDataset.waveType == 1
    
    fprintf('Applying DVF...\n');

    %Transform deqDataset to Fourier domain at dense sampling grid points
    %(inverse spherical Fourier transform)
    fprintf('Extracting %d deq transfer functions. This may take some time...\n',size(samplingGrid,1));
    [deqTF_L,deqTF_R] = supdeq_getArbHRTF(deqDataset,samplingGrid,[],[],'ak'); %Use AKtools...
    fprintf('%d deq transfer functions extracted...\n',size(samplingGrid,1))
    
    if parallax %Consider acoustic parallax
        
        fprintf('Considering acoustic parallax...\n');
    
        %Factor acoustic parallax effect into dense sampling grid
        %(1) - Transform sampling grid and ear positions to radiant in matlab coordinates system
        samplingGridRad(:,1:2) = samplingGrid(:,1:2)*pi/180;
        samplingGridRad(:,2) = pi/2 - samplingGridRad(:,2);
        earPosition_L = deqDataset.earPosition(1:2)*pi/180;
        earPosition_L(2) = pi/2 - earPosition_L(2);
        earPosition_R = deqDataset.earPosition(3:4)*pi/180;
        earPosition_R(2) = pi/2 - earPosition_R(2);
        %(2) - Transform sampling grid and ear positions to cartesian coordinates and apply parallax shift
        [x,y,z] = sph2cart(samplingGridRad(:,1),samplingGridRad(:,2),deqDataset.sourceDistance);
        [x_ear_L,y_ear_L,z_ear_L] = sph2cart(earPosition_L(:,1),earPosition_L(:,2),deqDataset.radius); 
        [x_ear_R,y_ear_R,z_ear_R] = sph2cart(earPosition_R(:,1),earPosition_R(:,2),deqDataset.radius);  
        xL = x-x_ear_L;
        xR = x-x_ear_R;
        yL = y-y_ear_L;
        yR = y-y_ear_R;
        zL = z-z_ear_L;
        zR = z-z_ear_R;
        [AzParL,ElParL] = cart2sph(xL,yL,zL);
        [AzParR,ElParR] = cart2sph(xR,yR,zR);
        %(3) - Transform back to degree in SH coordinates system
        AzParL = mod(AzParL,2*pi);
        AzParR = mod(AzParR,2*pi);
        AzParL = AzParL*180/pi;
        AzParR = AzParR*180/pi;
        ElParL=pi/2-ElParL;
        ElParR=pi/2-ElParR;
        ElParL = ElParL*180/pi;
        ElParR = ElParR*180/pi;
        %(4) - Store in new samplingGrids for left and right channel
        samplingGridParL = [AzParL,ElParL];
        samplingGridParR = [AzParR,ElParR];

        %Get HRTF from equalized HRTF dataset by inverse spherical Fourier transform at
        %sampling points of dense sampling grid 
        fprintf('Extracting %d HRTFs. This may take some time...\n',size(samplingGridParL,1));
        %Get left and right HRTFs separately
        [eqHRTF_L,~] = supdeq_getArbHRTF(eqHRTFdataset,samplingGridParL,'DEG',0,'ak'); %Use AKtools...
        [~,eqHRTF_R] = supdeq_getArbHRTF(eqHRTFdataset,samplingGridParR,'DEG',1,'ak'); %Use AKtools...
        fprintf('%d HRTFs extracted...\n',size(samplingGrid,1))
    
    else %parallax = false, do not consider acoustic parallax
        
        fprintf('Not considering acoustic parallax...\n');
        
        %Get HRTF from equalized HRTF dataset by inverse spherical Fourier transform at
        %sampling points of dense sampling grid
        fprintf('Extracting %d HRTFs. This may take some time...\n',size(samplingGrid,1));
        [eqHRTF_L,eqHRTF_R] = supdeq_getArbHRTF(eqHRTFdataset,samplingGrid,[],[],'ak'); %Use AKtools...
        fprintf('%d HRTFs extracted...\n',size(samplingGrid,1))
        
    end

    %Perform de-equalization (multiplication in frequency domain)
    HRTF_deequalized_L = zeros(size(samplingGrid,1),length(deqDataset.f));
    HRTF_deequalized_R = zeros(size(samplingGrid,1),length(deqDataset.f));
    for kk = 1:size(samplingGrid,1)
        HRTF_deequalized_L(kk,:)=eqHRTF_L(kk,:).*deqTF_L(kk,:);
        HRTF_deequalized_R(kk,:)=eqHRTF_R(kk,:).*deqTF_R(kk,:); 

        if ~mod(kk,100)
            fprintf('%d of %d HRTFs de-equalized...\n',kk,size(samplingGrid,1))
        end

    end

end

%% Store in structs

% Transform HRTFs to time domain --> get HRIRs
%Get both sided spectrum
HRTF_deequalized_L_double = conj([HRTF_deequalized_L(:,:), conj(fliplr(HRTF_deequalized_L(:,2:end-1)))])';
HRTF_deequalized_R_double = conj([HRTF_deequalized_R(:,:), conj(fliplr(HRTF_deequalized_R(:,2:end-1)))])';
%Transform back to time domain
HRIR_deequalized_L = real(ifft(HRTF_deequalized_L_double));
HRIR_deequalized_R = real(ifft(HRTF_deequalized_R_double));
%Cut if FFToversize of eqHRTFdataset > 1
if eqHRTFdataset.FFToversize > 1
    newLength = size(HRIR_deequalized_L,1)/eqHRTFdataset.FFToversize;
    HRIR_deequalized_L = HRIR_deequalized_L(1:newLength,:);
    HRIR_deequalized_R = HRIR_deequalized_R(1:newLength,:);
end

%Window if windowLength is not empty
if ~isempty(windowLength) %windowing
    disp('Applying windows...');
    
    %Apply normalization if desired
    if ~isempty(normalization)
        disp('Applying normalization...');
        
        %Get abs max of left/right channel
        maxL = max(max(abs(HRIR_deequalized_L)));
        maxR = max(max(abs(HRIR_deequalized_R)));
        
        %Calculate normalization factor
        if maxL >= maxR
            normFactor = normalization/maxL;
        elseif maxR > maxL
            normFactor = normalization/maxR;
        end
        
        %Apply to HRIRs
        HRIR_deequalized_L = HRIR_deequalized_L*normFactor;
        HRIR_deequalized_R = HRIR_deequalized_R*normFactor; 
    end
    
    %Windowing
    %Get window function in time domain
    win = supdeq_win(size(HRIR_deequalized_L,1),windowLength);
    %Apply to HRIRs
    for kk = 1:size(HRIR_deequalized_L,2)
        HRIR_deequalized_L(:,kk) = HRIR_deequalized_L(:,kk) .* win;
        HRIR_deequalized_R(:,kk) = HRIR_deequalized_R(:,kk) .* win;
    end
    
    %Store in struct denseHRIRdataset (window, normalization optionally)
    denseHRIRdataset.HRIR_L = HRIR_deequalized_L'; %Flip again to have same array settings as in SOFiA
    denseHRIRdataset.HRIR_R = HRIR_deequalized_R';
    denseHRIRdataset.fs = eqHRTFdataset.f(end)*2;
    denseHRIRdataset.samplingGrid = samplingGrid;
    denseHRIRdataset.deqEarDistance = deqDataset.earDistance;
    denseHRIRdataset.deqWaveType  = deqDataset.waveType;
    if deqDataset.waveType == 1
        denseHRIRdataset.sourceDistance = deqDataset.sourceDistance;
    end
    if ~isempty(normalization)
        denseHRIRdataset.normFactor = normFactor;
    end
    denseHRIRdataset.windowLength = windowLength;
    
    %Transform de-equalized windowed denseHRIRdataset to Fourier domain to get
    %denseHRTFdataset (with same FFToversize)
    HRTF_deequalized_win_L = fft(HRIR_deequalized_L',size(HRIR_deequalized_L,1)*eqHRTFdataset.FFToversize,2);
    HRTF_deequalized_win_L = HRTF_deequalized_win_L(:,1:end/2+1);
    HRTF_deequalized_win_R = fft(HRIR_deequalized_R',size(HRIR_deequalized_R,1)*eqHRTFdataset.FFToversize,2);
    HRTF_deequalized_win_R = HRTF_deequalized_win_R(:,1:end/2+1);
    
    %Store in struct denseHRTFdataset (window, normalization optionally)
    denseHRTFdataset.HRTF_L = HRTF_deequalized_win_L;
    denseHRTFdataset.HRTF_R = HRTF_deequalized_win_R;
    denseHRTFdataset.f      = eqHRTFdataset.f;
    denseHRTFdataset.fs     = eqHRTFdataset.f(end)*2;
    denseHRTFdataset.FFToversize = eqHRTFdataset.FFToversize;
    denseHRTFdataset.samplingGrid = samplingGrid;
    denseHRTFdataset.deqEarDistance = deqDataset.earDistance;
    denseHRTFdataset.deqWaveType  = deqDataset.waveType;
    if deqDataset.waveType == 1
        denseHRTFdataset.sourceDistance = deqDataset.sourceDistance;
    end
    if ~isempty(normalization)
        denseHRTFdataset.normFactor = normFactor;
    end
    denseHRTFdataset.windowLength = windowLength;
    
    % Transform de-equalized HRTFs to SH-domain --> get
    % denseHRTFdataset_sfd (window, normalization optionally)
    fprintf('Performing spherical Fourier transform with N = %d. This may take some time...\n',N);
    denseHRTFdataset_sfd = supdeq_hrtf2sfd(HRTF_deequalized_win_L,HRTF_deequalized_win_R,N,samplingGrid,fs,'ak'); %Use AKtools...
    denseHRTFdataset_sfd.FFToversize = eqHRTFdataset.FFToversize;
    denseHRTFdataset_sfd.deqEarDistance = deqDataset.earDistance;
    denseHRTFdataset_sfd.deqWaveType  = deqDataset.waveType;
    if deqDataset.waveType == 1
        denseHRTFdataset_sfd.sourceDistance = deqDataset.sourceDistance;
    end
    if ~isempty(normalization)
        denseHRTFdataset_sfd.normFactor = normFactor;
    end
    denseHRTFdataset_sfd.windowLength = windowLength;
    
    fprintf('Done with de-equalization...\n');
      
else %no windowing
 
%Apply normalization if desired
    if ~isempty(normalization)
        disp('Applying normalization...');
        
        %Get abs max of left/right channel
        maxL = max(max(abs(HRIR_deequalized_L)));
        maxR = max(max(abs(HRIR_deequalized_R)));
        
        %Calculate normalization factor
        if maxL >= maxR
            normFactor = normalization/maxL;
        elseif maxR > maxL
            normFactor = normalization/maxR;
        end
        
        %Apply to HRIRs
        HRIR_deequalized_L = HRIR_deequalized_L*normFactor;
        HRIR_deequalized_R = HRIR_deequalized_R*normFactor;   
        
        %Store in struct denseHRIRdataset (no window, normalization)
        denseHRIRdataset.HRIR_L = HRIR_deequalized_L'; %Flip again to have same array settings as in SOFiA
        denseHRIRdataset.HRIR_R = HRIR_deequalized_R';
        denseHRIRdataset.fs = eqHRTFdataset.f(end)*2;
        denseHRIRdataset.samplingGrid = samplingGrid;
        denseHRIRdataset.deqEarDistance = deqDataset.earDistance;
        denseHRIRdataset.deqWaveType  = deqDataset.waveType;
        if deqDataset.waveType == 1
            denseHRIRdataset.sourceDistance = deqDataset.sourceDistance;
        end
        denseHRIRdataset.normFactor = normFactor;
        
        %Transform de-equalized normalized denseHRIRdataset to Fourier domain to get
        %denseHRTFdataset (with same FFToversize)
        HRTF_deequalized_L = fft(HRIR_deequalized_L',size(HRIR_deequalized_L,1)*eqHRTFdataset.FFToversize,2);
        HRTF_deequalized_L = HRTF_deequalized_L(:,1:end/2+1);
        HRTF_deequalized_R = fft(HRIR_deequalized_R',size(HRIR_deequalized_R,1)*eqHRTFdataset.FFToversize,2);
        HRTF_deequalized_R = HRTF_deequalized_R(:,1:end/2+1);
        
        %Store in struct denseHRTFdataset (no window, normalization)
        denseHRTFdataset.HRTF_L = HRTF_deequalized_L;
        denseHRTFdataset.HRTF_R = HRTF_deequalized_R;
        denseHRTFdataset.f      = eqHRTFdataset.f;
        denseHRTFdataset.fs     = eqHRTFdataset.f(end)*2;
        denseHRTFdataset.FFToversize = eqHRTFdataset.FFToversize;
        denseHRTFdataset.samplingGrid = samplingGrid;
        denseHRTFdataset.deqEarDistance = deqDataset.earDistance;
        denseHRTFdataset.deqWaveType  = deqDataset.waveType;
        if deqDataset.waveType == 1
            denseHRTFdataset.sourceDistance = deqDataset.sourceDistance;
        end
        denseHRTFdataset.normFactor = normFactor;
        
        % Transform de-equalized HRTFs to SH-domain --> get
        % denseHRTFdataset_sfd (no window, normalization)
        fprintf('Performing spherical Fourier transform with N = %d. This may take some time...\n',N);
        denseHRTFdataset_sfd = supdeq_hrtf2sfd(HRTF_deequalized_L,HRTF_deequalized_R,N,samplingGrid,fs,'ak'); %Use AKtools...
        denseHRTFdataset_sfd.FFToversize = eqHRTFdataset.FFToversize;
        denseHRTFdataset_sfd.deqEarDistance = deqDataset.earDistance;
        denseHRTFdataset_sfd.deqWaveType  = deqDataset.waveType;
        if deqDataset.waveType == 1
            denseHRTFdataset_sfd.sourceDistance = deqDataset.sourceDistance;
        end
        denseHRTFdataset_sfd.normFactor = normFactor;
        fprintf('Done with de-equalization...\n');
    
        
    else % No normalization
        %Store in struct denseHRIRdataset (no window, no normalization)
        denseHRIRdataset.HRIR_L = HRIR_deequalized_L'; %Flip again to have same array settings as in SOFiA
        denseHRIRdataset.HRIR_R = HRIR_deequalized_R';
        denseHRIRdataset.fs = eqHRTFdataset.f(end)*2;
        denseHRIRdataset.samplingGrid = samplingGrid;
        denseHRIRdataset.deqEarDistance = deqDataset.earDistance;
        denseHRIRdataset.deqWaveType  = deqDataset.waveType;
        if deqDataset.waveType == 1
            denseHRIRdataset.sourceDistance = deqDataset.sourceDistance;
        end

        %Store in struct denseHRTFdataset (no window, no normalization)
        denseHRTFdataset.HRTF_L = HRTF_deequalized_L;
        denseHRTFdataset.HRTF_R = HRTF_deequalized_R;
        denseHRTFdataset.f      = eqHRTFdataset.f;
        denseHRTFdataset.fs     = eqHRTFdataset.f(end)*2;
        denseHRTFdataset.FFToversize = eqHRTFdataset.FFToversize;
        denseHRTFdataset.samplingGrid = samplingGrid;
        denseHRTFdataset.deqEarDistance = deqDataset.earDistance;
        denseHRTFdataset.deqWaveType  = deqDataset.waveType;
        if deqDataset.waveType == 1
            denseHRTFdataset.sourceDistance = deqDataset.sourceDistance;
        end

        % Transform de-equalized HRTFs to SH-domain --> get
        % denseHRTFdataset_sfd (no window, no normalization)
        fprintf('Performing spherical Fourier transform with N = %d. This may take some time...\n',N);
        denseHRTFdataset_sfd = supdeq_hrtf2sfd(HRTF_deequalized_L,HRTF_deequalized_R,N,samplingGrid,fs,'ak'); %Use AKtools...
        denseHRTFdataset_sfd.FFToversize = eqHRTFdataset.FFToversize;
        denseHRTFdataset_sfd.deqEarDistance = deqDataset.earDistance;
        denseHRTFdataset_sfd.deqWaveType  = deqDataset.waveType;
        if deqDataset.waveType == 1
            denseHRTFdataset_sfd.sourceDistance = deqDataset.sourceDistance;
        end
        fprintf('Done with de-equalization...\n');

    end
end

end

