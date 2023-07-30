function out = may2011(input,fs)
%MAY2011   Azimuthal localization of concurrent talkers
%   Usage: out = may2011(input,fs);
%
%   Input parameters:
%       input : binaural input signal [nSamples x 2 channels]
%
%       fs : sampling frequency in Hertz
%
%   Output parameters:
%     out : localization structure 
%
%   MAY2011 is a probabilistic localization model.
%
%   The out structure has the following fields:
%
%     .param      processing parameters
%
%     .azFrames   frame-based azimuth estimate    [Frames x 1]
%
%     .azimuth    time-frequency azimuth estimate [nFilter x nFrames]
%
%     .rangeAZ    azimuth grid                    [nAzDir x 1]
%
%     .prob       3D probability map              [nFilter x nFrames x nAzDir]
%
%     .loglik     3D log-likelihood map           [nFilter x nFrames x nAzDir]
%
%     .itd        interaural time difference      [nFilter x nFrames]
%
%     .ild        interaural level difference     [nFilter x nFrames]
%
%     .ic         interaural coherence             [nFilter x nFrames]
%
%
%   References:
%     T. May, S. van de Par, and A. Kohlrausch. A probabilistic model for
%     robust localization based on a binaural auditory front-end. IEEE Trans
%     Audio Speech Lang Proc, 19:1--13, 2011.
%     
%     N. Roman, D. L. Wang, and G. J. Brown. Speech segregation based on
%     sound localization. The Journal of the Acoustical Society of America,
%     114(4):2236--2252, 2003.
%     
%
%   See also: may2011_interpolateparabolic
%   may2011_classifygmm
%   may2011_xcorrnorm
%   may2011_findlocalpeaks
%   may2011_cbarlabel
%   may2011_gammatone
%   may2011_framedata
%   may2011_neuraltransduction
%   may2011_estazimuthgmm
%   may2011_fireprint
%   may2011_gammatoneinit
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/may2011.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Tobias May (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   ***********************************************************************

% Initialize persistent memory
persistent AMT_may2011PERC AMT_may2011PERmethod


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin ~= 2
    help(mfilename);
    error('Wrong number of input arguments!');
end


%% *******************  INITIALIZE AZIMUTH CLASSIFIER  ********************
% 
% 
% Block parameter 
blockSec  = 20e-3; % Block size in seconds
hopSec    = 10e-3; % Step size in seconds
winType   = 'rectwin';
methodGMM = 'default'; 

% Select GMM preset 
switch lower(methodGMM)
    case 'default'
        % Room 5 -90:5:90
        A.gmmModel = 'modeldata.mat';   
    otherwise
        error(['GMM module ''',lower(methodGMM),'''is not available.'])
end

% Check if localization module can be re-loaded from persistent memory
if ~isempty(AMT_may2011PERC) && ~isempty(AMT_may2011PERmethod) && isequal(methodGMM,AMT_may2011PERmethod)
    % Re-load GMMs
    C = AMT_may2011PERC;
else
    % Load localization module
    x=amt_load('may2011',A.gmmModel);
    C=x.C;
    
    % Store classifier to persistent memory
    AMT_may2011PERC = C; AMT_may2011PERmethod = methodGMM;
end

% Store classifier
A.GMM = C;

% Extract gammatone filterbank settings from the classifier
GFB = A.GMM.featSpace(1).GFB;
                 
% Adapt fs in GFB structure according to input signal
if ~isequal(GFB.fs,fs)
    GFB.fs = fs;
end

% Haircell processing
if isfield(C.featSpace(1),'neuralTransduction')
    A.haircellModel = C.featSpace(1).neuralTransduction;
else
    A.haircellModel = 'roman';
end

% Cross-correlation method
if isfield(C.featSpace(1),'xcorrMethod') && ...
   ~isempty(C.featSpace(1).xcorrMethod)
    A.xcorrMethod = C.featSpace(1).xcorrMethod;
else
    A.xcorrMethod = 'power';
end

% ITD parameter
A.bXcorrNorm    = true;
A.bXcorrDetrend = true;
A.maxDelay      = round(1e-3 * fs);
A.lags          = (-A.maxDelay:1:A.maxDelay).';
A.domain        = 'time';

% Resolution of GMM models
stepSizeGMM = 1;

% For short access
gmmFinal = C.classifier.gmmFinal;

% Build azimuth grid
azIdx  = 1:stepSizeGMM:length(C.azimuth);
azGrid = C.azimuth(1:stepSizeGMM:end);
azStep = diff(azGrid(1:2));

% Number of different azimuth directions 
nAz = length(azGrid);

% Determine length of input signal
nSamples = size(input,1);

% Block processing parameters
blockSamples = 2 * round(fs * blockSec / 2);
hopSamples   = 2 * round(fs * hopSec / 2); 
overSamples  = blockSamples - hopSamples;

% Calculate number of frames
nFrames = fix((nSamples-overSamples)/hopSamples);

% Allocate memory
prob   = zeros(GFB.nFilter,nFrames,nAz);
loglik = zeros(GFB.nFilter,nFrames,nAz);
feat   = zeros(nFrames,2);
itd    = zeros(GFB.nFilter,nFrames);
ild    = zeros(GFB.nFilter,nFrames);
ic     = zeros(GFB.nFilter,nFrames);


%% ************************  GAMMATONE FILTERING  *************************
%
% 
% Perform gammatone filtering
bm = may2011_gammatone(input,GFB);


%% ********************  BINAURAL FEATURE EXTRACTION  *********************
%
% 
% Loop over number of channels
for ii = 1 : GFB.nFilter
    
    % Neural transduction
    left  = may2011_neuraltransduction(bm(:,ii,1), fs, A.haircellModel);
    right = may2011_neuraltransduction(bm(:,ii,2), fs, A.haircellModel);
    
    % Framing
    frameL = may2011_framedata(left, blockSamples, hopSamples, winType);
    frameR = may2011_framedata(right,blockSamples, hopSamples, winType);
    
    % =====================================================================
    % ITD ESTIMATE
    % =====================================================================
    % 
    % Cross-correlation analysis
    switch lower(A.xcorrMethod)
        case 'power'
            % Time-based
            xcorr = may2011_xcorrnorm(frameL,frameR,A.maxDelay,A.bXcorrDetrend,...
                              A.bXcorrNorm);
    end
    
    % Maximum search per frame
    bM = transpose(argmax(xcorr,1)); 
    
    % Refine lag position by applying parabolic interpolation
    [delta,currIC] = may2011_interpolateparabolic(xcorr,bM);
    
    % Calculate ITD with fractional part
    feat(:,1) = A.lags(bM)/fs + delta * (1/fs);
    
    % =====================================================================
    % ILD ESTIMATE
    % =====================================================================
    %     
    % Calculate energy per frame
    bmEnergyL = sum(abs(frameL).^2)/blockSamples;
    bmEnergyR = sum(abs(frameR).^2)/blockSamples;
    
    % Compute ILD
    feat(:,2) = 10 * (log10(bmEnergyR + eps) -  log10(bmEnergyL + eps));
        
    % Compensate for square-root compression
    if isequal(A.haircellModel,'roman')
        % Reverse effect of sqrt() in decibel domain
        feat(:,2) = 2 * feat(:,2);
    end
    
    % =====================================================================
    % CREATE 3-DIMENSIONAL AZIMUTH-DEPENDENT LIKELIHOOD FUNCTION
    % =====================================================================
    %     
    % Classify azimuth
    prob(ii,:,:) = may2011_classifygmm(feat,gmmFinal{ii}(azIdx));
    
    % Normalize
    prob(ii,:,:) = prob(ii,:,:) ./ repmat(sum(prob(ii,:,:),3),[1 1 nAz]);
    
    % Store log-likelihood
    loglik(ii,:,:) = log(prob(ii,:,:));
    
    % Store ITD
    itd(ii,:) = feat(:,1);
    
    % Store ILD
    ild(ii,:) = feat(:,2);
    
    % Store IC
    ic(ii,:) = currIC;
end


%% ***********************  LOCALIZATION ESTIMATE  ************************
% 
% 
% =========================================================================
% TIME/FREQUENCY-BASED AZIMUTH ESTIMATE
% =========================================================================
% 
% 
% Maximum search
[tmp,maxIdx] = max(loglik,[],3); %#ok

% Warp to azimuth indices
azimuth = azGrid(maxIdx);

% =========================================================================
% FRAME-BASED AZIMUTH ESTIMATE
% =========================================================================
% 
% 
% Integrate log-likelihoods across channels
intChanGMM = squeeze(nanmean(loglik,1));

% Maximum search
[tmp,maxIdx] = max(intChanGMM,[],2); %#ok

% Warp to azimuth indices
azFrames = azGrid(maxIdx);

% Apply interpolation to increase azimuth resolution
delta = may2011_interpolateparabolic(exp(intChanGMM).',maxIdx);

% Find overall azimuth estimate
azFrames(:) = azFrames(:) + (azStep*delta(:));


%% ************************  CREATE RESULT STRUCT  ************************
%
% 
% Create output structure
out.param      = A;
out.azFrames   = azFrames;
out.azimuth    = azimuth;
out.rangeAZ    = azGrid;
out.prob       = prob;
out.loglik     = loglik;
out.itd        = itd;
out.ild        = ild;
out.ic         = ic;


function maxidx = argmax(input, dim)

if nargin < 2 || isempty(dim)
    [temp, maxidx] = max(input);
else
    [temp, maxidx] = max(input,[],dim);
end

function m=nanmean(x,dim)
m=mean(x(~isnan(x)),dim);


