function [bm,env,GFB] = may2011_gammatone(in,fs,fLow,fUp,nFilter,bEar,bAlign,bInfo)
%gammatone   Auditory filterbank.
% 
%USAGE
%   [BM,ENV,GFB] = may2011_gammatone(IN)
%   [BM,ENV,GFB] = may2011_gammatone(IN,GFB)
%   [BM,ENV,GFB] = may2011_gammatone(IN,FS,FLOW,FUP,NFILTER,bEAR,bALIGN,bINFO)
% 
%INPUT PARAMETER
%        IN : audio object. The gammatone parameter will be initialized
%             using the audio sampling frequency.
% 
%        IN : mono input signal [nSamples x nChannels]
%       GFB : gammatone parameter structure (see gammatoneInit) 
% 
%        IN : multi-channel input signal [nSamples x nChannels]
%        FS : sampling frequency in Hz
%      FLOW : center frequency of lowest auditory filter 
%             (default, FLOW = 80)
%       FUP : center frequency of highest auditory filter 
%             (default, FUP = 5e3)
%   NFILTER : number of auditory filters which will be spaced linear in 
%             the ERB domain. If NFILTER is not a scalar but a vector, the  
%             first value is assumed to represent the number of auditory 
%             channels and the following values represents the indices of 
%             the filters which should be processed. This can be useful if
%             a large number of channels is required as MATLAB might run  
%             out of memory for long signals if all filters should be 
%             computed in one step. 
%             (default, NFILTER = round(freq2erb(FS/2))
%      bEAR : Adjust gain coefficients of the auditory channels to 
%             incorporate middle ear effects. Note that this feature can 
%             only be used if the center frequencies of the auditory 
%             channels are above 20 Hz and below 12500 Hz.
%             (default, bEAR = true)
%    bALIGN : phase-aligned gammatone output (non-causal output)
%             (default, bALIGN = false)
%     bINFO : info flag printing gammatone parameters on the screen
%             (default, bINFO = false)
% 
%OUTPUT PARAMETER
%        BM : basilar membrane displacement [nSamples x nFilter]
%       ENV : instantaneous envelope        [nSamples x nFilter]
%       GFB : gammatone parameter structure 
% 
%NOTE 
%   If only the basilar membrane output "bm" is required, computing
%   BM = gammatone(...);
%   will be significantly faster then computing the envelope as well:
%   [BM,ENV] = gammatone(...);
% 
%REFERENCES
%   The MEX implementaion is based on the source code of the ohio-state 
%   university: www.cse.ohio-state.edu/pnl/shareware/roman-jasa03/
% 
%EXAMPLE
%   nSamples = 500;
%   % Initialize gammatone parameter structure
%   GFB = gammatoneInit(20e3);
%   % Filter impulse with gammatone filtering 
%   bm = gammatone([1; zeros(nSamples-1,1)],GFB);
%   % Plot result
%   waveplot(1:nSamples,GFB.cf,bm);
% 
%   See also gammatoneInit and isGFB.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/may2011_gammatone.php

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

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, 2008-2009
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.0.1   2008/01/21
%   v.0.2   2008/02/04 added phase alignment 
%   v.0.3   2008/02/09 fixed bug related to center frequency calculation
%   v.0.4   2008/05/02 embedded phase alignment into MEX file
%   v.0.5   2009/01/24 added audio object and multi-channel support
%   v.0.6   2009/02/03 check if envelope is required
%   v.0.7   2009/10/10 cleaned up
%   v.0.8   2009/10/12 added consistency check for channel selection
%   ***********************************************************************

% Check for proper input arguments
if nargin < 2 || nargin > 8
    help(mfilename);
    error('Wrong number of input arguments!');
end

% If two parameters are supplied, check if the second argument is a
% gammatone filterbank structure...
if nargin == 2 && isGFB(fs);
    % Set initialization flag to true
    bInitOK = true; 
else
    % Initialization not completed yet...
    bInitOK = false; 
end
   
% Check initialization flag
if bInitOK
    % Copy gammatone filterbank structure
    GFB = fs;
else
    % Set default values
    if nargin < 3 || isempty(fLow);    fLow    = 80;                    end
    if nargin < 4 || isempty(fUp);     fUp     = min(5e3,fs/2);         end
    if nargin < 5 || isempty(nFilter); nFilter = round(freq2erb(fs/2)); end
    if nargin < 6 || isempty(bEar);    bEar    = false;                 end
    if nargin < 7 || isempty(bAlign);  bAlign  = false;                 end
    if nargin < 8 || isempty(bInfo);   bInfo   = false;                 end
    
    % Initialize gammatone filterbank structure
    GFB = may2011_gammatoneinit(fs,fLow,fUp,nFilter,bEar,bAlign,bInfo);
end

% Check for consistent gammatone channel selection
if length(GFB.filter2Process) > 1 && ...
    GFB.filter2Process(1)     <   max(GFB.filter2Process(2:end))
    error(['Inconsistent channel selection. nFilter(1) represents ' ,...
           'the total number auditory filters which are linearly '  ,...
           'spaced on the ERB axis, whereas the remaining indices ' ,...
           'nFilter(2:end) correspond to the auditory channels '    ,...
           'which should be processed.'])
end

% Determine data size
[nSamples,nChannels] = size(in);

% TODO % Enable chunk-based processing by keeping track of filter states 

% Check number of output arguments
if nargout > 1
    % =====================================================================
    % Envelope required    
    % =====================================================================
    %
    % Allocate memory
    [bm,env]=deal(zeros(nSamples,length(GFB.filter2Process)-1,nChannels));
    
    % Loop over number of audio channels
    for ii = 1 : nChannels
        % replace the original MEX filename by the AMT MEX filename      
        if strcmp(GFB.fcnHandle,'gammatoneMEX'), GFB.fcnHandle='comp_may2011_gammatone'; end
        % Call gammatone filterbank routine (MEX-file implementation)
        [bm(:,:,ii),env(:,:,ii)] = feval(GFB.fcnHandle,in(:,ii),GFB.fs,...
                                         GFB.lowerFreq,GFB.upperFreq,...
                                         GFB.filter2Process,...
                                         GFB.bOuterMiddleEar,...
                                         GFB.bPhaseAlign,GFB.bInfo);
    end
else
    % =====================================================================
    % Compute basilar membrane response only
    % =====================================================================
    %
    % Allocate memory
    bm = zeros(nSamples,length(GFB.filter2Process)-1,nChannels);

     % Loop over number of audio channels
    for ii = 1 : nChannels
        % replace the original MEX filename by the AMT MEX filename
        if strcmp(GFB.fcnHandle,'gammatoneMEX'), GFB.fcnHandle='comp_may2011_gammatone'; end
        % Call gammatone filterbank routine (MEX-file implementation)
        bm(:,:,ii) = feval(GFB.fcnHandle,in(:,ii),GFB.fs,...
                           GFB.lowerFreq,GFB.upperFreq,... 
                           GFB.filter2Process,GFB.bOuterMiddleEar,...
                           GFB.bPhaseAlign,GFB.bInfo);
    end
end

% TODO % Adjust outer/middle ear gain here to allow various weightings



function out = isGFB(in)
%isGFB   Check if input is a gammatone filterbank structure. 
%   This is a small helper function in order to check if gammatone
%   filterbank is initialized properly. 
%
%USAGE
%   OUT = isGFB(IN)
%   
%INPUT ARGUMENTS
%    IN : input 
% 
%OUTPUT ARGUMENTS
%   OUT : true/false depending on whether IN is a gammatone structure
% 
%EXAMPLE
%   % Create gammatone structure
%   p = gammatoneInit(44.1e3);
%   % Check 
%   isGFB(p)
%ans = 
%      1
% 
%   See also gammatone.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, 2008 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.0.1   2008/05/11
%   ***********************************************************************


% Check for proper input arguments
if nargin ~= 1
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Initialize output
out = false; 

% Check if IN is a structure
if isstruct(in)
    % Required structure fields
    reqFields = {'fcnHandle' 'fs' 'lowerFreq' 'upperFreq' ...
                 'filter2Process' 'bOuterMiddleEar' 'bPhaseAlign' 'bInfo'};
    
    % Check if all required fields are present
    if all(isfield(in,reqFields))
        % Set flag to true
        out = true;
    end
end
    
