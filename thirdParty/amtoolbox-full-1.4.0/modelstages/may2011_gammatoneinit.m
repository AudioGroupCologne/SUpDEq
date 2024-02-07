function GFB = may2011_gammatoneinit(fs,lowFreq,upFreq,nFilter,bUseEar,bAlign,bInfo)
%MAY2011_GAMMATONEINIT   Initialize gammatone filterbank structure
%
%   Usage:
%     GFB = may2011_gammatoneinit(FS,FLOW,FUP,NFILTER,BEAR,BALIGN,BINFO)
% 
%   Input parameters:
%     FS      : sampling frequency in Hz
%     FLOW    : center frequency of lowest auditory filter 
%               (default, FLOW = 80)
%     FUP     : center frequency of highest auditory filter 
%               (default, FUP = 8e3)
%     NFILTER : number of auditory filters which will be spaced linear on the 
%               ERB scale. If NFILTER is not a scalar but a vector, the first 
%               value is assumed to represent the number of auditory channels
%               and the following values represents the indices of the 
%               filters which should be processed. This can be useful if a 
%               large number of channels is required as MATLAB might run out 
%               of memory for a longer exerpt if all filters are computed in 
%               one step. (default, NFILTER = round(freq2erb(FS/2))
%     BEAR :    Adjust gain coefficients of the auditory channels to 
%               incorporate middle ear effects. Note that this feature can 
%               only be used if the center frequencies of the auditory 
%               channels are above 20 hertz and below 12500 hertz.
%               (default, BEAR = true)
%     BALIGN :  time-aligned gammatone output (non-causal output)
%               (default, BALIGN = false)
%     BINFO :   info flag printing gammatone parameters on the screen
%               (default, BINFO = false)
% 
%   Output parameters:
%     GFB : gammatone parameter structure which can be passed as second
%           input argument to the function gammatone
%
%   MAY2011_GAMMATONEINIT initializes a gammatone structure suitable
%   for usage with may2011_gammatone.
%
%   Examples:
%
%     nSamples = 500;
%     % Initialize gammatone parameter structure
%     GFB = gammatoneInit(20e3);
%     % Filter impulse with gammatone filtering 
%     bm = gammatone([1; zeros(nSamples-1,1)],GFB);
%     % Plot result
%     waveplot(1:nSamples,GFB.cf,bm);
%
%   See also: gammatone
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/may2011_gammatoneinit.php


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


% Check for proper input arguments
if nargin < 1 || nargin > 7
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Set default values ...   
if nargin < 2 || isempty(lowFreq);  lowFreq = 80;                    end
if nargin < 3 || isempty(upFreq);   upFreq  = min(5e3,fs/2);         end
if nargin < 4 || isempty(nFilter);  nFilter = round(freq2erb(fs/2)); end
if nargin < 5 || isempty(bUseEar);  bUseEar = false;                 end
if nargin < 6 || isempty(bAlign);   bAlign  = false;                 end
if nargin < 7 || isempty(bInfo);    bInfo   = false;                 end

% First value corresponds to the total number of auditory filters
nTotalFilter = nFilter(1);

% Figure out how many auditory filters should be processed
if length(nFilter) == 1
    filter2Process = 1:nTotalFilter;
else
    filter2Process = nFilter(2:end);
end

% Transform frequencies to erb domain
lowerERB = freq2erb(lowFreq);
upperERB = freq2erb(upFreq);

% Calculate center frequencies being linear spaced on the ERB scale
cf = erb2freq(linspace(lowerERB,upperERB,nTotalFilter));

% Bandwidth correction factor (related to 4rth order gamma function)
bandwidthCorrection = 1.019;

% Bandwidth
bw = erb(cf) * bandwidthCorrection;

% Compute gammatone envelope delay in samples (for 4th order)
delay = (3 * fs)./(bw * 2 * pi);

% Phase compensation factor
pc = -cf*3./bw;

% Store gammatone-related parameter
GFB = struct('object','4th order gammatone filterbank',             ...
             'fcnHandle','comp_may2011_gammatone','fs',fs,                    ...
             'nFilter',nTotalFilter,                                ...
             'filter2Process',[nTotalFilter filter2Process],        ...
             'lowerFreq',lowFreq,'upperFreq',upFreq,'cf',cf,        ...
             'bw',bw,'delay',delay,'phaseCorrection',pc,            ...
             'bOuterMiddleEar',bUseEar,'bPhaseAlign',bAlign,'bInfo',bInfo);

        
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
    


function freqHz = erb2freq(erbf)
%erb2freq   Convert ERB rate to frequency in Hz.
%
%   Transformation is done according to Moore and Glasberg.
% 
%USAGE
%   FREQHZ = erb2freq(ERBF);
%
%INPUT ARGUMENTS
%     ERBF : ERB-warped frequencies which should be transformned to
%            frequency in Hz
%    
%OUTPUT ARGUMENTS
%   FREQHZ : frequency vector in Hz
%
%   See also erb and freq2erb.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, 2008 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.1.0   2008/04/02
%   ***********************************************************************

% Check for proper input arguments
if nargin ~= 1
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Convert erb to frequency
freqHz = (10.^(erbf/21.4)-1)/4.37e-3;


function erbf = freq2erb(freqHz)
%freq2erb   Convert frequency in Hz to ERB rate.
% 
%   Transformation is done according to Moore and Glasberg.
%
%USAGE
%     ERBF = freq2erb(FREQHZ);
%
%INPUT ARGUMENTS
%   FREQHZ : frequencies in Hz for which the erb width should be computed
%    
%OUTPUT ARGUMENTS
%     ERBF : auditory filter width 
%
%   See also erb and erb2freq.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, 2008 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.1.0   2008/05/08
%   ***********************************************************************

% Check for proper input arguments
if nargin ~= 1
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Convert frequency to erb
erbf = 21.4*log10(freqHz*0.00437 + 1.0);


function erbWidth = erb(freqHz)
%erb   Compute Equivalent Rectangular Bandwidth (ERB).
% 
%   The ERB is calculated according to Moore and Glasberg at each given
%   frequency.  
%
%USAGE
%   ERBWIDTH = ERB(FREQHZ);
%
%INPUT ARGUMENTS
%     FREQHZ : frequencies in hertz for which the erb width should be
%              computed
%    
%OUTPUT ARGUMENTS
%   ERBWIDTH : auditory filter width 
%
%   See also freq2erb and erb2freq.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, 2008 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.1.0   2008/05/08
%   ***********************************************************************

% Check for proper input arguments
if nargin ~= 1
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Compute ERB width
erbWidth = 24.7*(4.37e-3*freqHz+1);


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************


