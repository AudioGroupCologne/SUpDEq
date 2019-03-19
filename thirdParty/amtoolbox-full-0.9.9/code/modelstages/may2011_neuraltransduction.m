function audNerve = may2011_neuraltransduction(bm,fs,haircellMethod)

%   Developed with Matlab 7.5.0.342 (R2007b). Please send bug reports to:
%   
%   Author  :  Tobias May, 2008-2009
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :
%   v.1.0   2008/08/06
%   v.1.1   2008/11/05 added envelope compression from Binaural
%                      Cross-correlation Toolbox developed by A. Akeroyd
%   v.1.2   2009/02/10 added persistent memory for envelope lowpass
%   v.1.3   2009/02/11 delay compensation for envelope lowpass filter
%   **********************************************************************
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/may2011_neuraltransduction.php

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

% Initialize persistent memory
persistent PERfs PERlowpass


%% ***********************  CHECK INPUT ARGUMENTS  ************************

% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Set default values
if nargin < 3 || isempty(haircellMethod); haircellMethod = 'none'; end


%% ************************  NEURAL TRANSDUCTION  *************************
% 
% Wrapper for various models of neural transduction
% 
switch lower(haircellMethod)
    case 'none'
        % No processing ...
        audNerve = bm;
    case 'halfwave'
        % ----------------------
        % PALOMAEKI 2004
        % ----------------------
        % Half-wave rectification
        audNerve = max(bm,0);
    case 'roman'
        % ----------------------
        % ROMAN 2003
        % ----------------------
        % Half-wave rectification and square-root compression
        audNerve = sqrt(max(bm,0));
    case 'envelope'
        % Halfwave-rectification amd full envelope compression 
        %
        % The envelope compression itself is from Bernsten, van de Par
        % and Trahiotis (1996, especially the Appendix). The lowpass 
        % filtering is from Berstein and Trahiotis (1996, especially EQ 2 
        % on page 3781).

        % Define lowpass filter
        if (isempty(PERfs) || isempty(PERlowpass)) || ...
           (~isempty(PERfs) && ~isequal(PERfs,fs))
            % Define lowpass filter
            cutoff  = 425; %Hz
            order   = 4;
            lpf     = linspace(0, fs/2, 10000);
            f0      = cutoff * (1./ (2.^(1/order)-1).^0.5);
            lpmag   = 1./ (1+(lpf./f0).^2) .^ (order/2);
            lpf     = lpf ./ (fs/2);
            % Filter design
            lowpass = fir2(256, lpf, lpmag, hamming(257));
            
            % Store lowpass to persistent memory
            PERfs      = fs;
            PERlowpass = lowpass;
        else
            % Reload lowpass from persistent memory
            lowpass = PERlowpass;
        end        
        
        % Envelope compression using Weiss/Rose lowpass filter
        compress1 = 0.23;
        compress2 = 2.0;

        % ========================
        % Do the actual processing
        % ========================
       
        % Get envelope
        envelope = abs(hilbert(bm));
        % compress the envelope to a power of compression1, while 
        % maintaining the fine structure.
        compressedenvelope = (envelope.^(compress1 - 1)) .* bm;
        % rectify that compressed envelope
        rectifiedenvelope = max(compressedenvelope,0);
        % raise to power of compress2
        rectifiedenvelope = rectifiedenvelope.^compress2;
        % Zero-padd data to compensate for the delay
        [tmp,maxIdx] = max(abs(lowpass)); %#ok
        rectifiedenvelope = [rectifiedenvelope;zeros(maxIdx-1,size(bm,2))];
        % overlap-add FIR filter using the fft
        audNerve = fftfilt(lowpass, rectifiedenvelope);
        
        % Trim signal to its original length
        audNerve = audNerve(maxIdx:end,:);
    otherwise
        error(['Neural transduction method ''',haircellMethod,...
               ''' is not supported.'])
end


