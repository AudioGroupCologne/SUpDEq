function audNerve = may2011_neuraltransduction(bm,fs,haircellMethod)
%MAY2011_NEURALTRANSDUCTION calculates the auditory nerve response
%
%   Input parameters:
%     bm             : matrix containing signals after auditory filtering
%     fs             : sampling frequency [Hz]
%     haircellMethod : can be 'none', 'halfwave, 'roman', or 'envelope'
%
%   Output parameters:
%     audNerve       : matrix containing the signals after inner hair cell processing
%
%   MAY2011_NEURALTRANSDUCTION simulates the neural transduction in the inner ear
%   by means of lowpass filtering of the signal envelope
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/may2011_neuraltransduction.php


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



