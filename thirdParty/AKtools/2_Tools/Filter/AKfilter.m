% [y, varOut1] = AKfilter(x, type, f, g, fs, varIn1, varIn2, varIn3)
% 
% filters input signal x, e.g.
% y = AKfilter(x, 'hp', 50, 0, 44100, 4, 'butter')
% is a 4th order butter high pass at 50 Hz.
%
% For more examples see AKfilterDemo.m
% A python implementation of most filters can be found here:
% http://nbviewer.jupyter.org/github/sonible/nb/blob/master/iir_filter/index.ipynb
%
%
% I N P U T
% x            - Input signal to be filtered [n samples x M channels]
% type         - String that specifies the filter type. Different input
%                parameters are needed depending on the type:
%
%                High-pass, Low-pass: 
%                'hp', 'lp'   : varIn1=N, varIn2=design
%
%                Band-pass
%                'bp'         : varIn1=N, varIn2=design
%
%                Band-stop
%                'bs'         : varIn1=N, varIn2=design
%
%                Cross-over network
%                'xover'      : varIn1=N, varIn2=design
%
%                Parametric equalizer
%                'peq'        : varIn1=Q, varIn2=peq_design, varIn3=pre_warp
%
%                Lowshelve, Highshelve
%                'ls', 'hs'   : varIn1=N, varIn2=shelve_design
%
%                Fractional Octave filter bank
%                'fracOct'    : varIn1=N, varIn2=design, varIn3=frac
%
%                FFT brickwall Fractional Octave filter bank
%                'FFTfracOct' : varIn1=frac
%
%                FIR filters used in listening tests from Moore et al. [1]
%                'moore1989'  : varIn1=bw
%
% f            - Caracteristic frequency [Hz]:
%                Cut-off frequency if type is 'hp', or 'lp', or 'xover'
%
%                Lower and upper cut-off freq. if type is 'bp', 'bs'
%
%                Lower and upper limit of filterbank center frequencies if
%                type is 'fracOct', or 'FFTfracOct'
%
%                Mid-gain freq. if type is 'peq'
%
%                High-, mid-, or low-gain freq. if type is 'ls',
%                or 'hs' (depending on shelve_type)
%
% g            - Filter gain [dB]
%
% fs           - Sampling rate [Hz]
%
% N            - Filter order (must be 1 or 2 for 'ls', and 'hs')
%
% Q            - PEQ-Quality: f/bandWidth, bandWidth = f2-f1,
%                f1=f/2/Q*(sqrt(4*Q^2+1)-1), f2=fm/2/Q*(sqrt(4*Q^2+1)+1)
%
% design       - 'butter': butterworth filter
%                'LR': Linkwitz-Riley filter
%
% peq_design   - 'sym'   : symmetrical,
%                'hpl'   : half pad loss,
%                'constQ': constant Q
%
% shelve_design -'low' : f specifies the freq. close to 0 dB gain
%                'mid' : f specifies the freq. at g/2 dB
%                'high': f specifies the freq. close to g dB
%
% warp_type    - Strings specifying the pre-warping applied to Q for the PEQ:
%                'sin' sinus pre-warping
%                'tan' tangens pre-warping
%                'cos' cosinus pre-warping
%                'none' no pre-warping
%                The effect of warp_type is very similar if f << fs/2
%
% frac         - fractional Octave used for filter bank
%                (e.g. frac=3: third octave, frac=1: octave filter bank)
%
% bw           - bandwidth relativ to center frequency (c.f. [1])
%
%
% O U T P U T
% y            - filtered input signal
% varOut1      - exact center frequencies, and lower, and upper cut-off,
%                and bandwidth of filters if 'fracOct', or 'FFTfracOct'
%
% REFERENCES
% [1] Brian C J Moore, Simon R Oldfield, Gary J Dooley: "Detection and
%     discrimination of spectral peaks and notches at 1 and 8 kHz," J.
%     Acoust. Soc. Am., vol. 85(2):820-836, (1989).
%
% frank.schultz@uni-rostock.de, Telecommunication engineering, University
% Rostock
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
%
% v. 1.0 03/2015 initial dev.

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [y, varOut1] = AKfilter(x, type, f, g, fs, varIn1, varIn2, varIn3)

% ---- default arguments and input check ----
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('varIn1', 'var')
    varIn1 = false;
end
if ~exist('varIn2', 'var')
    varIn2 = false;
end
if ~exist('varIn3', 'var')
    varIn3 = false;
end

type = upper(type);

% ---- reshape ----
if ~any(strcmp(type, {'BP' 'BS' 'FRACOCT' 'FFTFRACOCT' 'SLANEY1998'}))
    % single values for f
    f = reshape(f, numel(f), 1);
else
    % f contains upper and lower cut-off frequencies
    if size(f,2)~=2 && size(f,1)==2
        f = f';
    elseif size(f,2)~=2 && size(f,1)~=2
        error('AKfilter:Input', 'fc must have values for upper and lower cut-off frequency')
    end
end
g      = reshape(g, numel(g), 1);
fs     = reshape(fs, numel(fs), 1);
varIn1 = reshape(varIn1, numel(varIn1), 1);
if ~iscell(varIn2)
    tmp       = varIn2;
    clear varIn2
    varIn2{1} = tmp;
end
varIn2 = reshape(varIn2, numel(varIn2), 1);
if isnumeric(varIn3)
    varIn3 = reshape(varIn3, numel(varIn3), 1);
elseif ~iscell(varIn3)
    tmp       = varIn3;
    clear varIn3
    varIn3{1} = tmp;
end
varIn3 = reshape(varIn3, numel(varIn3), 1);
clear tmp

% ---- Check if a single or multiple filters are applied ----

% number of filters
N_f = max([size(f,1) numel(g) numel(fs) numel(varIn1) numel(varIn2) numel(varIn3)]);

if N_f == 1
    % apply the same filter to all channels of x
    singleFilter = true;
else
    % apply different filters to all channels of x
    singleFilter = false;
    
    % repeat to match the size
    if size(f,1)~=N_f
        f     = repmat(f(1,:), [N_f 1]);
    end
    g(end+1:N_f,1)      = g(1);
    fs(end+1:N_f,1)     = fs(1);
    varIn1(end+1:N_f,1) = varIn1(1);
    varIn2(end+1:N_f,:) = varIn2(1,:);
    varIn3(end+1:N_f,:) = varIn3(1,:);
    
    if size(x,2) == 1
        x = repmat(x, [1 N_f]);
    elseif size(x,2)~=N_f
        error('AKfilter:Input', 'x must be single channel or match the number of the remaining input Parameters')
    end
end

% number of channels
N_c = size(x, 2);
if N_c == 1
    singleChannel = true;
else
    singleChannel = false;
end

% output matrix
y = zeros(size(x,1), max(N_f, N_c));

% ---- design filter ----
switch type
    % --------------------------------------------------------------- LP/HP
    case {'LP' 'HP'}
        
        % input variables
        if varIn1
            N = varIn1;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end
        
        if iscell(varIn2)
            filterType = varIn2;
        else
            error('AKfilter:Input', 'varIn2 not specified')
        end
        
        for m = 1:N_f
            % Check filter order in case of Linkwitz-Riley filter
            if strcmpi(filterType{m}, 'LR')
                if ~mod(N(m),2)
                    N(m) = ceil(N(m)/2);
                else
                    error('AKfilter:Input', 'N must be even for LR-filters')
                end
            end
            
            % design filter
            if strcmp(type, 'LP')
                Hd = fdesign.lowpass('N,F3dB', N(m), f(m)/fs(m)*2);
            else
                Hd = fdesign.highpass('N,F3dB', N(m), f(m)/fs(m)*2);
            end
            
            if strcmpi(filterType{m}, 'LR')
                Hd = design(Hd, 'butter');
                Hd = dfilt.cascade(Hd, Hd);
            else
                Hd = design(Hd, filterType{m});
            end
            
            % check stability
            if ~isstable(Hd)
                error('AKfilter:Input', [filterType{m} '-filter not stable for N=' num2str(N(m)) ' and fc=' num2str(f(m))])
            end
            
            % apply
            if singleFilter && ~singleChannel
                y = filter(Hd, x) * 10^(g/20);
            else
                y(:,m) = filter(Hd, x(:,m)) * 10^(g(m)/20);
            end
            
        end
    % -------------------------------------------------------- BP/BS/X-OVER
    case {'BP' 'BS' 'XOVER', 'X-OVER'}
        
        % input variables
        if varIn1
            N = varIn1;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end
        
        if iscell(varIn2)
            filterType = varIn2;
        else
            error('AKfilter:Input', 'varIn2 not specified')
        end
        
        % modify Input in case of x-over filters
        if any(strcmp(type, {'XOVER' 'X-OVER'}))
            f = repmat(f, [1 2]);
            y = repmat(y, [1 1 2]);
        end
        
        for m = 1:N_f
            % Check filter order in case of Linkwitz-Riley filter
            if strcmpi(filterType{m}, 'LR')
                if ~mod(N(m),2)
                    N(m) = ceil(N(m)/2);
                else
                    error('AKfilter:Input', 'N must be even for LR-filters')
                end
            end
            
            % check if cut-off frequecies are sorted correctly
            if f(m,1) > f(m,2)
                f(m,:) = fliplr(f(m,:));
            end
            
            % design filter
            if strcmp(type, 'BP')
                Hd1 = fdesign.highpass('N,F3dB', N(m), f(m,1)/fs(m)*2);
                Hd2 = fdesign.lowpass ('N,F3dB', N(m), f(m,2)/fs(m)*2);
            else
                Hd1 = fdesign.lowpass ('N,F3dB', N(m), f(m,1)/fs(m)*2);
                Hd2 = fdesign.highpass('N,F3dB', N(m), f(m,2)/fs(m)*2);
            end
            
            if strcmpi(filterType{m}, 'LR')
                Hd1 = design(Hd1, 'butter');
                Hd1 = dfilt.cascade(Hd1, Hd1);
                Hd2 = design(Hd2, 'butter');
                Hd2 = dfilt.cascade(Hd2, Hd2);
            else
                Hd1 = design(Hd1, filterType{m});
                Hd2 = design(Hd2, filterType{m});
            end
            
            % check stability
            if ~isstable(Hd1)
                error('AKfilter:Input', [filterType{m} '-filter not stable for N=' num2str(N(m)) ' and fc=' num2str(f(m,1))])
            elseif ~isstable(Hd2)
                error('AKfilter:Input', [filterType{m} '-filter not stable for N=' num2str(N(m)) ' and fc=' num2str(f(m,2))])
            end
            
            % apply
            if singleFilter && ~singleChannel
                if any(strcmp(type, {'XOVER', 'X-OVER'}))
                    y(:,:,1) = filter(Hd1, x) * 10^(g/20);
                    y(:,:,2) = filter(Hd2, x) * 10^(g/20);
                    % invert phase in case butter filter order is not even and
                    % we have a LR x-over
                    if ~mod(N-1,2) && strcmpi(filterType, 'LR')
                        y(:,:,2) = -y(:,:,2);
                    end
                    % invert phase in case butter filter order is of 2:4:6
                    if ~mod(N-2,4) && strcmpi(filterType, 'BUTTER')
                        y(:,:,2) = -y(:,:,2);
                    end
                elseif strcmp(type, 'BS')
                    y = (filter(Hd1, x) + filter(Hd2, x)) * 10^(g/20);
                else
                    Hd1 = dfilt.cascade(Hd1, Hd2);
                    y = filter(Hd1, x) * 10^(g/20);
                end
            else
                if any(strcmp(type, {'XOVER', 'X-OVER'}))
                    y(:,m,1) = filter(Hd1, x(:,m)) * 10^(g(m)/20);
                    y(:,m,2) = filter(Hd2, x(:,m)) * 10^(g(m)/20);
                    % invert phase in case butter filter order is not even and
                    % we have a LR x-over
                    if ~mod(N(m)-1,2) && strcmpi(filterType{m}, 'LR')
                        y(:,m,2) = -y(:,m,2);
                    end
                    % invert phase in case butter filter order is of 2:4:6
                    if ~mod(N(m)-2,4) && strcmpi(filterType{m}, 'BUTTER')
                        y(:,m,2) = -y(:,m,2);
                    end
                elseif strcmp(type, 'BS')
                    y(:,m) = (filter(Hd1, x(:,m)) + filter(Hd2, x(:,m))) * 10^(g(m)/20);
                else
                    Hd1 = dfilt.cascade(Hd1, Hd2);
                    y(:,m) = filter(Hd1, x(:,m)) * 10^(g(m)/20);
                end 
            end
            
        end
    % ----------------------------------------------------------------- PEQ
    case 'PEQ'
        
        % input variables
        q         = varIn1;
        peq_type  = varIn2;
        warp_type = varIn3;
        
        for m = 1:N_f
            
            % check input
            if ~any(strcmpi(warp_type{m}, {'sin', 'cos', 'tan', 'none'}))
                error('AKfilter:Input', [warp_type{m} ' is not a valid warp_type'])
            end
            
            % prewarp
            q_warp = AKprewarpQ(f(m), fs(m), q(m), warp_type{m});
            
            if any(strcmpi(peq_type{m}, {'I' 'con' 'cons' 'const' 'constQ'}))
                peq_type{m} = 'I';
            elseif any(strcmpi(peq_type{m}, {'II' 'symm' 'sym' 'symmetric' 'symmetrical'}))
                peq_type{m} = 'II';
            elseif any(strcmpi(peq_type{m}, {'III' 'half' 'halfpadloss' 'hpl'}))
                peq_type{m} = 'III';
            else
                error('AKfilter:Input', [peq_type{m} ' is not a valid peq_type'])
            end
            
            % design PEQ
            [b, a] = AKpeq2(f(m), fs(m), g(m), q_warp, peq_type{m});
            
            % apply PEQ
            if singleFilter && ~singleChannel
                y      = filter(b, a, x);
            else
                y(:,m) = filter(b, a, x(:,m));
            end
            
        end
    % ------------------------------------------------ LOWSHELVE/HIGHSHELVE
    case {'LS' 'HS'}
        
        % input variables
        if varIn1
            N = varIn1;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end

        if iscell(varIn2)
            shelve_type = varIn2;
        else
            error('AKfilter:Input', 'varIn2 not specified')
        end
        
        
        for m = 1:N_f

            % set shelve type
            if any( strcmpi(shelve_type{m}, {'high' 'I'}) )
                shelve_type{m} = 'I';
            elseif any( strcmpi(shelve_type{m}, {'low' 'II'}) )
                shelve_type{m} = 'II';
            elseif any( strcmpi(shelve_type{m}, {'mid' 'III'}) )
                shelve_type{m} = 'III';
            else
                error('AKfilter:Input', [shelve_type{m} ' is not a valid shelve_type'])
            end
            
            % Butterworth like quality
            q_z = 1/sqrt(2);
            q_p = q_z;
            
             % get filter coefficients
            if strcmpi(type, 'LS')
                if N(m) == 1
                    [b, a] = AKlowshelve1(f(m), fs(m), g(m), shelve_type{m});
                elseif N(m) == 2
                    [b, a] = AKlowshelve2(f(m), fs(m), g(m), q_z, q_p, shelve_type{m});
                else
                    error('AKfilter:Input', 'N must be 1 or 2')
                end
            else
                if N(m) == 1
                    [b, a] = AKhighshelve1(f(m), fs(m), g(m), shelve_type{m});
                elseif N(m) == 2
                    [b, a] = AKhighshelve2(f(m), fs(m), g(m), q_z, q_p, shelve_type{m});
                else
                    error('AKfilter:Input', 'N must be 1 or 2')
                end
            end
            
            % apply filter
            if singleFilter && ~singleChannel
                y      = filter(b, a, x);
            else
                y(:,m) = filter(b, a, x(:,m));
            end
            
            
        end
        
    % --------------------------------------- Fractional octave filter bank
    case 'FRACOCT'
        
        % input variables
        if varIn1
            N = varIn1;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end
        if iscell(varIn2)
            filterType = varIn2;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end
        if varIn3
            frac    = varIn3;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end
        if numel(unique(f(:,1)))~=1 || numel(unique(f(:,2)))~=1 || numel(unique(frac))~=1
            error('AKfilter:Input', 'f must contain two values, and varIn3 must be scalar for fractional octave filter bank')
        end
        
        for m = 1:N_f
            
            % get fractional octave frequencies
            f_oct = AKfractionalOctaves(frac(m), f(m,:), fs(m));
            
            % number of fiters
            C = numel(f_oct(:,1));
            
            % allocate output signal
            if m==1
                y = zeros(size(x,1), size(x,2), C);
            end
            
            % Check filter order in case of Linkwitz-Riley filter
            if strcmpi(filterType{m}, 'LR')
                if ~mod(N(m),2)
                    N(m) = ceil(N(m)/2);
                else
                    error('AKfilter:Input', 'N must be even for LR-filters')
                end
            end
            
            for c = 1:C
                
                % design filter
                Hd1 = fdesign.highpass('N,F3dB', N(m), f_oct(c,2)/fs(m)*2);
                Hd2 = fdesign.lowpass ('N,F3dB', N(m), f_oct(c,3)/fs(m)*2);
                
                if strcmpi(filterType{m}, 'LR')
                    Hd1 = design(Hd1, 'butter');
                    Hd1 = dfilt.cascade(Hd1, Hd1);
                    Hd2 = design(Hd2, 'butter');
                    Hd2 = dfilt.cascade(Hd2, Hd2);
                else
                    Hd1 = design(Hd1, filterType{m});
                    Hd2 = design(Hd2, filterType{m});
                end
                
                % check stability
                if ~isstable(Hd1)
                    error('AKfilter:Input', [filterType{m} '-filter not stable for N=' num2str(N(m)) ' and fc=' num2str(f(m,1))])
                elseif ~isstable(Hd2)
                    error('AKfilter:Input', [filterType{m} '-filter not stable for N=' num2str(N(m)) ' and fc=' num2str(f(m,2))])
                end
                
                % cascade
                Hd1 = dfilt.cascade(Hd1, Hd2);
                % scale to 0 dB at center frequency
                Hf_oct = abs(freqz(Hd1, [f_oct(c,1) 100], fs(m))); % MUST pass two values for f - second one is discarded
                
                % apply
                if singleFilter && ~singleChannel
                    y(:,:,c) = filter(Hd1, x)      * 10^(g   /20) / Hf_oct(1);
                else
                    y(:,m,c) = filter(Hd1, x(:,m)) * 10^(g(m)/20) / Hf_oct(1);
                end
            end
            
            varOut1 = f_oct;
            
        end
        
    % ------------------------- FFT brickwall fractional octave filter bank
    case 'FFTFRACOCT'
        
        % input variables
        if varIn1
            frac = varIn1;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end
        
        if numel(unique(f(:,1)))~=1 || numel(unique(f(:,2)))~=1 || numel(unique(frac))~=1
            error('AKfilter:Input', 'f must contain two values, and varIn3 must be scalar for fractional octave filter bank')
        end
        
        for m = 1:N_c
            [in, ~, fc] = AKfftFilterbank(x(:,m), fs(m), frac(m), f(m,:));
            if m == 1
                y = zeros(size(in, 1), numel(fs), size(in,3));
            end
            y(:,m,:) = in * 10^(g(m)/20);
        end
        
        varOut1      = fc(:,4:6);
        varOut1(:,4) = varOut1(:,3) - varOut1(:,2);
        
    % ---------------------------------------------------------- Moore 1989
    case 'MOORE1989'
        
        if varIn1
            bw = varIn1;
        else
            error('AKfilter:Input', 'varIn1 not specified')
        end
        
        % filter length
        N = size(x,1);
        % get filter
        b = AKmoore1989(f, g, bw, fs, N, N_f);
        % apply filter
        y = fftfilt(b, x);
        
    otherwise
        error([type 'is not suported for input argument type'])
end
