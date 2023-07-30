function [sii,A,Z] = hauth2020_sii(E,N,T,I,G)
%HAUTH2020_SII calculates the SII according to ANSI S3.5-1997
%
%   Usage: 
%     sii = hauth2020_sii(E,N,T,I,G)
%
%   Input parameters:
%     E :   Speech Spectrum Level (Section 3.6 in the standard)
%     N :   Equivalent Noise Spectrum Level (Section 3.15 in the standard)
%     T :   Equivalent Hearing Threshold Level [dBHL] (Section 3.23 in the standard)
%     I :   Band Importance function (Section 3.1 in the standard)
%     G :   Insertion Gain [dB] (Section 3.28 in the standard)
%
%   Output parameters:
%     sii : Speech Intelligibility Index
%
%   HAUTH2020_SII calculates the Speech Intelligibility Index
%   as required by the model hauth2020 and according to ANSI S3.5-1997.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hauth2020_sii.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Christopher F. Hauth (2020)
%   #Author: Dr. Thomas Brand (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. handling of input arguments, setting of general parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5 
    help(mfilename); 
    return; 
end

if T == 0
    T = zeros(size(E)); % No hearing loss
elseif length(T) ~= length(E)
    error('Equivalent Hearing Threshold Level: Vector size incorrect');
end % if T == 0

if G == 0
    G = zeros(size(E)); % No insertion gain
elseif length(G) ~= length(E)
    error('Insertion Gain: Vector size incorrect');
end % if G == 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. get frequency bands and variables according to number of frequency bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch length(E)
    case 6
        freqb.name = 'octave'; %octave frequency bands
        freqb      = local_frequencybands(freqb,I);
    case 17
        freqb.name = 'equal_cb'; % equal-contributing critical frequency bands
        freqb      = local_frequencybands(freqb,I);
    case 18
        freqb.name = 'one_third_octave'; % one-third octave frequency bands
        freqb      = local_frequencybands(freqb,I);
    case 21
        freqb.name = 'critical_band'; % critical frequency bands
        freqb      = local_frequencybands(freqb,I);
    case 30
        freqb.name = 'gammatone'; % gammatone
        freqb      = local_frequencybands(freqb,I);
        G = 0;      % don't use insertion gain
    otherwise
        freqb.name = 'critical_band'; 
        freqb      = local_frequencybands(freqb,I);
        c_bw = freqb.bands(:,2)-freqb.bands(:,1);
        c_f0 = freqb.bands(:,3);
        c_th = freqb.bands(:,4);
        c_sp = freqb.bands(:,5);
        c_im = freqb.importance;
        clear freqb;
        freqb.name = 'gammatone'; 
        g_f0 = G(:,1);
        g_bw = G(:,2);
        freqb.importance = interp1(c_f0,c_im./c_bw,g_f0,'linear','extrap').*g_bw;
        freqb.importance(freqb.importance<0) = 0;
        freqb.importance = freqb.importance/sum(freqb.importance);
        freqb.bands(:,3) = g_f0(:);
        freqb.bands(:,1) = -g_bw/2+sqrt(g_bw.^2/4+g_f0.^2);
        freqb.bands(:,2) = +g_bw/2+sqrt(g_bw.^2/4+g_f0.^2);
        freqb.bands(:,4) = interp1(c_f0,c_th,g_f0,'linear','extrap');
        freqb.bands(:,5) = interp1(c_f0,c_sp,g_f0,'linear','extrap');
        G = 0;
        %error('Equivalent speech Spectrum level: Vector size must be 6,17,18 or 21 (or 30)');
end % switch length(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. calculate SII according to  ANSI S3.5-1997
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equivalent Speech Spectrum Level (5.1.3, Eq. 17): 
% speech spectrum level + insertion gain
E = E(:) + G(:);    

% Self-Speech Masking Spectrum (4.3.2.1 Eq. 5)
V = E - 24;

% maximum of equivalent noise spectrum level and self-speech masking
% spectrum level (4.3.2.2)
B = max(V,N(:)+G(:));

% Calculate slope parameter Ci for spread of masking (4.3.2.3 )
switch freqb.name
    case 'octave'
    case 'equal_cb'
    case 'one_third_octave'
        % Eq. 7
        C = 0.6 .* (B + 10 * log10(freqb.bands(:,2)) - 6.353) - 80;
    case 'critical_band'
        % Eq. 6
        C = 0.6.* (B + 10 * log10(freqb.bands(:,2) - freqb.bands(:,1))) - 80;
    case 'gammatone'
        C = 0.6.* (B + 10 * log10(freqb.bands(:,2) - freqb.bands(:,1))) - 80;
end % switch freqb.

% Initialize Equivalent Masking Spectrum Level (4.3.2.4) 
% (mod by rainerb, 2007-05-22)
Z = B(:);

% Calculate Equivalent Masking Spectrum Level (4.3.2.5)
switch freqb.name
    case 'octave'
    case 'equal_cb'
    case 'one_third_octave'
        % Eq. 9
        for i = 2:18    % for lowest band no effect of upwards spread of masking
            Z(i) = 10 * log10(10.^(0.1 * N(i)) + sum(10.^(0.1 * (B(1:(i-1)) + 3.32 .* ...
                C(1:(i-1)) .* log10(0.89 ...
                * freqb.bands(i,3) ./ freqb.bands(1:(i-1),3))))));
        end;
    case 'critical_band'
        % Eq. 8 
        for i = 2:21    % for lowest band no effect of upwards spread of masking
            Z(i) = 10 * log10(10.^(0.1 * N(i)) + sum(10.^(0.1 * (B(1:(i-1)) + 3.32 .* ...
                C(1:(i-1)) .* ...
                log10(freqb.bands(i,3) ./ ...
                freqb.bands(1:(i-1),2))))));
        end;
    case 'gammatone'
%% rainerb, 2008-06-19: spread of masking is already included in gammatone
end % switch freqb.

% Equivalent Internal Noise Spectrum Level (4.4 Eq. 10):
% reference internal noise spectrum level (4th column of freqb.bands) + dB HL
X = freqb.bands(:,4) + T(:);

% Disturbance Spectrum Level (4.5):
% maximum of equivalent internal noise spectrum level and eqivalent
% masking spectrum level
D = max(Z,X);

% Level Distortion Factor (4.6 Eq. 11)
L = 1 - (E - freqb.bands(:,5) - 10) ./ 160;
L = ((L <= 1) .* L ) + (L > 1);
L = 1; % How does distortion fit to using external noise for simulating the hearing threshold in BSIM?

% Band Audibility Function
% temporary variable K (4.7.1 Eq. 12)
K = (E - D + 15) / 30;
K = K .* (K > 0);               % limit K to [0,1]
K = ((K <= 1) .*K ) + (K > 1);  % limit K to [0,1]

%  band audibility function (7.7.2 Eq. 13):
% product of distortion factor and temporal variable
A = L .* K;

% Speech Intelligibility Index (4.8 Eq. 14)
S = sum(freqb.importance .* A);

sii = S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/////////////////////////////////////////////////////////////////////////////////%
%% private functions
%/////////////////////////////////////////////////////////////////////////////////%


function freqb = local_frequencybands(freqb,I)
%this function selects the frequency bands, center frequencies etc. needed for calculation

switch freqb.name
    case 'octave'
        freqb.bands          = ...
            [
            177     354     250     -3.9    34.75; ...  % center freq.,  reference internal noise spectrum level, standard speech spectrum level for normal speech
            353     707     500     -9.7    34.27; ...
            707     1414    1000    -12.5   25.01; ...
            1414    2828    2000    -17.7   17.32; ...
            2828    5657    4000    -25.9    9.33; ...
            5657    11314   8000    -7.1    1.13];
        freqb.importance     = [0.06 0.17 0.24 0.26 0.21 0.05]';
    case 'equal_cb'
        freqb.bands          = [];
        freqb.importance     = [];
    case 'one_third_octave'
        freqb.bands          = ...
            [
            143   179   160    0.60  32.41
            179   224   200   -1.70  34.48
            224   281   250   -3.90  34.75
            281   355   315   -6.10  33.98
            355   447   400   -8.20  34.59
            447   561   500   -9.70  34.27
            561   710   630  -10.80  32.06
            710   894   800  -11.90  28.30
            894  1118  1000  -12.50  25.01
            1118  1414  1250  -13.50  23.00
            1414  1789  1600  -15.40  20.15
            1789  2236  2000  -17.70  17.32
            2236  2806  2500  -21.20  13.18
            2806  3550  3150  -24.20  11.55
            3550  4472  4000  -25.90   9.33
            4472  5612  5000  -23.60   5.31
            5612  7099  6300  -15.80   2.59
            7099  8980  8000   -7.10   1.13
            ];
        % Band importance functions:
        % 1:	standard of Table 3
        % 2:	various nonsense syllable tests where most English phonemes occur equally often
        % 3:	CID-22
        % 4:	NU6
        % 5:	Diagnostic Rhyme test
        % 6:	short passages of easy reading material
        % 7:	SPIN
        % 8: Easy speech (Pavlovic 87, JASA 82,413-422, Table II)
        freqb.importance     = ...
            [
            % 				0.00083	0		0.0365	0.0168	0		0.0114	0	    0.0114
            % 				0.0095	0		0.0279	0.013	0.024	0.0153	0.0255	0.0153
            % 				0.015	0.0153	0.0405	0.0211	0.033	0.0179	0.0256	0.0179
            % 				0.0289	0.0284	0.05	0.0344	0.039	0.0558	0.036	0.0558
            % 				0.044	0.0363	0.053	0.0517	0.0571	0.0898	0.0362	0.0898
            % 				0.0578	0.0422	0.0518	0.0737	0.0691	0.0944	0.0514	0.0944
            % 				0.0653	0.0506	0.0514	0.0658	0.0781	0.0709	0.0616	0.0709
            % 				0.0711	0.0584	0.0575	0.0644	0.0751	0.066	0.077	0.0660
            % 				0.0818	0.0667	0.0717	0.0664	0.0781	0.0628	0.0718	0.0628
            % 				0.0844	0.0774	0.0873	0.0802	0.0811	0.0672	0.0718	0.0672
            % 				0.0882	0.0893	0.0902	0.0987	0.0961	0.0747	0.1075	0.0747
            % 				0.0898	0.1104	0.0938	0.1171	0.0901	0.0755	0.0921	0.0755
            % 				0.0868	0.112	0.0928	0.0932	0.0781	0.082	0.1026	0.0820
            % 				0.0844	0.0981	0.0678	0.0783	0.0691	0.0808	0.0922	0.0808
            % 				0.0771	0.0867	0.0498	0.0562	0.048	0.0483	0.0719	0.0483
            % 				0.0527	0.0728	0.0312	0.0337	0.033	0.0453	0.0461	0.0453
            % 				0.0364	0.0551	0.0215	0.0177	0.027	0.0274	0.0306	0.0274
            % 				0.0185	0		0.0253	0.0176	0.024	0.0145	0	    0.0145
            0.0008  0.0000  0.0365  0.0168  0.0000  0.0114  0.0000  0.0114  0.0000
            0.0095  0.0000  0.0279  0.0130  0.0240  0.0153  0.0255  0.0153  0.0000
            0.0150  0.0153  0.0405  0.0211  0.0330  0.0179  0.0256  0.0179  0.0000
            0.0289  0.0284  0.0500  0.0344  0.0390  0.0558  0.0360  0.0558  0.0000
            0.0440  0.0363  0.0530  0.0517  0.0571  0.0898  0.0362  0.0898  0.0000
            0.0578  0.0422  0.0518  0.0737  0.0691  0.0944  0.0514  0.0944  1.0000
            0.0653  0.0506  0.0514  0.0658  0.0781  0.0709  0.0616  0.0709  0.0000
            0.0711  0.0584  0.0575  0.0644  0.0751  0.0660  0.0770  0.0660  0.0000
            0.0818  0.0667  0.0717  0.0664  0.0781  0.0628  0.0718  0.0628  0.0000
            0.0844  0.0774  0.0873  0.0802  0.0811  0.0672  0.0718  0.0672  0.0000
            0.0882  0.0893  0.0902  0.0987  0.0961  0.0747  0.1075  0.0747  0.0000
            0.0898  0.1104  0.0938  0.1171  0.0901  0.0755  0.0921  0.0755  0.0000
            0.0868  0.1120  0.0928  0.0932  0.0781  0.0820  0.1026  0.0820  0.0000
            0.0844  0.0981  0.0678  0.0783  0.0691  0.0808  0.0922  0.0808  0.0000
            0.0771  0.0867  0.0498  0.0562  0.0480  0.0483  0.0719  0.0483  0.0000
            0.0527  0.0728  0.0312  0.0337  0.0330  0.0453  0.0461  0.0453  0.0000
            0.0364  0.0551  0.0215  0.0177  0.0270  0.0274  0.0306  0.0274  0.0000
            0.0185  0.0000  0.0253  0.0176  0.0240  0.0145  0.0000  0.0145  0.0000
            ];
        if 1 == max(size(I))
            freqb.importance = freqb.importance(:,I);
        else
            if 18 == length(I(:))
                freqb.importance = I(:)/sum(I(:));
            else
                error('wrong band importance function');
            end
        end
    case 'critical_band'
        freqb.bands          = ...
            [
            100   200   150    1.50  31.44;  % lower and upper band limit, center freq
            200   300   250   -3.90  34.75;  % [Hz], internal noise spectrum
            300   400   350   -7.20  34.14;  % level [dB] and Standard
            400   510   450   -8.90  34.58;  % speech spectrum level
            510   630   570  -10.30  33.17;  % at normal vocal effort
            630   770   700  -11.40  30.34;
            770   920   840  -12.00  27.59;
            920  1080  1000  -12.50  25.01;
            1080  1270  1170  -13.20  23.52;
            1270  1480  1370  -14.00  22.28;
            1480  1720  1600  -15.40  20.15;
            1720  2000  1850  -16.90  18.29;
            2000  2320  2150  -18.80  16.37;
            2320  2700  2500  -21.20  13.80;
            2700  3150  2900  -23.20  12.21;
            3150  3700  3400  -24.90  11.09;
            3700  4400  4000  -25.90   9.33;
            4400  5300  4800  -24.20   5.84;
            5300  6400  5800  -19.00   3.47;
            6400  7700  7000  -11.70   1.78;
            7700  9500  8500   -6.00  -0.14
            ];
        freqb.importance     = ...
            [
            0.0130;      % Band importance for SPIN
            0.0478;
            0.0451;
            0.0470;
            0.0523;
            0.0591;
            0.0591;
            0.0503;
            0.0503;
            0.0556;
            0.0699;
            0.0625;
            0.0602;
            0.0684;
            0.0638;
            0.0605;
            0.0534;
            0.0394;
            0.0291;
            0.0132;
            0.0000
            ];
    case 'gammatone'
        freqb.bands          = ...
            [
             126.33       166.81      146.02       1.7149      31.308     %  18.778        30    
             166.8        211.89      188.74      -0.59196     32.722     %  14.657        31    
             211.9        262.13      236.34      -3.1624      34.298     %  11.52         33.66 
             262.14       318.09      289.36      -5.1989      34.51      %   9.065        34.341
             318.09       380.43      348.42      -7.1479      34.15      %   7.0957       34.337
             380.43       449.87      414.21      -8.2916      34.423     %   5.43         34.371
             449.87       527.22      487.5       -9.3375      34.139     %   4.0375       34.4  
             527.23       613.4       569.15     -10.29        33.18      %   3.0553       33.49 
             613.4        709.39      660.1      -11.062       31.209     %   2.2938       31.978
             709.39       816.32      761.41     -11.663       29.134     %   1.9362       29.834
             816.32       935.44      874.27     -12.107       27.037     %   1.9485       27.753
             935.44       1068.1     1000        -12.5         25.01      %   2.2          25.655
            1068.1       1216        1140.1      -13.077       23.782     %   2.5922       24.211
            1216         1380.6      1296.1      -13.704       22.738     %   2.7157       22.995
            1380.6       1564.1      1469.9      -14.608       21.355     %   2.0205       21.775
            1564.1       1768.4      1663.5      -15.781       19.678     %   1.0715       20.176
            1768.4       1996        1879.2      -17.085       18.103     %  -0.38433      18.634
            1996         2249.6      2119.4      -18.606       16.566     %  -1.9165       16.952
            2249.6       2532.1      2387.1      -20.426       14.629     %  -3.5223       15.016
            2532.1       2846.8      2685.2      -22.126       13.064     %  -4.7983       13.308
            2846.8       3197.3      3017.3      -23.599       11.947     %  -5.8713       12.281
            3197.3       3587.8      3387.3      -24.857       11.118     %  -6.2442       11.329
            3587.8       4022.8      3799.4      -25.566        9.9184    %  -6.1472       10.347
            4022.8       4507.4      4258.6      -25.35         8.2019    %  -5.1692       8.9881
            4507.4       5047.2      4770        -24.264        5.9709    %  -3.328        7.065 
            5047.2       5648.5      5339.7      -21.394        4.5609    %  -0.40942      5.1768
            5648.5       6318.4      5974.4      -17.939        3.2244    %   3.4962       3.8291
            6318.4       7064.6      6681.4      -13.638        2.2287    %   7.2948       2.6833
            7064.6       7895.8      7469         -9.9178       1.1797    %  11.001        1.8332
            7895.8       8821.8      8346.3       -6.5841       0.056736  %  13.725        0     
            ];
        freqb.importance     = ...
            [
            0.0076    0.0047        0.0050
            0.0108    0.0120        0.0655
            0.0149    0.0217        0.0225
            0.0311    0.0263        0.0035
            0.0573    0.0283        0.0596
            0.0650    0.0304        0.0761
            0.0649    0.0334        0.0159
            0.0532    0.0377        0.0023
            0.0473    0.0411        0.0000
            0.0418    0.0440        0.0743
            0.0380    0.0451        0.0617
            0.0366    0.0419        0.0062
            0.0355    0.0406        0.0178
            0.0385    0.0438        0.0004
            0.0405    0.0509        0.0058
            0.0402    0.0562        0.0435
            0.0397    0.0502        0.0153
            0.0396    0.0488        0.0394
            0.0416    0.0518        0.0827
            0.0424    0.0513        0.0485
            0.0406    0.0473        0.0209
            0.0355    0.0434        0.0013
            0.0280    0.0382        0.0305
            0.0236    0.0320        0.0259
            0.0235    0.0244        0.0342
            0.0195    0.0208        0.0183
            0.0150    0.0162        0.0425
            0.0139    0.0108        0.0658
            0.0103    0.0058        0.0172
            0.0036    0.0010        0.0347
            ];
        if 1 == max(size(I))
            freqb.importance = freqb.importance(:,I);
        else
            if 30 == length(I(:))
                freqb.importance = I(:)/sum(I(:));
            else
                error('wrong band importance function');
            end
        end

end % switch fb
return


