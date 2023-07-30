function outsig = sig_transposedtone(siglen,fc,fm,fs,varargin)
%SIG_TRANSPOSEDTONE  Transposed tone test stimuli
%   Usage:  ts = sig_transposedtone(fc,fm,dur,fs);
%
%   Input parameters:
%     siglen  : Length of signal
%     fc      : Vector of carrier frequencies (Hz)
%     fm      : Vector of modulation frequencies (Hz)
%     fs      : Sampling frequency (Hz)
% 
%   Output parameters:
%     outsig  : transposed tone (column vector)
%
%   SIG_TRANSPOSEDTONE(siglen,fc,fm,dur,fs) generates a transposed tone test
%   stimuli as defined in Kolrausch et. al (1997).
%
%   By default, the output is normalized to have an RMS value of 1, but this
%   can be changed by passing any of the flags from the normalize function.
% 
%   Some example parameters as used in a study by Santurette:
%   
%     siglen = 44100;
%     fc     = 5000;
%     fm     = 435;
%     fs     = 44100;
%     outsig = sig_transposedtone(fc,fm,dur,fs);
%
%   References:
%     A. Kohlrausch, R. Fassel, M. Heijden, R. Kortekaas, S. Par, A. Oxenham,
%     and D. Puschel. Detection of tones in low-noise noise: Further evidence
%     for the role of envelope fluctuations. Acta Acustica united with
%     Acoustica, 83(4):659--669, 1997.
%     
%     A. Oxenham, J. Bernstein, and H. Penagos. Correct tonotopic
%     representation is necessary for complex pitch perception. Proceedings
%     of the National Academy of Sciences, 101(5):1421--1425, 2004.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_transposedtone.php


%   #Author: Sébastien Santurette (2009)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin<4
  error('Too few input parameters.');
end;

definput.import={'normalize'};
definput.importdefaults={'rms'};
[flags,keyvals]=ltfatarghelper({},definput,varargin);


% time vector
t = (0:siglen-1)/fs;

% number of tones
n = length(fc);
                       
outsig = zeros(siglen,1);
for ii = 1:n
  % carrier tone
  carrier = cos(2*pi*fc(ii)*t);        

  % "modulator"   
  hrsine = sin(2*pi*fm(ii)*t);         

  % half-wave rectify the modulator
  hrsine = max(0,hrsine);             
  
  % Compute coefficients for low-pass filtering
  fcuth = .2*fc(ii)/(fs/2);
  [b,a] = butter(4,fcuth,'low');

  % Low-pass filter the modulator
  hrsine = filter(b,a,hrsine);

  % Add the carrier*modulator to the transposed tone
  outsig = outsig+(carrier.*hrsine)';
end;

outsig=normalize(outsig,flags.norm);



