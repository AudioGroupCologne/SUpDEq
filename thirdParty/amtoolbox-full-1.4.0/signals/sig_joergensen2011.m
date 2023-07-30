function [ir,fs]=sig_joergensen2011(duration)
%SIG_JOERGENSEN2011  Return a simulated room impulse response for Joergensen models
%   Usage: ir=sig_joergensen2011(duration);
%          [ir,fs]=sig_joergensen2011(duration);
%
%   SIG_JOERGENSEN2011(duration) returns a simulated impulse
%   response of length duration, where the duration is measured in
%   seconds. Only four durations (T30) are possible to select: 0.4, 0.7,
%   1.3, and 2.3 seconds.
%
%   [ir,fs]=SIG_JOERGENSEN2011(duration) also returns the
%   sampling frequency, fs = 44100 Hz.
%
%   The impulse responses were created using the ODEON room acoustic
%   software version 10 (Christensen, 2009). The simulated room was shaped
%   like an auditorium with maximal dimensions of 28x16x10m
%   (length-width-height).
%
%   The source and the receiver were horizontally with a fixed distance of
%   5m, and placed approximately the center of the room. All surfaces had
%   the same coefficient, which was adjusted individually frequency such
%   that the room had similar reverberation (T30) in the octave bands from
%   63 to 8000 Hz.
%
%   The corresponding clarity (C50), defined as the ratio of the energy of
%   first 50 ms of the impulse response to the energy of the part, was 0.60,
%   -2.9, -6.6, and -8.0 dB for the four different lengths, respectively.
%
%   See also: joergensen2011 joergensen2013
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_joergensen2011.php


%   #Author: Peter L. Sondergaard

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

switch(duration)
 case {0.4,0.7,1.3,2.3}
   [ir,fs]=amt_load('sig_joergensen2011',['simulatedimpulseresponses_',num2str(duration),'s.wav']);
 otherwise
   error('%s: Unsupported duration.',upper(mfilename));  
end;


