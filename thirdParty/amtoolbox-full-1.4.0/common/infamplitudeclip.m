function xout = infamplitudeclip(xin,options)
%INFAMPLITUDECLIP  Perform infinite amplitude clipping on signal
%   Usage:  xout = infamplitudeclip(xin);
%          xout = infamplitudeclip(xin,'norm');
%
%   INFAMPLITUDECLIP(xin) performs infinite amplitude clipping on the
%   input signal. This is a signal modification that preserves the
%   zero-crossings of the signal, but sets the amplitude to either +1 or
%   -1 depending on the sign of the signal. This type of modification was
%   used in "Licklider and Pollack, 1948".
%
%   INFAMPLITUDECLIP(xin,'norm') or INFAMPLITUDECLIP(xin,'rms') will do
%   the same, but scale the signal so as to preserve the l^2 of the
%   signal (its RMS value).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/infamplitudeclip.php


%   #Author : Peter Soendergaard (2009)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   REFERENCES:
%   J. Licklider and I. Pollack. Effects of differentiation, integration,
%   and infinite peak clipping upon the intelligibility of speech. The
%   Journal of the Acoustical Society of America, 20:42-52, 1948.

  
  error(nargchk(1,2,nargin));

  l2 = norm(xin);
  
  xout = sign(xin);
  
  if nargin>1
    switch(lower(options))
     case {'rms','norm'}
      xout=xout/norm(xout)*l2;
     otherwise
      error('Unknown option: %s',options);
    end;
    
  end;
  
  
  
  


