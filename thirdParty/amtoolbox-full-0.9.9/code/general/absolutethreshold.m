function [t,table]=absolutethreshold(freq,varargin)
%ABSOLUTETHRESHOLD  Absolute threshold at specified frequencies
%   Usage:  t=absolutethreshold(freq);
%           t=absolutethreshold(freq,...);
%
%   ABSOLUTETHRESHOLD(freq) will return the absolute threshold of hearing
%   in dB SPL at the frequencies given in freq. The output will have the
%   same shape as the input.
%
%   [outvals,table]=ABSOLUTETHRESHOLD(...) additionally returns the
%   frequencies table defining the choosen standard. The first column of
%   the table contains frequencies in Hz, and the second column contains the
%   absolute threshold at the given frequency.
%
%   ABSOLUTETHRESHOLD takes the following optional parameters:
%
%     'iso226_2003'  Free field binaural as given by the ISO 226:2003
%                    standard.
%
%     'iso389_2005'  Diffuse field binaural as given by the ISO 389-7:2005
%                    standard.
%
%     'map'          The ISO 226:2003 standard coverted to minimal
%                    audible pressure using the method from Bentler &
%                    Pavlovic 1989.
%
%     'er3a'         Using ER-3A insert earphones. This is described in
%                    the ISO 389-2:1994(E) standard
%
%     'er2a'         Using ER-2A insert earphones. This is described in
%                    Han & Poulsen (1998)
%
%     'hda200'       Using HDA200 circumaural earphone. This is
%                    decribed in the ISO 389-8:2004 standard.
%
%   Absolute thresholds for the ER2A and HDA-200 are provided up to 16
%   kHz by the ISO 389-5:2006 standard.
%
%   The default is to use the 'iso226_2003' setting.
%
%   Demos: demo_absolutethreshold
%
%   References:
%     R. A. Bentler and C. V. Pavlovic. Transfer Functions and Correction
%     Factors used in Hearing Aid Evaluation and Research. Ear and Hearing,
%     10:58-63, 1989.
%     
%     L. Han and T. Poulsen. Equivalent threshold sound pressure levels for
%     Sennheiser HDA 200 earphone and Etymotic Research ER-2 insert earphone
%     in the frequency range 125 Hz to 16 kHz. Scandinavian Audiology,
%     27(2):105-112, 1998.
%     
%     ISO 226:2003. Acoustics - Normal equal-loudness-level contours.
%     International Organization for Standardization, Geneva, Switzerland,
%     2003.
%     
%     ISO 389-2:1994(E). Acoustics - Reference zero for the calibration of
%     audiometric equipment - Part 2: Reference equivalent threshold sound
%     pressure levels for pure tones and insert earphones. International
%     Organization for Standardization, Geneva, Switzerland, 1994.
%     
%     ISO 389-5:2006. Acoustics - Reference zero for the calibration of
%     audiometric equipment - Part 5: Reference equivalent threshold sound
%     pressure levels for pure tones in the frequency range 8 kHz to 16 kHz.
%     International Organization for Standardization, Geneva, Switzerland,
%     2006.
%     
%     ISO 389-8:2004. Acoustics - Reference zero for the calibration of
%     audiometric equipment - Part 8: Reference equivalent threshold sound
%     pressure levels for pure tones and circumaural earphones. International
%     Organization for Standardization, Geneva, Switzerland, 2004.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/absolutethreshold.php

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
  
% AUTHOR : Peter Søndergaard based on data collected by Claus Elberling
  
if nargin<1
  error('Too few input parameters.');
end;

definput.import = {'absolutethreshold'};
definput.importdefaults = {'iso226_2003'};

[flags,kv,table]  = ltfatarghelper({'table'},definput,varargin);

t=spline(table(:,1),table(:,2),freq);


