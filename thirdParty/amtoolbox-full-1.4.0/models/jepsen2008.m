function [outsig, fc, mfc] = jepsen2008(insig, fs, varargin);
%JEPSEN2008   Auditory model from Jepsen et. al. 2008
%   Usage: [outsig, fc] = jepsen2008(insig,fs);
%          [outsig, fc] = jepsen2008(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%
%   *Warning:* This code cannot be verified. It has not been possible to
%   tell from the desciption in the original paper nor from personal
%   communication with the original authors what the correct parameter set
%   used for the model is. This code is kept here as a reminder of the
%   structure of the model, and may reappear in a future work if a
%   verified parameter set can be established. The status of this piece
%   of code is "not even wrong": http://en.wikipedia.org/wiki/Not_even_wrong.
%
%   JEPSEN2008(insig,fs) computes the internal representation of the signal insig
%   sampled with a frequency of fs Hz as described in Jepsen, Ewert and
%   Dau (2008).
%  
%   [outsig,fc]=JEPSEN2008(...) additionally returns the center frequencies of
%   the filter bank.
%
%   The full Jepsen et al. (2008) model consists of the following stages:
% 
%     1) a heaphone filter to simulate the effect of a standard set of
%        headphones.
%
%     2) a middle ear filter to simulate the effect of the middle ear, and
%        to convert to stapes movement.
%
%     3) lopezpoveda2001 - Dual resonance non-linear filterbank.
%
%     4) an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz.
%
%     5) an expansion stage
%
%     6) an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%     7) a modulation filterbank.
%
%   Any of the optinal parameters for LOPEZPOVEDA2001, IHCENVELOPE and
%   ADAPTLOOP may be optionally specified for this function. They will be
%   passed to the corresponding functions.
%
%   See also: lopezpoveda2001, ihcenvelope, adaptloop, modfilterbank, dau1997
%
%   References:
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am.,
%     124(1):422--438, 2008.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/jepsen2008.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Untrusted
%   #Author: Torsten Dau 
%   #Author: Morten Løve Jepsen 
%   #Author: Peter L. Søndergaard (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
% ------ Checking of input parameters ------------

error(['This code of this function is incorrect. Please see the description ' ...
       'in the help text.']);

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import={'lopezpoveda2001','ihcenvelope','adaptloop'};
definput.importdefaults={'jepsen2008'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

%% Headphone filter
hp_fir = headphonefilter(fs);
outsig = filter(hp_fir,1,insig);

%% lopezpoveda2001 and compensation for middle-ear
[outsig, fc] = lopezpoveda2001(outsig, fs, 'argimport',flags,keyvals);
outsig = gaindb(outsig,50);

%% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

%% Expansion stage
outsig = outsig.^2;

%% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

%% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);


