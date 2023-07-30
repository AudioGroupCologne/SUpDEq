function inoutsig = adaptloop(inoutsig,fs,varargin)
%ADAPTLOOP   Adaptation loops
%   Usage: outsig = adaptloop(insig,fs,limit,minspl,tau);
%          outsig = adaptloop(insig,fs,limit,minlvl,tau);
%          outsig = adaptloop(insig,fs,limit,minspl);
%          outsig = adaptloop(insig,fs,limit,minlvl);
%          outsig = adaptloop(insig,fs,limit);
%          outsig = adaptloop(insig,fs);
%
%   ADAPTLOOP(insig,fs,limit,minspl,tau) applies non-linear adaptation to an
%   input signal insig sampled at a sampling frequency of fs Hz. 
%   limit (in arbitrary units) is used to limit the overshoot of the output.
%   minspl determines the lowest audible SPL of the signal (in dB).
%   minlvl is minspl but expressed as a linear amplitude to be directly passed
%   to the core of adaptloop. In order to be recognized as a valid minlvl, 
%   it needs to be 0 < minlvl < 1.
%   tau is a vector with time constants involved in the adaptation loops. 
%   The number of adaptation loops is determined by the length of tau.
%
%   ADAPTLOOP(insig,fs,limit,minspl) does as above, but uses the values for
%   tau determined in Dau. et al (1996a).
%
%   ADAPTLOOP(insig,fs,limit) does as above with an minspl of 0 dB.
%
%   ADAPTLOOP(insig,fs) does as above with an overshoot limit of limit=10.
%
%   ADAPTLOOP takes the following flags at the end of the line of input
%   arguments:
%
%     'adt_dau1997'        Default. This consists of 5 adaptation loops with
%                      an overshoot limit of 10 and a minimum SPL of
%                      0 dB. The adaptation loops have an 
%                      exponential delay. 
%
%     'adt_dau1996'        This is as in adt_dau1997 but without any
%                      overshoot limiting.  
%
%     'adt_puschel1988'    This consists of 5 adaptation loops without
%                      overshoot limiting. The adapation loops have a linear spacing.
%
%     'adt_breebaart2001'  As 'adt_puschel1998'
%
%     'adt_relanoiborra2019' As 'adt_dau1997' but with minspl of -34 dB. 
%     
%     'dim',d          Do the computation along dimension d of the input. 
%
%   See also: auditoryfilterbank, lopezpoveda2001
%
%   Demos: demo_adaptloop
%
%   References:
%     H. Relaño-Iborra, J. Zaar, and T. Dau. A speech-based computational
%     auditory signal processing and perception model. J. Acoust. Soc. Am.,
%     146(5), 2019.
%     
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074--1088, August 2001.
%     
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     J. Acoust. Soc. Am., 99(6):3615--3622, 1996a.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. II. Spectral and temporal integration. J.
%     Acoust. Soc. Am., 102:2906--2919, 1997b.
%     
%     D. Pueschel. Prinzipien der zeitlichen Analyse beim Hoeren. PhD thesis,
%     Universitaet Goettingen, 1988.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/adaptloop.php


%   #Author: Stephan Ewert (1999-2004) and Morten L. Jepsen: Original version
%   #Author: Peter L. Søndergaard (2009-2013): adapted to AMT
%   #Author: Piotr Majdak (2013-2021): adapted to AMT 1.0
%   #Author: Alejandro Osses (2021): bug fixes for AMT 1.0.x
%   #Author: Piotr Majdak (2021): adaptation for AMT 1.1

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% ------ Checking of input parameters and default parameters ---------

if nargin<2
  error('Too few input parameters.');
end;

definput.import = {'adaptloop'};
definput.keyvals.dim=[];
[~,keyvals,limit,minspl,tau]  = ltfatarghelper({'limit','minspl','tau'},definput,varargin);

if isfield(keyvals,'minlvl'), 
  warning('minlvl depracated, use minspl instead'); 
  if ~isempty(keyvals.minlvl), minspl=keyvals.minlvl; end % for backwards compatibility in case somebody provides 'minlvl' instead 'minspl'.
end
  
if ~isnumeric(minspl) || ~isscalar(minspl)
  error('%s: minlvl must be a scalar either as linear amplitude (0<minlvl<1) or in dB re 10 �Pa (otherwise).',upper(mfilename));
end;

% 
if minspl>0 && minspl<1, 
  % minspl is provided as minlvl (i.e., a linear amplitude), use it as it is.
  minlvl_lin=minspl;
else
  % minspl is provided as SPL in dB re 10 �Pa (the reference level of adaptloop)
  % convert it to linear amplitude minlvl_lin
  minlvl_lin=scaletodbspl(minspl,[],100); 
end
% amt_disp(['minspl is ' num2str(minspl) ' thus minlvl_lin is ' num2str(minlvl_lin)],'debug');

if ~isnumeric(tau) || ~isvector(tau) || any(tau<=0)
  error('%s: tau must be a vector with positive values.',upper(mfilename));
end;

if ~isnumeric(limit) || ~isscalar(limit) 
  error('%s: "limit" must be a scalar.',upper(mfilename));
end;  

% Note that the implementation of the adaptation loops assumes
% that any SPL is re 10 �Pa, i.e., dboffset=100.
% This needs to be considered when SPLs (in dB) are provided/interpret
% with this implementation of the adaptation lookps.
 
[inoutsig,~,~,~,dim,permutedsize,order]=assert_sigreshape_pre(inoutsig,[],keyvals.dim,upper(mfilename));
  inoutsig=comp_adaptloop(inoutsig,fs,limit,minlvl_lin,tau);
inoutsig=assert_sigreshape_post(inoutsig,dim,permutedsize,order);



