function inoutsig = ihcenvelope(inoutsig,fs,varargin)
%IHCENVELOPE   Inner hair cell envelope extration
%   Usage:  outsig=ihcenvelope(insig,fs,methodname);
%
%   IHCENVELOPE(insig,fs,methodname) extract the envelope of an input signal
%   insig sampled with a sampling frequency of fs Hz. The envelope
%   extraction is performed by half-wave rectification followed by low pass
%   filtering. This is a common model of the signal transduction of the
%   inner hair cells.
%
%   The parameter methodname describes the kind of low pass filtering to
%   use. The name refers to a set of papers where in this particular
%   method has been utilized or studied. The options are
%
%     'ihc_bernstein'  Compute the Hilbert envelope, compress the envelope
%                      by raising it to the power .2, combine the envelope
%                      with the original fine-structure, half-wave rectify it, 
%                      square it and low-pass filter it with a cut-off
%                      frequency of 425 Hz. This method is defined in
%                      Bernstein (1999). Note that this method includes both a
%                      compression and an expansion stage.
%
%     'ihc_breebaart2001'  Use a 5th order filter with a cut-off frequency of 770
%                      Hz. This method is given in Breebaart (2001). Page
%                      94 in Breebart's thesis.
%
%     'ihc_filter_order',n
%                      Filter order for the Breebaart filter, default: 5.
%
%     'ihc_dau1996'    Use a 1st-order Butterworth filter with a cut-off
%                      frequency of 1000 Hz. This method has been used in all
%                      models deriving from the original 1996 model by 
%                      Dau et. al. These models are mostly monaural in nature.
%
%     'ihc_king2019'   Use a 1st-order Butterworth filter with a cut-off
%                      frequency of 1500 Hz. 
%
%     'hilbert'        Use the Hilbert envelope instead of the half-wave
%                      rectification and low pass filtering. This is not a
%                      releastic model of the inner hair envelope extraction
%                      process, but the option is included for
%                      completeness. The Hilbert envelope was first suggested
%                      for signal analysis in Gabor (1946).
%
%     'ihc_lindemann'  Use a 1st order Butterworth filter with a cut-off
%                      frequency of 800 Hz. This method is defined in the
%                      Lindemann (1986a) paper.
%
%     'ihc_meddis'     Use the Meddis inner hair cell model.
%
%     'minlvl'         Set all values in the output equal to minlvl.
%                      This ensures that the output is non-negative and
%                      that further processing is not affected by
%                      unnaturally small values. The default value of []
%                      means to not do this.
%
%     'dim',d          Work along dimension d.
%
%
%   References:
%     L. Bernstein, S. van de Par, and C. Trahiotis. The normalized
%     interaural correlation: Accounting for NoSπ thresholds obtained with
%     Gaussian and l̈ow-noisem̈asking noise. J. Acoust. Soc. Am.,
%     106:870--876, 1999.
%     
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074--1088, August 2001.
%     
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     J. Acoust. Soc. Am., 99(6):3615--3622, 1996a.
%     
%     D. Gabor. Theory of communication. J. IEE, 93(26):429--457, 1946.
%     
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608--1622, 1986.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/ihcenvelope.php


%   #Author: The AMT Team (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  

% ------ Checking of input parameters --------------------------------

if nargin<2
  error('Too few input parameters.');
end;


if ~isnumeric(inoutsig)
  error('%s: The input signal must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import = {'ihcenvelope'};
definput.keyvals.dim=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% ------ Computation -------------------------------------------------

[inoutsig,siglen,dummy,nsigs,dim,permutedsize,order]=assert_sigreshape_pre(inoutsig,[],keyvals.dim, ...
                                                  upper(mfilename));

if flags.do_ihc_undefined
  error(['%s: you must supply a flag to designate the IHC model to ' ...
         'use.'],upper(mfilename));
end;

if flags.do_ihc_bernstein1999
  % The computational trick mentioned in the Bernstein paper is used
  % here: Instead of raising the envelope to power .23 and combine with its
  % TFS, we raise it to power -.77, and combine with the original
  % signal. In this way we avoid computing the fine structure.
  inoutsig=max(abs(hilbert(inoutsig)).^(-.77).*inoutsig,0).^2;
  cutofffreq=425;
  [b, a] = butter(2, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
end;

if flags.do_ihc_breebaart2001
  inoutsig = max( inoutsig, 0 );
  % due to the successive application of the filter, the given 2000 Hz
  % correspond to a cut off-frequency of 770 Hz after the five iterations
  cutofffreq=2000;
  [b, a] = butter(1, cutofffreq*2/fs);
  for ii=1:keyvals.ihc_filter_order
    inoutsig = filter(b,a, inoutsig);
  end;
end;

if flags.do_ihc_dau1996
  inoutsig = max( inoutsig, 0 );
  cutofffreq=1000;
  [b, a] = butter(1, cutofffreq*2/fs);
  for ii=1:keyvals.ihc_filter_order
    inoutsig = filter(b,a, inoutsig);
  end
end;

if flags.do_hilbert
  inoutsig = abs(hilbert(inoutsig));
end;

if flags.do_ihc_king2019
    inoutsig = max( inoutsig, 0 );
    cutofffreq=1000;
    [b, a] = butter(1, cutofffreq*2/fs);
    for ii=1:keyvals.ihc_filter_order
        inoutsig = filter(b,a, inoutsig);
    end
end;

if flags.do_ihc_lindemann1986
  inoutsig = max( inoutsig, 0 );
  cutofffreq=800;
  [b, a] = butter(1, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
end;

if flags.do_ihc_meddis1990
  inoutsig = comp_meddishaircell(inoutsig, fs);
end;

if flags.do_ihc_relanoiborra2019
  % Calculate filter coefficients for the low-pass filtering following
  % half-wave rectification. 2nd order butterworth
  cutofffreq=1000;
  [b, a] = butter(2, cutofffreq*2/fs);
  % 'haircell' envelope extraction. Part 1: Half-wave rectification
  inoutsig = max( inoutsig, 0 );
  % 'haircell' envelope extraction. Part 2: Low-pass filtering
  inoutsig = filter(b,a, inoutsig);   
end;

if ~isempty(keyvals.ihc_minlvl)
  inoutsig = max( inoutsig, keyvals.ihc_minlvl );
end;

inoutsig=assert_sigreshape_post(inoutsig,dim,permutedsize,order);


