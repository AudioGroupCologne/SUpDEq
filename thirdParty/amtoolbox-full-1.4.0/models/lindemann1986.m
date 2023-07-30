function [crosscorr,t,ild,cfreq] = lindemann1986(insig,fs,varargin)
%LINDEMANN1986 Binaural activity map based on cross-correlation
%   Usage: [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f,M_f,T_int,N_1)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f,M_f,T_int)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f,M_f)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s,w_f)
%          [crosscorr,t] = lindemann1986(insig,fs,c_s)
%          [crosscorr,t] = lindemann1986(insig,fs)
%
%   Input parameters:
%       insig       : binaural signal for which the cross-correlation
%                     should be calculated
%       fs          : sampling rate (Hz)
%
%   Output parameters:
%       crosscorr   : A matrix containing the cross-correlation signal
%                     for every frequency channel fc and every time step n.
%                     The format of this matrix is output(n,m,fc), where m*
%                     denotes the correlation (delay line) time step.
%       t           : time axis for the time steps n in crosscorr
%       ild         : interaural level difference (ILD) for every freqeuncy
%                     channel fc*
%       cfreq       : center frequencies of every frequency channel
%
%   LINDEMANN1986(insig,fs) calculates a binaural activity map for the given
%   insig using a cross-correlation (delay-line) mechanism. The calculation
%   is done for every frequency band in the range 5-40 Erb.
%
%   Lindemann has extended the delay line model of Jeffres (1948) by a
%   contralateral inhibition, which introduce the ILD to the model.  Also
%   monaural detectors were extended, to handle monaural signals (and some
%   stimuli with a split off of the lateralization image). Hess has
%   extended the output from the Lindemann model to a binaural activity map
%   dependent on time, by using a running cross-correlation function.
%   This has been done here by starting a new running cross-correlation
%   every time step T_int.  A detailed description of these cross-
%   correlation steps is given in the LINDEMANN1986_BINCORR function.
%
%   The steps of the binaural model to calculate the result are the
%   following:
%
%   1) The given stimulus is filtered using an erb bank to
%      get 36 frequency bands containing a stimulus waveform.
%
%   2) In a second step the auditory nerve is simulated by extracting the
%      envelope using a first order low pass filter with a cutoff frequency
%      of 800 Hz and half-wave rectification.
%
%   3) Calculation of the cross-correlation between the left and right
%      channel.  This is done using the model described in Lindemann
%      (1986a) and Hess (2007). These are extensions to the delay line model
%      of Jeffres (1948).
%
%   You may supply any flags or key/value pairs of the AUDITORYFILTERBANK,
%   IHCENVELOPE or LINDEMANN1986_BINCORR at the end of the line of input
%   arguments.
%
%   Examples:
%   ---------
%
%   This example shows how to the binaural activity map for one frequency
%   channel of the Lindemann binaural model for a sinusoid with a binaural
%   modulation rate of 2 Hz. :
%     
%     fs = 44100; % Sampling rate    
%     f = 500;    % Frequency of the sinusoid
%     mf = 2;     % Binaural modulation frequency
%
%     % Generate 1~s binaural modulated sinusoid
%     sig = sig_lindemann1986(f,mf,fs);
%
%     % Model parameter (Note: T_int (ms) should be a multiple of 1000/f == 2)
%     % Begin of the storage of the cross-correlation is set to 1, because we have a
%     % non-stationary signal
%
%     % Calculate binaural cross-correlation
%     [cc,t] = lindemann1986(sig,fs,'T_int',6);
%
%     % Plot frequency channel 11, due to round(freqtoerb(500))==11
%     plot_lindemann1986(cc,t,'fc',f);
%
%   See also: lindemann1986_bincorr, plot_lindemann1986, gammatone, ufilterbankz
%             ihcenvelope itd2angle itd2angle_lookuptable data_lindemann1986
%             sig_lindemann1986 demo_lindemann1986 lindemann1986_centroid
%             wierstorf2013_estimateazimuth exp_lindemann1986
%
%   Demos: demo_lindemann1986
%
%   References:
%     W. Gaik. Combined evaluation of interaural time and intensity
%     differences: Psychoacoustic results and computer modeling. J. Acoust.
%     Soc. Am., 94:98--110, 1993.
%     
%     W. Hess. Time-Variant Binaural-Activity Characteristics as Indicator of
%     Auditory Spatial Attributes. PhD thesis, Ruhr-Universitaet Bochum,
%     2007.
%     
%     L. Jeffress. A place theory of sound localization. Journal of
%     comparative and physiological psychology, 41(1):35--39, 1948.
%     
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608--1622, 1986.
%     
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. II. The law of the first wave front. J.
%     Acoust. Soc. Am., 80:1623--1630, 1986.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/lindemann1986.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Author: Wolfgang Hess (before 2021)
%   #Author: Hagen Wierstorf (2012)
%   #Author: Piotr Majdak (2017): various adaptations for the AMT
%   #Author: Clara Hollomey (2021): various adaptations for the AMT 1.0

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Checking of input  parameters ---------------------------------
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;
if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end
if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

% Parse the command line and load default parameters
% For default values see lindemann1986a page 1613
% NOTE: I modified the default value for T_int from 10 to 5.
definput.import={'auditoryfilterbank','ihcenvelope','lindemann1986_bincorr'};
% Highest and lowest frequency to use for the erbfilterbank (this gives us
% 36 frequency channels, channel 5-40)
definput.importdefaults = ...
    {'flow',erbtofreq(5),'fhigh',erbtofreq(40),'ihc_lindemann1986'};
[flags,keyvals,c_s,w_f,M_f,T_int,N_1]  = ...
    ltfatarghelper({'c_s','w_f','M_f','T_int','N_1'},definput,varargin);


%% ------ Computation ---------------------------------------------------

% ------ Erb Bank -------------------------------------------------------
% Apply the auditory filterbank
% NOTE: Lindemann uses a bandpass filterbank after Duifhuis (1972) and
% Blauert and Cobben (1978).
[inoutsig,cfreq] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

% ------ ILD ------------------------------------------------------------
% Calculate the interaural level difference (ILD) for every frequency channel
% NOTE: this was not part of the original Lindemann model
ild = dbspl(inoutsig(:,:,2))-dbspl(inoutsig(:,:,1));

% ------ Cross-correlation computation ---------------------------------
% Extract the envelope, apply a half-wave rectification and calculate a
% running cross-correlation for every given frequency band
% ------ Haircell simulation -------
% Half-wave rectification and envelope extraction
inoutsig = ihcenvelope(inoutsig,fs,'argimport',flags,keyvals);
% ------ Cross-correlation ------
% Calculate the cross-correlation after Lindemann (1986a).
[crosscorr,t] = lindemann1986_bincorr(inoutsig,fs,c_s,w_f,M_f,T_int,N_1);


