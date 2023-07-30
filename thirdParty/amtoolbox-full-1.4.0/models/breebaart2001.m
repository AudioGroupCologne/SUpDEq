function [ei_map,fc,outsigl,outsigr] = breebaart2001(insig,fs,tau,ild,varargin)
%BREEBAART2001   Binaural masking level differences
%   Usage: [ei_map, fc, outsigl, outsigr] = breebaart2001(insig,fs);
%          [ei_map, fc, outsigl, outsigr] = breebaart2001(insig,fs,...);
%          [ei_map, fc] = breebaart2001(insig,fs);
%          [ei_map, fc] = breebaart2001(insig,fs,...);
%
%   Input parameters:
%        insig  : input acoustic signal.
%        fs     : sampling rate.
%        tau    : characteristic delay in seconds (positive: left is leading)
%        ild    : characteristic ILD in dB (positive: left is louder)
%
%   Output parameters:
%       ei_map  : EI-cell representation
%       outsigl : internal monaural representation of the left ear
%       outsigr : internal monarual respresnetaion of th right ear
%       fc      : center frequencies of the filterbank
%  
%   BREEBAART2001(insig,fs,tau,ild) computes the EI-cell
%   representation of the signal insig sampled with a frequency of fs*
%   Hz as described in Breebaart (2001) of the signal insig sampled with 
%   a frequency of fs Hz. The parameters tau and ild define the 
%   sensitivity of the EI-cell.
%
%   The input must have dimensions '[time x left/right channel x signal no]'.
% 
%   [ei_map,fc]=BREEBAART2001(...) additionally 
%   returns the center frequencies of the filter bank.
%
%   [ei_map,fc,ml,mr]=BREEBAART2001(...) additionally 
%   returns the center frequencies of the filter bank and the internal
%   monaural representations.
%  
%   The Breebaart 2001 model consists of the following stages:
%
%   1) an outer and middle ear transfer function
%   
%   2) a gammatone filter bank with 1-erb spaced filters.
%
%   3) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 770 Hz.
%
%   4) an adaptation stage modelling nerve adaptation by a cascade of 5
%      loops.
%
%   5) an excitation-inhibition (EI) cell model.
%
%   Parameters for AUDITORYFILTERBANK, IHCENVELOPE, ADAPTLOOP and
%   BREEBAART2001_EICELL can be passed at the end of the line of input arguments.
%
%   Examples
%   --------
%
%   The following code sets up a simple test example for the binaural output:
%
%     % Setup parameters
%     fs      = 44100;            % Sampling rate
%     T       = 0.3;              % Duration
%     Spl1    = 75;               % SPL of input signal 1
%     Spl2    = 75;               % SPL of input signal 2
%     rho     = 0;                % normalized correlation of signals
%     tau     = 0;
%     ild     = 0;
%
%     % Generate signals:
%     t  = [0:1/fs:T];
%     n1 = scaletodbspl(randn(length(t),1),Spl1);
%     n2 = scaletodbspl(randn(length(t),1),Spl2);
%     x1 = n1*sqrt((1+rho)/2) + n2*sqrt((1-rho)/2);
%     x2 = n1*sqrt((1+rho)/2) - n2*sqrt((1-rho)/2);
%
%     % Run the model and plot it
%     [ei_map,fc] = breebaart2001([x1,x2], fs, tau, ild);
%     plotfilterbank(ei_map,1,fc,fs,'audtick','lin');
%
%   See also: breebaart2001_eicell auditoryfilterbank ihcenvelope adaptloop breebaart2001_outmiddlefilter
%             modfilterbank data_breebaart2001 sig_breebaart2001 demo_breebaart2001
%             breebaart2001_centralproc amt_emuexp exp_osses2021 exp_breebaart2001
%             exp_osses2022 dietz2011 osses2021
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074--1088, August 2001.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/breebaart2001.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Peter L. Soendergaard (2011)
%   #Author: Martina Kreuzbichler (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

  
% ------ Checking of input parameters ------------

if nargin<4
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import = {'auditoryfilterbank','ihcenvelope','adaptloop','breebaart2001_eicell'};
definput.importdefaults={'fhigh',8000,'ihc_breebaart2001','adt_breebaart2001'};

[flags,keyvals,~,~,~]  = ltfatarghelper({'flow', 'fhigh', ...
                    'basef'},definput,varargin);

% ------ Checking of output parameters ------------
% to ensure backwards compatibility
do_mono = 0;
if nargout > 2
    do_mono = 1;
end

% ------ do the computation -------------------------

%% Outer- and middle ear transfer function
for earnumber = 1:2
    outsig(:,earnumber) = breebaart2001_outmiddlefilter(insig(:,earnumber),fs);
end

%% Apply the auditory filterbank
% only if spacing is 1 filter per ERB and order is 4
[outsig, fc] = auditoryfilterbank(outsig,fs,'argimport',flags,keyvals);

% % Spacing is two filters per ERB -> spacing in ERB is 0.5
% order = 3; % filter order
% spacing = 0.5;
% flow = 20;
% keyvals.basef = (exp(0.11)-1)/0.00437;
% fc = audspacebw(flow,fhigh,spacing,keyvals.basef,'erb');
% [bgt,agt] = gammatone(fc,fs,order,'complex');
% outsig = 2*real(ufilterbankz(bgt,agt,outsig)); 

%% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

%% non-linear adaptation loops, 
% lowpass filter for monaural output and  breebaart2001_eicell for binaural output

outsignal = adaptloop(outsig,fs,'argimport',flags,keyvals);

if do_mono == 1
    % Calculate filter coefficients for the 10 ms lowpass filter.
    % This filter places a pole /very/ close to the unit circle.
    mlp_a = exp(-(1/0.01)/fs);
    mlp_b = 1 - mlp_a;
    mlp_a = [1, -mlp_a];

    % Apply the low-pass modulation filter.
    outsig = filter(mlp_b,mlp_a,outsignal);
    outsigl = outsig(:,:,1);
    outsigr = outsig(:,:,2);
end

[siglen,nfreqchannels,~,nsignals] = size(outsignal);

ei_map = zeros(siglen, nfreqchannels, nsignals);
for k=1:nsignals
  for g=1:nfreqchannels
    ei_map(:,g,k) = breebaart2001_eicell(squeeze(outsignal(:,g,:,k)),fs,tau,ild,'rc_a',keyvals.rc_a);
  end
end


