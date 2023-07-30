function [outsig_mfb, fc_mod, outsig_afb, fc] = relanoiborra2019_featureextraction(insig, fs, varargin)
%RELANOIBORRA2019_FEATUREEXTRACTION Creates internal representation based on Relano-Iborra et al. (2019)
%   Usage: [out, varargout] = relanoiborra2019_featureextraction(insig, fs, varargin)
%          [out, varargout] = relanoiborra2019_featureextraction(insig, fs, flow, fhigh, varargin)
%
%
%   Input parameters:
%     insig      :  signal to be processed
%     fs         :  Sampling frequency
%     flow       :  lowest center frequency of auditory filterbank
%     fhigh      :  highest center frequency of auditory filterbank
%     sbj        :  subject profile for drnl definition
%
%   Output parameters:
%     out        : correlation metric structure inlcuding
%
%   The out struct contains the following fields:
%
%     .dint         correlation values for each modulation band
%
%     .dsegments    correlation values from each time window and mod. band.
%
%     .dfinal       final (averaged) correlation
%
%
%     This script builds the internal representations of the template and target signals
%     according to the CASP model (see references).
%     The code is based on previous versions of authors: Torsten Dau, Morten
%     Loeve Jepsen, Boris Kowalesky and Peter L. Soendergaard
%
%     The model has been optimized to work with speech signals, and the
%     preprocesing and variable names follow this principle. The model is
%     also designed to work with broadband signals. In order to avoid undesired
%     onset enhancements in the adaptation loops, the model expects to recive a
%     prepaned signal to initialize them.
%
%   References:
%     H. Relaño-Iborra, J. Zaar, and T. Dau. A speech-based computational
%     auditory signal processing and perception model. J. Acoust. Soc. Am.,
%     146(5), 2019.
%     
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am., 124(1),
%     2008.
%     
%
%   See also: relanoiborra2019
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/relanoiborra2019_featureextraction.php


%   #Author: Helia Relano Iborra (March 2019): v4.0 provided to the AMT team
%   #Author: Clara Hollomey (2021): adapted to the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Stats M-Signal

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Auditory filtering:
if isoctave
   warning(['Currently this model is only fully functional under MATLAB.']); 
end

definput.import={'relanoiborra2019'}; % load defaults from arg_relanoiborra2019
[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Calculate fc's
if flags.do_erbspace,
  fc = erbspace(kv.flow, kv.fhigh, kv.fcnt); % method from Relano-Iborra et al. (2019)
end
if flags.do_erbspacebw,
  fc = erbspacebw(kv.flow, kv.fhigh, kv.bwmul, kv.basef); % method from Osses et al. (2022)
end
nchannels = size(fc,2);

%% Outer- and middle-ear filtering
%C.H.============================================================================    
b_hp = headphonefilter(fs); % calc headphone filtercoeffs
b_me = middleearfilter(fs);

insig = filter(b_hp,1,insig); % Outer-ear filterring
insig = filter(b_me,1,insig); % middle-ear-ear filterring

% Pre-allocate memory for auditory bands
outsig_afb = zeros(length(insig),nchannels);

for n = 1:nchannels % Filter
    outsig_afb(:,n) = relanoiborra2019_drnl(insig',fc(n),fs, kv.subject); 
end
outsig_afb = outsig_afb * 10^(50/20); % Gain to compensate for the Outer/middle ear attenuation

%% Add noise to channels (freq specific):

if flags.do_internalnoise,
  int_noise_lvl = [-68.15;-73.05;-74.75;-75.25;-77.45;-77.75;-77.65;-78.15;...
      -79.15;-79.75;-79.25;-78.35;-77.75;-77.35;-77.05;-76.45;-75.85;-75.37;...
      -74.85;-74.25;-73.65;-73.65;-74.15;-75.35;-74.65;-75.15;-76.15;-73.05;...
      -74.55;-75.45;-76.45;-77.35;-77.95;-80.35;-79.75;-78.15;-78.85;-79.55;...
      -82.25;-79.25;-82.75;-81.45;-83.05;-85.05;-87.45;-90.15;-92.55;-91.55;...
      -90.85;-90.15;-89.45;-88.65;-89.65;-90.55;-90.95;-92.95;-94.45;-96.15;...
      -97.85;-99.65];   % Internal noise levels to account for NH thresholds (single channel)

  broadband_compensation = -3.7; % parameter to compensate for use of full filter bank (not needed for single-channel simulations) 

  int_noise_lvl = int_noise_lvl + broadband_compensation;

  for m =1:size(outsig_afb, 2) %For each channel
      %int_noise = wgn(size(outsig_afb, 1),1, int_noise_lvl(m));
      int_noise = randn(size(outsig_afb, 1),1);
      int_noise = scaletodbspl(int_noise, int_noise_lvl(m),100);
      outsig_afb(:, m) = outsig_afb(:, m) + int_noise;
  end
end

%% 'Haircell' envelope extraction
if flags.do_ihc,
  outsig_afb = ihcenvelope(outsig_afb, fs, 'ihc_relanoiborra2019');
end
% outsig_ihc = outsig_ihc * 10^(50/20); % Gain to compensate for the Outer/middle ear attenuation

%% Expansion (Auditory nerve)
if flags.do_an,
  outsig_afb = outsig_afb.^2;
  outsig_afb = adaptloop(outsig_afb, fs,'argimport',flags,kv); 
end

%% Modulation processing:

% set lowest mf as constant value. The multiplication by 0 is just an easy
% way to get an array of zeros of the correct size.
mflow = fc.*0;

mfhigh= 1500; % Set maximum modulation center freq. to 1.5 kHz (v1.0 - HRI)

[fc_mod, outsig_mfb] = relanoiborra2019_mfbtd(outsig_afb,mflow,mfhigh,1,fs);	% MFB incl 150 LP

outsig_mfb = local_mfbtdpp(outsig_mfb,fc_mod); % Post processing of modulation subbands

%%% Piotr and Clara: 
%     1. relanoiborra2019_mtbtd + local_mfbtdpp is equivalent to
%        modfilterbank.m. What you lost when you reversed this change is the 
%        limitation of mfc bands (fc/4). 
%     2. Please note that my previous code relanoiborra2019_debug.m not only used
%        modfilterbank.m but I also adapted the decision back-end to read the
%        cell-formatted outputs. It is 'incorrect' to believe that the limitation
%        of bands is a part of a decision back-end.
%     3. mfhigh = 1500 Hz does not limit the filter bank to mfc of 1500 Hz, 
%        actually it is still 1000 Hz... try and see!
%        (Note by Alejandro on 6 June 2021)

function out = local_mfbtdpp(in,f)
% mfbtdpp.m - post processing for modulation filterbank 'mfbtd.m' (Dau et al. 1997) output.
%				  Gets the real part for centerfrequencies <= 10 Hz and the absolute value
%				  otherwise.
%
% Usage: out = mfbtdpp(in,inf,fs)
%
% in    = input matrix from mfbtd.m.
%   f   = center frequency info vector
% fs	  = sampling rate in Hz
%
% [out1,out2, ...,outn] = output matrix
%
% copyright (c) 1999 Universitaet Oldenburg
% changed by MLJ 25. april 2007
% adapted by Piotr Majdak (2021) to the AMT
out=in;

for i=1:length(f) % v2 MJ 17. oct 2006
   if f(i) <= 10
          out(:,:,i) = 1*real(out(:,:,i));
   else
       
      out(:,:,i) = 1/sqrt(2)*abs(out(:,:,i));
   end
end




