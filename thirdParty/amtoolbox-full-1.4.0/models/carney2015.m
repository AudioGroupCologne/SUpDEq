function [ic_sout_BE,ic_sout_BS,cn_sout,kv] = carney2015(an_sout, BMF, fs, varargin)
%CARNEY2015 Brainstem processing
%
%   Usage: [ic_sout_BE,ic_sout_BS,cn_sout] = carney2015(an_sout, BMF, fs)
%     [ic_sout_BE,ic_sout_BS,cn_sout, keyvals] = carney2015(an_sout, BMF, fs)
%     [..] = carney2015(an_sout, BMF, fs, varargin)
%
%   Input parameters:
%     an_sout : auditory nerve output [time CF]
%     BMF     : best masking frequency (Hz)
%     fs      : sampling frequency (Hz)
%
%   Output parameters:
%     ic_sout_BE : output of the band-enhanced cell [time CF]
%     ic_sout_BS : output of the band-suppressed cell [time CF]
%     cn_sout    : cochlear nucleus output [time CF]
%     keyvals    : parameters used in calculations, as key-value pairs:
%                  keyvals.t_cn : time vector for plotting CN responses,
%                  see arg_carney2015 for the list of other parameters
%
%   This function implements the LPBR model from Carney et al. (2015), 
%   which is an extension of the SFIE model from Nelson and Carney (2004). As in 
%   the Nelson and Carney (2004) model, it calculates the output ic_sout_BE*
%   of the BP IC cell. Further, it also calculates the output ic_sout_BS*
%   which is the output of the band-supressive IC cell. Further, 
%   the output cn_sout of the CN cell is provided. 
%   
%
%   Additional input parameters:
%
%     'tau_ex_cn',texc   CN excitation time constant (in s), see Eq. 2 in Nelson et al. (2004)
%
%     'tau_inh_cn',tic   CN inhibition time constant (in s)
%
%     'cn_delay',D       CN disynaptic inhibition delay (in s)
%
%     'Sinh_cn',Sc       CN excitatory strength
%
%     'afamp_cn',ac      CN alpha function area --> changes RATE of output cell
%
%     'tau_ex_ic',texc   IC excitation time constant (in s), see Eq. 2 in Nelson et al. (2004)
%
%     'tau_inh_ic',tic   IC inhibition time constant (in s)
%
%     'ic_delay_inh',D   IC inhibition delay (in s)
%
%     'afamp_ic',ai      IC alpha function area --> changes RATE of output IC BE cell
%
%     'Sinh_ic',Si       IC inhibitory strength
%
%     'inh_str_bs',isb   IC inhibitory strength, BS cell, see Carney et al. (2015)
%
%     'tau_inh_bs',tib   IC inhibition, from BE to BS cell
%
%     'ic_delay_bs',idb  IC local delay from BE to BS cell (in s)
%
%     'Aex',ae           rate Scalar for BS cell
%
%
%   See also: demo_carney2015 carney2015_fitaudiogram carney2015_getalphanorm
%             carney2015_generateneurogram zilany2014 bruce2018 
%
%   References:
%     L. Carney, T. Li, and J. McDonough. Speech coding in the brain:
%     Representation of vowel formants by midbrain neurons tuned to sound
%     fluctuations. eNeuro, 20(2), 2015.
%     
%     P. C. Nelson and L. Carney. A phenomenological model of peripheral and
%     central neural responses to amplitude-modulated tones. J. Acoust. Soc.
%     Am., 116(4), 2004.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/carney2015.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Authors: University of Rochester (UR EAR) team
%   #Authors: Clara Hollomey (2020): integration in the AMT
%   #Authors: Piotr Majdak (2021): integration for the AMT 1.0
%   #Authors: Alejandro Osses (2021): extensions for the AMT 1.1

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

% CN MODEL PARAMETERS:
% Detault parameter options by Carney, Li, and McDonough (2015):
definput.keyvals.tau_ex_cn  = 0.5e-3; % CN exc time constant
definput.keyvals.tau_inh_cn = 2.0e-3; % CN inh time constant
definput.keyvals.cn_delay   = 1.0e-3; % "disynaptic inhibition delay" (all ANFs excitatory)

definput.keyvals.Sinh_cn = 0.6;       % re: excitatory strength == 1
definput.keyvals.afamp_cn = 1.5;      % alpha function area --> changes RATE of output cell

% IC MODEL PARAMETERS:  
% BMF-dependent SFIE parameters
definput.keyvals.tau_ex_ic  = 1/(10*BMF); %[0.00025 0.0005 0.001 0.002]; % Time constant excitation in seconds
definput.keyvals.tau_inh_ic = definput.keyvals.tau_ex_ic*1.5;	 % Time constant inhibition in seconds
definput.keyvals.ic_delay_inh = definput.keyvals.tau_ex_ic*2;% Delay of inhibition in seconds
definput.keyvals.afamp_ic = 1;     % alpha function area --> changes RATE of output IC BE cell 
definput.keyvals.Sinh_ic  = 0.9;   % inhibitory strength

% BS parameters (LPBR model from Carney et al., 2015)
definput.keyvals.inh_str_bs = 4;  
definput.keyvals.tau_inh_bs = definput.keyvals.tau_inh_ic; %1.0e-3; % relatively long inhibition, from BE to BS
definput.keyvals.ic_delay_bs = 1.0e-3;  % Delay from BE to BS cell (local)
definput.keyvals.Aex = 0.5; % Rate Scalar for BS cell; note that this is effectively multiplied by afamp_ic (for Table in eNeuro)

definput.flags.ic_hwr = {'ic_hwr','no_ic_hwr'}; % Provide optioon to disable the halve-way rectification of the CN and IC outputs. 

[flags,kv]  = ltfatarghelper({},definput,varargin);
% CN model:
% Generate frequency-domain equivalent of alpha functions
[B1, A1] = carney2015_getalphanorm(kv.tau_ex_cn, fs, 1);
[B2, A2] = carney2015_getalphanorm(kv.tau_inh_cn, fs, 1);

factor = kv.afamp_cn*(1/fs);
numCF=size(an_sout,2);

cn_ex  = factor*[filter(B1, A1, an_sout); zeros(fs*kv.cn_delay,numCF)];
cn_inh = factor*[zeros(fs*kv.cn_delay,numCF); kv.Sinh_cn*filter(B2, A2,an_sout)];

% final CN model response:
if flags.do_ic_hwr
    cn_sout = (((cn_ex-cn_inh) + abs(cn_ex-cn_inh))/2); % subtract inhibition from excitation and half-wave-rectify
else
    cn_sout = cn_ex-cn_inh; % No half-wave rectification
end
kv.t_cn = ((0:length(cn_sout)-1)/fs)'; % time vector for plotting CN responses

% IC Model: (SFIE; Bandpass MRF)
[B3, A3] = carney2015_getalphanorm(kv.tau_ex_ic, fs, 1);
[B4, A4] = carney2015_getalphanorm(kv.tau_inh_ic, fs, 1);

factor = kv.afamp_ic*(1/fs);

ic_lp_ex1  = factor*[(filter(B3, A3, cn_sout)); zeros(floor(fs*kv.ic_delay_inh),numCF)];
ic_lp_inh1 = factor*[zeros(floor(fs*kv.ic_delay_inh),numCF); kv.Sinh_ic*(filter(B4, A4, cn_sout))];

if flags.do_ic_hwr
    % Exc min inh and half wave rectification
    ic_sout_BE = (((ic_lp_ex1-ic_lp_inh1) + abs(ic_lp_ex1-ic_lp_inh1))/2); % half-wave rectified; standard SFIE model
else
    ic_sout_BE = ic_lp_ex1-ic_lp_inh1;
end

%  Band-suppressed cell (see Carney et al., 2015)
[B5, A5] = carney2015_getalphanorm(kv.tau_inh_bs, fs, 1);

ic_bs_ex = kv.Aex * [ic_lp_ex1; zeros(floor(fs*kv.ic_delay_bs),numCF)]; % add zeros at end to match lengths  
ic_bs_inh = [zeros(floor(fs*kv.ic_delay_bs),numCF); kv.Aex*kv.inh_str_bs*(1/fs)*(filter(B5, A5,ic_sout_BE))];

if flags.do_ic_hwr
    % Exc min inh and half wave rectification
    ic_sout_BS = (((ic_bs_ex-ic_bs_inh) + abs(ic_bs_ex-ic_bs_inh))/2); % half-wave rectified
else
    % Only exc min inh
    ic_sout_BS = ic_bs_ex-ic_bs_inh;
end


