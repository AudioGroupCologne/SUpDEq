function [outsig, fc] = lopezpoveda2001(insig,fs,varargin)
%LOPEZPOVEDA2001  Dual resonance non-linear (DRNL) filterbank
%   Usage: outsig = lopezpoveda2001(insig,fs);
%
%   LOPEZPOVEDA2001(insig,fs) computes the Dual Resonance Non-Linear (lopezpoveda2001)
%   filterbank of the input signal insig sampled at fs Hz with channels
%   specified by the center frequencies in fc. The lopezpoveda2001 is described in
%   the paper Lopez-Poveda and Meddis (2001). The lopezpoveda2001 models the basilar
%   membrane non-linearity.
%
%   This version of the lopezpoveda2001 incoorperate the middleear filter used in
%   Lopez-Poveda and Meddis (2001).
%
%   The lopezpoveda2001 takes a lot of parameters which vary over frequency. Such a
%   parameter is described by a 1 x2 vector [b a] and indicates
%   that the value of the parameter at frequency fc can be calculated by
%
%       10^(b+a*log10(fc));
%
%   The parameters are:
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 80 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 8000 Hz.
%
%     'basef',basef  Ensure that the frequency basef is a centre frequency
%                    in the filterbank. The default value of [] means
%                    no default.
%
%     'middleear'    Perform middleear filtering before the actual lopezpoveda2001
%                    is applied using the middleear filter specified in
%                    Lopez-Poveda and Meddis (2001), and compensate for
%                    the effect of the filter after lopezpoveda2001 filtering. 
%                    This is the default.
%
%     'no_middleear' No middleear filtering. Be carefull with this setting,
%                    as another scaling must then be perform to convert the
%                    input to stapes movement.
%  
%     'bothparts'    Compute both the linear and the non-linear path of
%                    the lopezpoveda2001. This is the default.
%
%     'linonly'      Compute only the linear path.
%
%     'nlinonly'     Compute only the non-linear path.
%
%     'lin_ngt',n    Number of cascaded gammatone filter in the linear
%                    part, default value is 2.
%
%     'lin_nlp',n    Number of cascaded lowpass filters in the linear
%                    part, default value is 4
%
%     'lin_gain',g   Gain in the linear part, default value is [4.20405 ...
%                    .47909].
%
%     'lin_fc',fc    Centre frequencies of the gammatone filters in the
%                    linear part. Default value is [-0.06762 1.01679].
%
%     'lin_bw',bw    Bandwidth of the gammatone filters in the linear
%                    part. Default value is [.03728  .78563]
%
%     'lin_lp_cutoff',c
%                    Cutoff frequency of the lowpass filters in the
%                    linear part. Default value is [-0.06762 1.01679 ]
%
%     'nlin_ngt_before',n
%                    Number of cascaded gammatone filters in the
%                    non-linear part before the broken stick
%                    non-linearity. Default value is 3.
%
%     'nlin_ngt_after',n
%                    Number of cascaded gammatone filters in the
%                    non-linear part after the broken stick
%                    non-linearity. The default value of [] means to use the 'before'
%                    value.
%
%     'nlin_nlp',n   Number of cascaded lowpass filters in the
%                    non-linear part. Default value is 3.
%
%     'nlin_fc_before',fc
%                    Center frequencies of the gammatone filters in the
%                    non-linear part before the broken stick
%                    non-linearity. Default value is [-0.05252 1.01650].
%
%     'nlin_fc_after',fc
%                    Center frequencies of the gammatone filters in the
%                    non-linear part after the broken stick
%                    non-linearity. The default value of [] means to use the 'before'
%                    value.
%
%     'nlin_bw_before',bw
%                    Bandwidth of the gammatone filters in the
%                    non-linear part before the broken stick
%                    non-linearity. Default value is [-0.03193 .77426 ].
%
%     'nlin_bw_after',w
%                    Bandwidth of the gammatone filters in the non-linear
%                    part after the broken stick non-linearity. The default
%                    value of [] means to use the 'before' value.
%
%     'nlin_lp_cutoff',c
%                    Cutoff frequency of the lowpass filters in the
%                    non-linear part. Default value is [-0.05252 1.01650 ].
%
%     'nlin_a',a     The a coefficient for the broken-stick non-linearity. Default
%                    value is [1.40298 .81916 ].
%
%     'nlin_b',b     The b coefficient for the broken-stick non-linearity. Default
%                    value is [1.61912 -.81867].
%
%     'nlin_c',c     The c coefficient for the broken-stick non-linearity. Default
%                    value is [-.60206 0].
%
%     'nlin_d',d     The d coefficient for the broken-stick non-linearity. Default
%                    value is 1.
%
%   The output from LOPEZPOVEDA2001 can be conveniently visualized using the
%   plotfilterbank function from LTFAT.
%
%   See also: middleearfilter, headphonefilter, adaptloop, dbspl, middleearfilter, 
%             data_lopezpoveda2001, demo_lopezpoveda2001, exp_lopezpoveda2001, 
%             lopezpoveda2001, baumgartner2013, dau1997, osses2021
%
% 
%   References:
%     E. Lopez-Poveda and R. Meddis. A human nonlinear cochlear filterbank.
%     J. Acoust. Soc. Am., 110:3107--3118, 2001.
%     
%     R. Meddis, L. O'Mard, and E. Lopez-Poveda. A computational algorithm
%     for computing nonlinear auditory frequency selectivity. J. Acoust. Soc.
%     Am., 109:2852--2861, 2001.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/lopezpoveda2001.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Qualified
%   #Requirements: M-Signal
%   #Author: Morten Løve Jepsen
%   #Author: Marton Marschall (2008): bugfixes
%   #Author: Katarina Egger (2011): adaptation to AMT
%   #Author: Peter L. Søndergaard: cleanup

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% In comparison to the original code, the gammatone filters are computed
% by convolving the filter coefficients, instead of performing multiple
% runs of 'filter'. The lowpass filtering is still performed by mutliple
% runs through 'filter', as the linear-part lowpass filter turned out to
% be unstable when the coefficients was convolved.

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

  
% Import the parameters from the arg_drnl.m function.
definput.import={'lopezpoveda2001'};

[flags,kv]=ltfatarghelper({'flow','fhigh'},definput,varargin);

% Expand the 'after' settings, if they are empty.

if isempty(kv.nlin_ngt_after)
  kv.nlin_ngt_after=kv.nlin_ngt_before;
end;

if isempty(kv.nlin_fc_after)
  kv.nlin_fc_after=kv.nlin_fc_before;
end;

if isempty(kv.nlin_bw_after)
  kv.nlin_bw_after=kv.nlin_bw_before;
end;


% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(kv.flow, kv.fhigh, kv.bwmul, kv.basef);

%% Convert the input to column vectors
nchannels  = length(fc);

% Change f to correct shape.
[insig,siglen,nsigs,wasrow,remembershape]=comp_sigreshape_pre(insig,'lopezpoveda2001',0);

outsig=zeros(siglen,nchannels,nsigs);

% The current implementation of the lopezpoveda2001 works with
% dboffset=100, so we must change to this setting.
% The output is always the same, so there is no need for changing back.
 
% Obtain the dboffset currently used.
dboffset=dbspl(1);

% Switch signal to the correct scaling.
insig=gaindb(insig,dboffset-100);

%% Apply the middleear filter
if flags.do_middleear
  me_fir = middleearfilter(fs);
  insig = filter(me_fir,1,insig);  
end;

if flags.do_jepsen2008
  me_fir = middleearfilter(fs,'jepsen2008');
  insig = filter(me_fir,1,insig);
end;

%% ---------------- main loop over center frequencies

% Handle the compression limiting in the broken-stick non-linearity
if ~isempty(kv.compresslimit)
  fclimit=min(fc,kv.compresslimit);
else
  fclimit=fc;
end;

% Sanity checking, some center frequencies may go above the Nyquest
% frequency
% Happens for lin_lp_cutoff

for ii=1:nchannels

  % -------- Setup channel dependant definitions -----------------

  lin_fc        = polfun(kv.lin_fc,fc(ii));
  lin_bw        = polfun(kv.lin_bw,fc(ii));
  lin_lp_cutoff = polfun(kv.lin_lp_cutoff,fc(ii));
  lin_gain      = polfun(kv.lin_gain,fc(ii));
  
  nlin_fc_before = polfun(kv.nlin_fc_before,fc(ii));
  nlin_fc_after  = polfun(kv.nlin_fc_after,fc(ii));
  
  nlin_bw_before = polfun(kv.nlin_bw_before,fc(ii));
  nlin_bw_after  = polfun(kv.nlin_bw_after,fc(ii));
  
  nlin_lp_cutoff = polfun(kv.nlin_lp_cutoff,fc(ii));
  
  % Expand a and b using the (possibly) limited fc
  nlin_a = polfun(kv.nlin_a,fclimit(ii));

  % b [(m/s)^(1-c)]
  nlin_b = polfun(kv.nlin_b,fclimit(ii));
      
  nlin_c = polfun(kv.nlin_c,fc(ii));
  
  % Compute gammatone coefficients for the linear stage
  [GTlin_b,GTlin_a] = gammatone(lin_fc,fs,kv.lin_ngt,lin_bw/audfiltbw(lin_fc),'classic');
  
  % Compute coefficients for the linear stage lowpass, use 2nd order
  % Butterworth.
  [LPlin_b,LPlin_a] = butter(2,lin_lp_cutoff/(fs/2));

  % Compute gammatone coefficients for the non-linear stage
  %[GTnlin_b_before,GTnlin_a_before] = coefGtDRNL(nlin_fc_before,nlin_bw_before,...
  %                                               kv.nlin_ngt_before,fs);
  [GTnlin_b_before,GTnlin_a_before] = gammatone(nlin_fc_before,fs,kv.nlin_ngt_before,...
                                                nlin_bw_before/audfiltbw(nlin_fc_before),'classic');

  [GTnlin_b_after,GTnlin_a_after] = gammatone(nlin_fc_after,fs,kv.nlin_ngt_after,...
                                                nlin_bw_after/audfiltbw(nlin_fc_after),'classic');

  % Compute coefficients for the non-linear stage lowpass, use 2nd order
  % Butterworth.
  [LPnlin_b,LPnlin_a] = butter(2,nlin_lp_cutoff/(fs/2));

  % -------------- linear part --------------------------------

  if flags.do_bothparts || flags.do_linonly
    % Apply linear gain
    y_lin = insig.*lin_gain; 
    
    % Gammatone filtering
    y_lin = filter(GTlin_b,GTlin_a,y_lin);    
    
    % Multiple LP filtering
    for jj=1:kv.lin_nlp
      y_lin = filter(LPlin_b,LPlin_a,y_lin);
    end;
  else
    y_lin = zeros(size(insig));
  end;

  % -------------- Non-linear part ------------------------------

  if flags.do_bothparts || flags.do_nlinonly
    
    % GT filtering before
    y_nlin = filter(GTnlin_b_before,GTnlin_a_before,insig);
    
    % Broken stick nonlinearity
    if kv.nlin_d~=1
      % Just to save some flops, make this optional.
      y_nlin = sign(y_nlin).*min(nlin_a*abs(y_nlin).^kv.nlin_d, ...
                                 nlin_b*(abs(y_nlin)).^nlin_c);
    else
      y_nlin = sign(y_nlin).*min(nlin_a*abs(y_nlin), ...
                                 nlin_b*(abs(y_nlin)).^nlin_c);
    end;
    
    % GT filtering after
    y_nlin = filter(GTnlin_b_after,GTnlin_a_after,y_nlin);
    
    % then LP filtering
    for jj=1:kv.nlin_nlp
      y_nlin = filter(LPnlin_b,LPnlin_a,y_nlin);
    end;
  else
    y_nlin = zeros(size(insig));
  end;
  
  outsig(:,ii,:) = reshape(y_lin + y_nlin,siglen,1,nsigs);    
    
end;
 
function outpar=polfun(par,fc)
  %outpar=10^(par(1)+par(2)*log10(fc));
  outpar=10^(par(1))*fc^par(2);



