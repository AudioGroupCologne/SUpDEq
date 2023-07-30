function [crosscorr,t] = lindemann1986_bincorr(insig,fs,varargin)
%LINDEMANN1986_BINCORR Cross-correlation between two input signals a la Lindemann
%   Usage: crosscorr = bincorr(insig,fs,c_s,w_f,M_f,T_int,N_1)
%
%   Input parameters:
%     insig      : signal with t xnfc x2 (right, left channel)
%     fs         : sampling rate
%
%   Output parameters:
%     crosscorr  : output matrix containing the correlations                     
%     t          : time axis corresponding to the n time samples in crosscorr
%
%   LINDEMANN1986_BINCORR(insig,fs) is an implementation of the Lindemann
%   cross-correlation algorithm to simulate a binaural delay line. The
%   output is a matrix of size n xm xamplitude, where
%   n = length(t)/fs.
%
%   The cross-correlation is calculated using:
%
%                  t
%     CC(tau,t) = int R(l-tau/2) * L(k+tau/2) exp(-(t-k)/T_int) dk
%                -inf
%
%   where T_{int} denotes an integration time constant and R, L the right and 
%   left input signal.
%
%   LINDEMANN1986_BINCORR takes the following key/value pairs at the end of
%   the command line:
%
%     'c_s',c_s     Stationary inhibition factor, 0 <= c_s <= 1 
%                   (0.0 = no inhibition). Default value is 0.3.
%
%     'w_f',w_f     Monaural sensitivity at the end of the delay line, 
%                   0 <= w_f < 1. Default value is 0.035.
%
%     'M_f',M_f     Determines the decrease of the monaural sensitivity along 
%                   the delay line. Default value is 6
%
%     'T_int',t_int
%                   integration time window (ms). This is the memory of the 
%                   correlation process with exp(-1/T_int). Also this
%                   determines the time steps in the binaural activity map,
%                   because every time step T_int a new running
%                   cross-correlation is started, so every T_int we have a new
%                   result in crosscorr. You can set T_int = inf if you like
%                   to have no memory effects, then you will get only one
%                   time step in crosscorr. Default value is 5 ms.
%
%     'N_1',N_1     Sample at which the first running cross-correlation should
%                   be started to avoid onset effects (see Lindemann (1986a) p.
%                   1614). Default: 1, 17640 (*200/f*fs with f=500 and fs=44100) 
%                   for 'stationary' (see above)  
%
%     'stationary'  will set the default values of N_1=17640 and T_int=inf, use
%                   this for stationary input signals.
%
%   The key/values can also be specified first in the line of arguments
%   in the following order: c_s, w_f, M_f, T_int, N_1.
% 
%   See also: lindemann1986
%
%   References:
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608--1622, 1986.
%     
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. II. The law of the first wave front. J.
%     Acoust. Soc. Am., 80:1623--1630, 1986.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lindemann1986_bincorr.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Author: Hagen Wierstorf (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%
% --- Used abbreviations ---
% 
%   l,L    : left signal
%   r,R    : right signal
%   m      : discrete time steps on the delay line (tau)
%   n      : discrete time (t)
%   M      : number of discrete steps on the delay line, 
%             length(delay line) = length(-M:M)
%   w_r    : monaural sensitivity in dependence of l
%   w_l    : monaural sensitivity in dependence of r
%   c_s    : stationary inhibition factor
%   w_f    : monaural sensitivity at the end of the delay line
%   M_f    : decrease of monaural sensitivity along the delay line
%   T_int  : integration time window (see above)
%   N_1    : lower summation boundary for the running cross-correlation
%   N_2    : upper summation boundary for the running cross-correlation
%   cc     : running cross-correlation
% 

%% ------ Checking of input parameters -----------------------------------
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'lindemann1986_bincorr'};
[flags,keyvals,c_s,w_f,M_f,T_int,N_1]  = ...
    ltfatarghelper({'c_s','w_f','M_f','T_int','N_1'},definput,varargin);


if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

if ~isnumeric(c_s) || ~isscalar(c_s) || c_s<0 || c_s>1
    error('%s: 0 <= c_s <= 1',upper(mfilename));
end

if ~isnumeric(w_f) || ~isscalar(w_f) || w_f<0 || w_f>=1
    error('%s: 0 <= w_f < 1',upper(mfilename));
end

if ~isnumeric(M_f) || ~isscalar(M_f) || M_f<=0
    error('%s: M_f has to be a positive scalar!',upper(mfilename));
end

if ~isnumeric(T_int) || ~isscalar(T_int) || T_int<=0
    error('%s: T_int has to be a positive scalar!',upper(mfilename));
end

if ~isnumeric(N_1) || ~isscalar(N_1) || N_1<=0
    error('%s: N_1 has to be a positive scalar!',upper(mfilename));
end


%% ------ Computation ---------------------------------------------------- 

siglen = size(insig,1);
nfc = size(insig,2);   % number of frequency channels

% Ensure 0 <= insig <= 1, so that 0 <= r,l <= 1 (see lindemann1986a, eq. 4)
insig = insig ./ (max(insig(:))+eps);

% Integration time of summing cross-correlation
T_int = round(T_int/1000 * fs);

% Check if the given signal length is long enough to provide an output for the
% given start N_1 of the running cross-correlation (see lindemann1986a, p. 1614)
if T_int<Inf && siglen<=N_1+T_int
    error('%s: siglen has to be longer than N_1+T_int==%i', ...
        upper(mfilename),N_1+T_int);
elseif T_int==Inf && siglen<=N_1
    error('%s: siglen has to be longer than N_1==%i for T_int==Inf', ...
        upper(mfilename),N_1);
end

% ------ Time steps on the delay line ------------------------------------
% -M:M are the time steps of the delay.
% Maximum delay (in samples): 1ms (this is for 500 Hz pure tones (T/2). In 
% common this should be different for every frequency band, see 
% Lindemann (1986a) S. 1613, eq. 21.)
% NOTE:
% The delay-line needs a sampling rate of fs*2. Therefore the signals are
% not doubled and filled with 0 as in Lindemann (1986a), page 1610 and 
% 1611, but the delay line time is halfed. If the signals are really filled with
% zeros instead of doubling every entry Lindemanns inhibition process won't work
% anymore!
M = round(fs/2 / 1000);
% Length of the delay line
ndl = length(-M:M);


% ------ Monaural sensitives of the correlator ---------------------------
% The following equations are from Lindemann (1986a) equation 9 on page
%>1611. Obviously is w_r(m) = w_l(-m).

% Monaural sensitivities for the left ear
w_r = w_f .* exp(-(2*M:-1:0)./M_f);
w_r = w_r';
% Duplicate columns of w_r for every band
w_r = w_r(:,ones(1,nfc));

% Monaural sensitivities for the right ear
w_l = w_f .* exp(-(0:2*M)./M_f);
w_l = w_l';
% Duplicate columns of w_l for every band
w_l = w_l(:,ones(1,nfc));


% ------ Calculate the cross-correlation ---------------------------------

% Prepare the left and right delay line signal
l = zeros(ndl,nfc);
r = zeros(ndl,nfc);


% Set upper summation index for running cross-correaltion
% See lindemann1986a, eq. 24
N_2 = setN_2(N_1,T_int,siglen);

% Generate time axis
t = (N_2:(N_2-N_1):siglen)'/fs;

% Memory preallocation
crosscorr = zeros( floor( (siglen-N_1)/(N_2-N_1) ),ndl,nfc );
cc = zeros(ndl,nfc);

ii = 1; % crosscorr index
for n = 1:siglen
    
    % ------ Inhibition --------------------------------------------------
    % Stationary inhibition after Lindemann (1986a, p. 1612 eq. 11a):
    % r(m+1,n+1) = r(m,n) * [1 - c_s l(m,n)]
    % l(m-1,n+1) = l(m,n) * [1 - c_s r(m,n)]
    % The length of r and l are the same as the length of the delay line
    % (ndl = length(-M:M)). Also the signals are mirror-inverted, 
    % because of their different directions passing the delay line. l starts
    % on the right sight and r on the left side of the delay line. If you 
    % want to mirror the axes of the delay line, you have to exchange the
    % input direction of r and l into the delay line.
    r_mn = r;
    l_mn = l;
    r = [ insig(n,:,2); r_mn(1:ndl-1,:) .* (1-c_s.*l_mn(1:ndl-1,:)) ];
    l = [ l_mn(2:ndl,:) .* (1-c_s.*r_mn(2:ndl,:)); insig(n,:,1) ];
    % TODO: the spectral shape is also needed for some application, so we should
    % check if the following solution can be fixed to work as it should.
    % The normalization with max([r l]+eps) is done, to avoid a
    % dependence on the amplitude for the inhibition mechanism (see Gaik
    % 1993, p. 109). This preserves the spectral shape across frequency
    % channels. Therefore the inhibition depends now only on its stationary 
    % inhibition factor c_s.
    %l = [ l(2:ndl,:) .* (1 - (c_s .* r(2:ndl,:) ./ ...
	%	      max(max([l r]+eps)) ) ); insig(n,:,1) ];
    %r = [ insig(n,:,2);   r(1:ndl-1,:) .* ...
	%	      (1 - (c_s .* l(1:ndl-1,:) ./ max(max([l r]+eps))) ) ];
 
    % ------ Monaural sensitivity and trading ----------------------------
    %
    % Monaural sensitivities after Lindemann (1986a, p. 1611 eq. 6a + 6b):
    % r'(m,n) = r(m,n) * [1 - w_l(m)] + w_l(m)
    % l'(m,n) = l(m,n) * [1 - w_r(m)] + w_r(m)
    %
    % Trading factors after Gaik (1993, p. 106).
    %>Multiplication by a level factor to reach an ILD of 0 for a given
    %>ITD, if this ILD corresponded to a natural ILD for this ITD.
    %
    %NOTE: trading will be implemented later
    %R(m) = r(m) .* trading(m,1) .* (1-w_l(m))  +  w_l(m);
    %L(m) = l(m) .* trading(m,2)  .* (1-w_r(m))  +  w_r(m);
    R = r .* (1-w_l)  +  w_l;
    L = l .* (1-w_r)  +  w_r;
    
    % ------ Cross-correlation -------------------------------------------
    % Calculate running cross-correlation (e.g. Lindemann 1986a, eq. 10 + 
    % eq. 24). 
    %              __
    %             \   N_2
    %   cc(m,n) = /__ n=N_1  R(m,n) * L(m,n) * exp(-(N_2-n)/T_int)
    %
    % NOTE: For simplicity with stationary signals Lindemann uses only this
    % formula which can be managed by setting T_int=Inf:
    %              __                           
    %             \   N_2                       
    %   cc(m,n) = /__ n=N_1  R(m,n) * L(m,n)    
    %                                          
    if n>=N_1 && n<=N_2
        cc = cc + R.*L .* exp( -(N_2-n) / T_int );
    end
    % If we have reached the upper summation index, store the result in
    % crosscorr and start a new summation
    if n==N_2
        crosscorr(ii,:,:) = cc;
        cc = zeros(ndl,nfc);
        ii = ii+1;
        N_1 = N_2+1;
        N_2 = setN_2(N_1,T_int,siglen);
    end

end


% ------ Subfunctions ----------------------------------------------------
function N_2 = setN_2(N_1,T_int,siglen)
% SETN_2 Sets the summation boundary N_2 to a new value
%
if T_int==Inf
    N_2 = siglen;
else
    N_2 = N_1 + T_int;
end



