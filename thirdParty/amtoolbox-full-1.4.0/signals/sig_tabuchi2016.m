function [WaveMout, nPhase] = sig_tabuchi2016(CondName, sfs, f0, n1, n2, spl, durms, C, varargin)
% SIG_TABUCHI2016 Generate schroeder-phase harmonic complex with Roex frequency weighting
%
%   Input parameters:
%     CondName     : can be OnFreqPre or OffFreqPre
%     sfs          : sampling frequency
%     f0           : fundamental frequency
%     n1           : the lowest frequency component (f0*n1 kHz)
%     n2           : the highest frequency component (f0*n2 kHz)
%     spl          : level [dB SPL]
%     durms        : duration [ms]
%     C            : C (curvature) value
%
%   Output parameters:
%     WaveMout     : output signal
%     nPhase       : phase index
%
%   this function generates a Roex-weighted Schroeder complex
%
%   Optional parameters:
%
%     'leadms',leadms      leading silence [ms]
%
%     'trailms',trailms    trailing silence [ms]
%
%     'pref',pref          reference pressure [pa]
%
%   See also: tabuchi2016
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_tabuchi2016.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Verified
%   #Requirements: 
%   #Author: Hisaaki Tabuchi
%   #Author: Clara Hollomey (adaptations for AMT)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.

ListenerStr = 'GrandMeanSd';

if strcmp(CondName,'OffFreqPre')
    CondStr = 'OffFreq';
elseif strcmp(CondName,'OnFreqPre')
    CondStr = 'OnFreq';
else
    error('CondName is not defined.')
end

if spl == 90 && strcmp(CondName,'OnFreqPre')
    f_wei_lin_vec = [];
else
    data = amt_load('tabuchi2016', ['Attenu_Roex_' ListenerStr '_', ...
    CondStr '_' int2str(spl), 'dB.mat'], 'f_wei_lin_vec');
    f_wei_lin_vec = data.f_wei_lin_vec;
end

definput.keyvals.leadms = 0;
definput.keyvals.trailms = 0;
definput.keyvals.pref=0.00002;

[~,kv]  = ltfatarghelper({},definput,varargin);

% synthesize
smpl = 1/sfs;
tt = (smpl:smpl:durms)'; % waveform should be a column vector %

% Schroeder harmonic complexes: Lentz et al (2001)
nPhase = NaN(n2-n1+1,1);
wav1 = zeros(length(tt),1);
MatFlag = 1;
for harmnumLoop = n1:n2 
    nPhase(MatFlag,1) = C * pi * harmnumLoop * (harmnumLoop+1)/(n2-n1+1);
    if isempty(f_wei_lin_vec) == 0 % frequency weighting by Roex filter
        wav1 = wav1 + f_wei_lin_vec(MatFlag) * cos(( 2 * pi * harmnumLoop * f0 * tt)+ nPhase(MatFlag,1) );
    else % no frequency weighting
        wav1 = wav1 + cos(( 2 * pi * harmnumLoop * f0 * tt)+ nPhase(MatFlag,1) );        
    end
    MatFlag = MatFlag + 1;
end


% leading zeros
nz1 = round(kv.leadms/smpl);
wav1 = cat(1,zeros(nz1,1),wav1);

% trailing zeros
nz2 = round(kv.trailms/smpl);
wav1 = cat(1,wav1,zeros(nz2,1));

% scale amplitude 
% check amplitude, excluding leading/trailing silence
mxx=max(abs(wav1));
n1=find( wav1>(0.02*mxx),1,'first' );
n2=find( wav1>(0.02*mxx),1,'last' );

ss=sum(wav1(n1:n2).*wav1(n1:n2));
rrr = sqrt( ss/length(wav1(n1:n2)) );
if rrr>0
    dbtmp=20*log10( rrr/kv.pref );
    dbdif=spl-dbtmp;
    sfact=10^(dbdif/20);
    WaveMout=sfact*wav1;
else
    WaveMout=wav1;
end


end

