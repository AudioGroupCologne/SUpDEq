%DEMO_ZAKARAUSKAS1993 script for localization model according to ZAKARAUSKAS et al.(1993)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_zakarauskas1993.php


%   #Author: Robert Baumgartner (2010)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

b = amt_load('zakarauskas1993', 'hrtf_M_KEMAR-gardner 22kHz.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global settings
fs=b.stimPar.SamplingRate;       % sampling frequency
plotflag=1;     % switch for plots
lat=00;         % lateral angle of sagital plane
delta=0;        % lateral variability in degree
% model settings
q10=10;         % relativ bandwidth of filter bands;  default: 10
space=1.10;     % spacing factor of filter bank (distance of bands)
bal= 1;         % balance of left to right channel;   default: 1
fstart =1000;   % start frequency of filter bank;     default: 1kHz
fend   =11000;  % end frequency of filter bank;       default: 11kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx=find(b.posLIN(:,4)>=-(delta+0.01)/2+lat & b.posLIN(:,4)<=(delta+0.01)/2+lat); % median plane with lateral delta-variation; +0.01 because of rounding errors in posLIN
[pol,polidx]=unique(b.posLIN(idx,5));   % sorted polar angles
sagir=double(b.hM(:,idx,:));            % unsorted impulse responses on sagital plane
ir=sagir(:,polidx,:);                   % sorted

paper = zakarauskas1993_inputfilter( 'imp','paper' );
shock = zakarauskas1993_inputfilter( 'imp','shock' );

h = waitbar(0,'Please wait...');

do = 1;
[dp,hitsdp,errdp,aedp] = zakarauskas1993( ir,pol,paper,q10,do,bal,fstart,fend,fs,space ); 
waitbar(1/4)
[ds,hitsds,errds,aeds] = zakarauskas1993( ir,pol,shock,q10,do,bal,fstart,fend,fs,space ); 
waitbar(2/4)

do = 2;
[cp,hitscp,errcp,aecp] = zakarauskas1993( ir,pol,paper,q10,do,bal,fstart,fend,fs,space );
waitbar(3/4)
[cs,hitscs,errcs,aecs] = zakarauskas1993( ir,pol,shock,q10,do,bal,fstart,fend,fs,space );
waitbar(4/4)
close(h)

% plots
if plotflag==1
    plot_zakarauskas1993_dif( ds,cs,pol,find(pol>=190,1,'first'),'word "shock"' );
    plot_zakarauskas1993_responsepattern( errds,errcs,pol,'word "shock"' );
    % result table
    disp('_______________________________________');
    disp('               Operator  Paper   Shock');
    disp(sprintf('#exact/%d:         D      %d      %d',size(cp,1),hitsdp,hitsds));
    disp(sprintf('#exact/%d:         C      %d      %d',size(cp,1),hitscp,hitscs));
    disp(sprintf('Average error:     D    %5.2f�  %5.2f�',aedp,aeds));
    disp(sprintf('Average error:     C    %5.2f�  %5.2f�',aecp,aecs));
    disp('_______________________________________');
end

