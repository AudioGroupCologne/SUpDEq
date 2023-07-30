function stim = zakarauskas1993_inputfilter( in,flag )
%ZAKARAUSKAS1993_INPUTFILTER filters input signal with source spectra given in Fig.7 of Zakarauskas et al.(1993)
%
%   Usage:
%     [ stim ] = zakarauskas1993_inputfilter( in,flag )
%
%   Input parameters:
%     in:       input signal, or set to 'imp' if it should be an impulse
%     flag:     'paper' for using the spectrum of crumbling a paper (Fig.7 dotted)
%               'shock' for using the spectrum of word "shock" (Fig.7 solid)
%
%   Filters the input signal with source spectra given in Fig.7 of Zakarauskas et al.(1993)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/zakarauskas1993_inputfilter.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Robert Baumgartner (2010), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if strcmp('imp',in)==1
    in=zeros(11000,1);in(1)=1;
end

fs=22000;
len=100;
fi=linspace(100,fs/2,1024)/(fs/2);

switch flag
    case 'paper'
        fp  = [0.0   0.080   0.0909    0.1045    0.1364    0.2000    0.2909    0.4091    0.8182    1.0];
        apdB= [0     0       18         19.5        19      18       15        16.8      3.5       0.0];
        ap=10.^(apdB/20);
        hp=interp1(fp,ap,fi);
        bp = firls(len,fi,hp);
        stim = filter(bp,1,in);
    case 'shock'
        fw  = [0.0  0.05     0.080   0.0909    0.0936    0.1018    0.1455    0.2364    0.2818    0.3000    0.5364   0.7273    0.9091  1.0];
        awdB= [0     0       0       13.0       12.0       13.0     7.0        16       12       14        11         1         0.5   0.0];
        aw=10.^(awdB/20);
        hw=interp1(fw,aw,fi);
        bw = firls(len,fi,hw);
        stim = filter(bw,1,in);
end
end


