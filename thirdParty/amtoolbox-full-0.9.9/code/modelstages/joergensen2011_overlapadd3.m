function [ReconstructedSignal xx]=OverlapAdd3(XNEW,yphase,windowLen,ShiftLen)

%Y=OverlapAdd(X,A,W,S);
%Y is the signal reconstructed signal from its spectrogram. X is a matrix
%with each column being the fft of a segment of signal. A is the phase
%angle of the spectrum which should have the same dimension as X. if it is
%not given the phase angle of X is used which in the case of real values is
%zero (assuming that its the magnitude). W is the window length of time
%domain segments if not given the length is assumed to be twice as long as
%fft window length. S is the shift length of the segmentation process ( for
%example in the case of non overlapping signals it is equal to W and in the
%case of %50 overlap is equal to W/2. if not givven W/2 is used. Y is the
%reconstructed time domain signal.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/joergensen2011_overlapadd3.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

wnd=hanning(windowLen,'periodic');

[FreqRes FrameNum]=size(XNEW);

Spec=XNEW.*exp(1j*yphase);


 if mod(windowLen,2) %if FreqResol is odd
     Spec=[Spec;flipud(conj(Spec(2:end,:)))];
 else
     Spec=[Spec;flipud(conj(Spec(2:end-1,:)))];
 end

% envelope=fades(length(wnd),512,2);


sig=zeros((FrameNum-1)*ShiftLen+windowLen,1);


for i=1:FrameNum
    start=(i-1)*ShiftLen+1;
    spec=Spec(:,i);
    xx=real(ifft(spec,windowLen));
%     xx=xx.*envelope;
    xxx(:,i)=xx;
    sig(start:start+windowLen-1)=sig(start:start+windowLen-1)+xx;
end


ReconstructedSignal=sig;


