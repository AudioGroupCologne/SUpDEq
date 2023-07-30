function [output_Y, output_N] = joergensen2011_sepsub(input,noise,W,padz,SP,factor)
%JOERGENSEN2011_SEPSUB estimates the clean signal output_Y
%
%   Usage: [output Nzeros] = joergensen2011_sepsub(input,noise,W,padz,SP,factor,fs)
%
%   Input parameters:
%     input     : Vector containg the noisy input (time) signal.
%     noise     : Vector containg the noise (time) signal that was added to the clean signal to create the noisy signal.
%     W         : Frame length  
%     padz      : zero padding (pad with padz/2 from the left and padz/2 from the right )
%     SP        : Shift percentage (overlap)
%     factor    : the over-subtraction factor
% 
%   Output parameters:
%     output    : the estimated "clean" time signal
%
%   This function calculates an estimate of the clean signal Y_hat from the
%   noisy signal (signal) and noise alone using spectral subtraction as
%   defined by Berouti et al. (1979).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/joergensen2011_sepsub.php


%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Untrusted
%   #Requirements: M-Signal M-Stats
%   #Author: Søren Jørgensen (2011)
%   #Author: Peter L. Sondergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% Typical values:

% fs=44100; %sampling frequency
% W=1024; % frame length
% padz=1024; %zero padding (pad with padz/2 from the left and padz/2 from the right )
% % Note that (W+padz) is the final frame window and hence the fft length (it is normally chose as a power of 2)
% SP=0.5; %Shift percentage is 50%

if nargin < 6
    error('Insufficient amount of arguments');
    return;
end

stim = zeros(2,length(input));
stim(1,:) = input;
stim(2,:) = noise;

%  CALCULATE BASIC VALUES

wnd=(hanning(W)); %create hanning window with length = W
% fr=ceil((2*RT*fs-2*W)/W+1); %fr is the number of past frames that the method will take into account in order to extract the reverberation spectrum estimation

for k = 1:2
    % CUT THE THE APPROPRIATE SIGNAL FRAMES
    
    wnd1=wnd(:); %make it a column vector
    L=length(stim(k,:));
    SP1=fix(W.*SP);
    N1=fix((L-W)/SP1 +1); %number of segments
    Index=(repmat(1:W,N1,1)+repmat((0:(N1-1))'*SP1,1,W))';
    hw=repmat(wnd1,1,N1);
    y_tmp =stim(k,:);
    y_tmp=y_tmp(Index).*hw;
    
    % PAD WITH ZEROS
    
    pad=zeros((padz/2),length(y_tmp(10,:)));
    y_pad=[pad' y_tmp' pad'];
    y_tmp=y_pad';
    y(k,:,:) = y_tmp;
    clear y_pad y_tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FREQUENCY DOMAIN

% signal:
y_tmp = reshape(y(1,:,:),size(y,2),size(y,3));
Y=(fft(y_tmp));
YY=(Y(1:round(end/2)+1,:)); % Half window (exploite the symmmetry)
YPhase=angle(Y(1:round(end/2)+1,:)); % Phase
Y1=abs(Y(1:round(end/2)+1,:)); % Spectrum
Y2=abs(Y(1:round(end/2)+1,:)).^2;% Power Spectrum
numberOfFrames=size(Y,2);
FreqResol=size(Y,1);

%  noise:
y_tmp = reshape(y(2,:,:),size(y,2),size(y,3));
Y_N=(fft(y_tmp));
YY_N=(Y_N(1:round(end/2)+1,:)); % Half window (exploite the symmmetry)
Y_NPhase=angle(Y_N(1:round(end/2)+1,:)); % Phase
Y_N1=abs(Y_N(1:round(end/2)+1,:)); % Spectrum
Y_N2=abs(Y_N(1:round(end/2)+1,:)).^2;% Power Spectrum

for i=3:numberOfFrames
    
%   The noise "estimate" is simply the average of the noise power spectral
%   density in the frame:
    P_N=mean(Y_N2((1:length(Y_N2(:,1))/1),i));   % Power of the noise frames is mean across frequency
%   Subtraction of the noise estimate from the noisy signal: 
    Y_hat(:,i)=Y2(:,i) - factor*P_N; % subtraction
    PN_hat(:,i)=Y_N2(:,i) - factor*P_N; % subtraction for noise alone 
%   Truncating the result to be <= 0:
    Y_hat(:,i)=max(Y_hat(:,i),0); % Makes the minimum equal zeros
    zero_INDEXES = Y_hat(:,i) == 0;
    PN_hat(zero_INDEXES,i) = 0;
   %  
end
% Y_hat = Y_hat(:,3:end);
 Nzeros = find(Y_hat(:,3:end) == 0);
 if isempty(Nzeros), Nzeros = 0; end
 Nzeros = length(Nzeros);    
% Combining the estimated power spectrum with the original noisy phase, and
% add the frames using an overlap-add technique
% 
output_Y=joergensen2011_overlapadd3(Y_hat.^(1/2),YPhase,(W+padz),SP*W);
output_N=joergensen2011_overlapadd3(PN_hat.^(1/2),Y_NPhase,(W+padz),SP*W);


