function [ y ] = fade( x,t_fade,Fs )
%FADE smoothen a signals on- and offset
%
%   Input parameters:
%     x     : input signal
%     t_fade: fading time
%     Fs: sampling frequency
%
%   Output parameters:
%     y   : output signal
%
%     takes a sequency and adds a rise and fall to its begining
%     and its end to minimize the pectral splatter.
%
%   See also: sig_whitenoiseburst
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/fade.php


%   #Author: The AMT Team (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%%  Create the 'fade in' and 'fade out'
t_max=(length(x)-1)/Fs;
t_in=0:1/Fs:t_fade;

f_fade=(0.01*25)/t_fade; % the best fading frequency for t_fade=0.01s turns to be 26 Hz. Find the appropriate f_fade for other t_fade values.

%%
fade_in=(sin(2*pi*f_fade.*t_in)).^2;
%t_out=t_max:-1/Fs:(t_max-t_fade);
t_out=t_fade:-1/Fs:0;
fade_out=(sin(2*pi*f_fade.*t_out)).^2;
%% multiply with the original signal
y=zeros(size(x));y=x;
y(1:length(t_in))=fade_in.*x(1:length(t_in));
y(length(y)-length(t_out)+1:length(y))=fade_out.*x(length(x)-length(t_out)+1:length(x));
%% display
% figure(2)
% plot(t_in,fade_in);
% figure(3)
% plot(length(y)-length(t_out)+1:length(y),fade_out);
% figure(4)
% plot(0:1/Fs:(length(y)-1)/Fs, y);

end


