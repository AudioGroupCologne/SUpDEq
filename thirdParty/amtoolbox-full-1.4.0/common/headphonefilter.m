function b = headphonefilter(fs,order)
%HEADPHONEFILTER  Combined headphone and outer ear filter
%   Usage:  b=headphonefilter(fs,order);
%           b=headphonefilter(fs);
%           b=headphonefilter;
%
%   HEADPHONEFILTER(fs,order) computes the filter coefficients of a FIR
%   filter or order order approximating the combined effect of headphones
%   and the outer ear. The data describes a generic set of headphones,
%   originally from Pralong et al. (1996)
%
%   HEADPHONEFILTER(fs) does the same with a FIR filter of order 512.
%
%   HEADPHONEFILTER without any input arguments returns a table
%   describing the frequency response of the headphone filter. First
%   column of the table contain frequencies and the second column
%   contains the amplitude of the frequency.
%
%   The following code displays the magnitude response of the filter:
%
%     fs=16000;
%     x=erbspace(0,fs/2,100);
%     b=headphonefilter(fs);
%     H=freqz(b,1,x,fs);
%     semiaudplot(x,10*log10(abs(H).^2));
%     xlabel('Frequency (Hz)');
%     ylabel('Magnitude (dB)');
%
%   See also: middleearfilter, data_pralong1996, data_lopezpoveda2001
%
%   References:
%     E. Lopez-Poveda and R. Meddis. A human nonlinear cochlear filterbank.
%     J. Acoust. Soc. Am., 110:3107--3118, 2001.
%     
%     D. Pralong and S. Carlile. The role of individualized headphone
%     calibration for the generation of high fidelity virtual auditory space.
%     J. Acoust. Soc. Am., 100:3785--3793, 1996.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/headphonefilter.php


%   #Author: Morten L. Jepsen: Original version
%   #Author: Peter L. Søndergaard (2009-2013): adapted to AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



if nargin==1
  order = 512;    % desired FIR filter order
end;

eardrum_data = data_pralong1996;

if nargin==0
  b = eardrum_data;
else
  
  if fs<=20000
    % In this case, we need to cut the table because the sampling
    % frequency is too low to accomodate the full range.
    
    indx=find(eardrum_data(:,1)<fs/2);
    eardrum_data=eardrum_data(1:indx(end),:);
  end;  

  % Extract the frequencies and amplitudes, and put them in the format
  % that fir2 likes.
  freq=[0;...
        eardrum_data(:,1).*(2/fs);...
        1];
  ampl=[0;...
        eardrum_data(:,2);...
        0];
  
  b = fir2(order,freq,ampl);

end;



