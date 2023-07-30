function [b,a] = middleearfilter(fs,varargin)
%MIDDLEEARFILTER   Middle ear filter
%   Usage: [b,a]=middleearfilter(fs,varargin);
%          [b,a]=middleearfilter(fs);
%          [b,a]=middleearfilter;
%
%   MIDDLEEARFILTER(fs) computes the filter coefficients of a FIR (or IIR)
%   filter approximating the effect of the middle ear.
%
%   The following parameter and flags can be specified additionally:
%
%     'order',order  Sets the filter order of the computed FIR filter.
%                    Default value is 512.
%
%     'minimum'      Calculates a minimum phase filter. This is the default.
%
%     'zero'         returns a filter with zero phase. Since Matlab shifts the
%                    symmetric impulse response due to no negative indices.
%                    This results in a linear phase and hence a delay in the 
%                    signal chain.
%
%     'lopezpoveda2001'  Use data from Lopez-Poveda and Meddis (2001). These
%                    data are in turn derived from Goode et al. (1994).
%                    This is the default. 
%
%     'jepsen2008'  Use the data originally used in the Jepsen et al. (2008).
%
%     'verhulst2012' IIR filter approximating the middle ear transfer function
%                    based on Puria2003 (M1 filter) as used by Verhulst 2012.
%
%     'verhulst2015' IIR filter approximating the middle ear transfer function
%                    based on Puria2003 (M1 filter) as used by Verhulst 2015.
%
%     'verhulst2018' IIR filter approximating the middle ear transfer function
%                    based on Puria2003 (M1 filter) as used by Verhulst 2018.
%
%     'zilany2009'   Second-order cascade IIR filters approximating the 
%                    middle ear transfer function described by Ibrahim 
%                    (2012, Appendix).
%
%     'zilany2009cat' Second-order cascade IIR filters approximating the 
%                    middle ear transfer function described by Zilany et al.
%                    (2006, their Eq. 1--3)
%
%   MIDDLEEARFILTER without any input arguments returns a table describing
%   the frequency response of the middle ear filter. First column of the
%   table contain frequencies and the second column contains the amplitude
%   (stapes peak velocity in m/s at 0dB SPL) of the frequency like in figure
%   2b) of Lopez-Poveda and Meddis (2001).
%
%   MIDDLEEARFILTER is meant to be used in conjunction with the LOPEZPOVEDA2001
%   function, as the output is scaled to make lopezpoveda2001 work. If you are not
%   using the lopezpoveda2001, you probably do not want to call this function. The
%   following code displays the magnitude response of the filter:
%
%     fs=16000;
%     x=erbspace(0,fs/2,100);
%     b=middleearfilter(fs);
%     H=freqz(b,1,x,fs);
%     semiaudplot(x,10*log10(abs(H).^2));
%     xlabel('Frequency (Hz)');
%     ylabel('Magnitude (dB)');
%
%   See also:  data_lopezpoveda2001, lopezpoveda2001
% 
%   References:
%     R. Ibrahim. The role of temporal fine structure cues in speech
%     perception. Ph.d., McMaster University, 2012.
%     
%     R. Goode, M. Killion, K. Nakamura, and S. Nishihara. New knowledge
%     about the function of the human middle ear: development of an improved
%     analog model. The American journal of otology, 15(2):145--154, 1994.
%     
%     E. Lopez-Poveda and R. Meddis. A human nonlinear cochlear filterbank.
%     J. Acoust. Soc. Am., 110:3107--3118, 2001.
%     
%     M. Zilany and I. Bruce. Modeling auditory-nerve responses for high
%     sound pressure levels in the normal and impaired auditory periphery. J.
%     Acoust. Soc. Am., 120:1446--1466, 2006.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/middleearfilter.php


%   #Author: Peter L. Søndergaard 
%   #Author: Katharina Egger 
%   #Author: Alejandro Osses (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Check input options --------------------------------------------

a = 1; % if a == 1 then it is a FIR filter (not the case for the Verhulst et
       %    et al middle-ear filters, based on Puria2003)

% Define input flags
definput.flags.filtertype = {'lopezpoveda2001','jepsen2008', ...
                              'verhulst2012','verhulst2015','verhulst2018', ...
                              'zilany2009','zilany2009cat'};
definput.flags.phase = {'minimum','zero'};
definput.keyvals.order = 512;

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_lopezpoveda2001
  
  data = data_lopezpoveda2001('fig2b', 'no_plot');
  
  if nargin==0
    b = data;
  else
    if fs<=20000
      % In this case, we need to cut the table because the sampling
      % frequency is too low to accomodate the full range.
      indx=find(data(:,1)<fs/2);
      data = data(1:indx(end),:);
    else
      % otherwise the table will be extrapolated towards fs/2
      % data point added every 1000Hz
      lgth = size(data,1);
      for ii = 1:floor((fs/2-data(end,1))/1000)
        data(lgth+ii,1) = data(lgth+ii-1,1) + 1000;
        % 1.1 corresponds to the decay of the last amplitude values = approx. ratio
        % between amplitudes of frequency values seperated by 1000Hz
        data(lgth+ii,2) = data(lgth+ii-1,2) / 1.1; 
      end
    end;  
    
    % for the function fir2 the last data point has to be at fs/2
    lgth = size(data,1);
    if data(lgth,1) ~= fs/2
      data(lgth+1,1) = fs/2;
      data(lgth+1,2) = data(lgth,2) / (1+(fs/2-data(lgth,1))*0.1/1000);
    end
    
    % Extract the frequencies and amplitudes, and put them in the format
    % that fir2 likes.
    freq=[0;...
          data(:,1).*(2/fs);...
         ];
    ampl=[0;...
          data(:,2);...
         ];
    
    b = fir2(kv.order,freq,ampl);
    
    b = b / 20e-6;      % scaling for SPL in dB re 20uPa
    
    if flags.do_minimum
      X = fft(b);
      Xmin = abs(X) .* exp(-1j*imag(hilbert(log(abs(X)))));
      b = real(ifft(Xmin));        
    end
    
  end
  
end;

if flags.do_jepsen2008
    
  stapes_data = [...
      50,	 48046.39731;...
      100, 24023.19865;...
      200, 12011.59933;...
      400,  6005.799663;...
      600, 3720.406871;...
      800,  2866.404385;...
      1000, 3363.247811;...
      1200, 4379.228921;...
      1400, 4804.639731;...
      1600, 5732.808769;...
      1800, 6228.236688;...
      2000, 7206.959596;...
      2200, 9172.494031;...
      2400, 9554.681282;...
      2600, 10779.64042;...
      2800, 12011.59933;...
      3000, 14013.53255;...
      3500, 16015.46577;...
      4000, 18017.39899;...
      4500, 23852.82136;...
      5000, 21020.29882;...
      5500, 22931.23508;...
      6000, 28027.06509;...
      6500, 28745.70779;...
      7000, 32098.9;...
      7500, 34504.4;...
      8000, 36909.9;...
      8500, 39315.4;...
      9000, 41720.9;...
      9500, 44126.4;...
      10000,46531.9;...
                ];
  
  % We need to find inverse because original data is stapes impedance and we
  % need stapes velocity.
  stapes_data (:,2) = 1./stapes_data(:,2); 
  
  if nargin==0
    b = stapes_data;
  else
    
    if fs<=20000
      % In this case, we need to cut the table because the sampling
      % frequency is too low to accomodate the full range.
      
      indx=find(stapes_data(:,1)<fs/2);
      stapes_data=stapes_data(1:indx(end),:);
    end;  
    
    % Extract the frequencies and amplitudes, and put them in the format
    % that fir2 likes.
    freq=[0;...
          stapes_data(:,1).*(2/fs);...
          1];
    ampl=[0;...
          stapes_data(:,2);...
          0];
    
    b = fir2(kv.order,freq,ampl);
    
    % See the figure text for figure 1, Lopez (2001).
    b = b/max(abs(fft(b)))*1e-8*10^(104/20); 
    
  end;      
  
end;
  
if flags.do_verhulst2012 || flags.do_verhulst2015 || flags.do_verhulst2018
    
    %%% PuriaM1 filter (Puria2003, Fig. 2A and 3A)
    if flags.do_verhulst2012
        % Source: AMT implementation and originally in fortran code
        fc1 =  100; % as in cochlear_model.py, ../model2012-fortran/SourceCode/PuriaM1.f90
        fc2 = 3000; % as in cochlear_model.py, ../model2012-fortran/SourceCode/PuriaM1.f90
        puria_gain = 2*gaindb(1,18); % puria_gain=10**(18./20.)*2. % The factor of 2 is explained 
                    % as: 'VoltageDivisionGain=2d0 !2 resistor matching network(Zme+Zoch)'
    end
    if flags.do_verhulst2015
        % Source: ../model2015/cochlear_model.py
        fc1 =  600; % Hz
        fc2 = 3000; % Hz
        puria_gain = 2*gaindb(1,18); % same gain as in verhulst2012
    end
    if flags.do_verhulst2018
        fc1 = 600; % as in cochlear_model2018.py
        fc2 = 4000; % as in cochlear_model2018.py
        puria_gain = gaindb(1,18); % puria_gain=10**(18./20.)
    end
    [b,a] = butter(1,[fc1 fc2]/(fs/2),'bandpass');
    b = b*puria_gain;
    
end  
  
if flags.do_zilany2009
    % The following humanised middle ear filter is a digital implementation
    % described by Rasha Ibrahim's thesis from 2012. She based this implementation
    % on the work described by Pascal et al. (JASA 1998)  */
    %
    % There are three resulting (2nd-order) filters. For an fs=100 kHz, the
    % transfer functions to be obtained are (see Ibrahim, 2012, her Eq. A1, A2, and A3):
    %
    %            0.9979 - 1.9408 z^-1 + 0.9429 z^-2
    % H1(z-1) = ------------------------------------
    %            1.0000 - 1.9395 z^-1 + 0.9420 z^-2
    %
    %            0.9984 - 1.9226 z^-1 + 0.9415 z^-2
    % H2(z-1) = ------------------------------------
    %            1.0000 - 1.9244 z^-1 + 0.9379 z^-2
    %
    %                  0.0286 + 0.0302 z^-1 + 0.0016 z^-2
    % H3(z-1) = 0.5 x ------------------------------------
    %                  1.0000 - 1.6748 z^-1 + 0.7847 z^-2
    % Note that the 0.5 is not defined in Rasha's thesis, but in the code 
    % that was included in the AMT toolbox version of zilany2009's model.
    
    fp = 1e3; % prewarping frequency 1 kHz
    C  = 2*pi*fp/tan(pi*fp/fs);
    
    
    m11=1/(C^2+5.9761e3*C+2.5255e7);
    m12=(-2*(C^2)+2*2.5255e7);
    m13=(C^2-5.9761e3*C+2.5255e7);
    m14=(C^2+5.6665e3*C);
    m15=-2*(C^2);
    m16=(C^2-5.6665e3*C);
    m21=1/(C^2+6.4255e3*C+1.3975e8);
    m22=(-2*(C^2)+2*1.3975e8);
    m23=(C^2-6.4255e3*C+1.3975e8);
    m24=(C^2+5.8934e3*C+1.7926e8); 
    m25=-2*(C^2)+2*1.7926e8;	
    m26=C^2-5.8934e3*C+1.7926e8;
    m31=1/(C^2+2.4891e4*C+1.2700e9);
    m32=(-2*(C^2)+2*1.27e9);
    m33=(C^2-2.4891e4*C+1.27e9);
    m34=(3.1137e3*C+6.9768e8);
    m35=2*6.9768e8;
    m36=(-3.1137e3*C+6.9768e8);
    megainmax=2;
    a(1,1:3) =     [1   m11*m12 m11*m13];
    b(1,1:3) = m11*[m14     m15     m16];

    a(2,1:3) =     [1   m21*m22 m21*m23];
    b(2,1:3) = m21*[m24     m25     m26]; 

    a(3,1:3) = [1 m31*m32 m31*m33];
    b(3,1:3) = m31*[m34 m35 m36]/megainmax;    
end

if flags.do_zilany2009cat
    % The following middle ear filter of the cat is a digital implementation
    % described by Rasha Ibrahim's thesis from 2012. She based this implementation
    % on the work described by Pascal et al. (JASA 1998)  */
    %
    % There are three resulting (2nd-order) filters. The obtained 2nd order
    % transfer functions at an fs=500 kHz can be found in (Zilany et al, 2006, 
    % their Eq. 1--3):
    
    fp = 1e3; % prewarping frequency 1 kHz
    C  = 2*pi*fp/tan(pi*fp/fs);
    
    m11 = C/(C + 693.48);
    m12 = (693.48-C)/C;
    m13 = 0.0;
    m14 = 1.0;
    m15 = -1.0;
    m16 = 0.0;
    m21 = 1/(C^2 + 11053*C + 1.163e8);  
    m22 = -2*(C^2) + 2.326e8;     
    m23 = C^2 - 11053*C + 1.163e8; 
    m24 = C^2 + 1356.3*C + 7.4417e8;    
    m25 = -2*(C^2) + 14.8834e8;   
    m26 = C^2 - 1356.3*C + 7.4417e8;
    m31 = 1/(C^2 + 4620*C + 909059944); 
    m32 = -2*(C^2) + 2*909059944; 
    m33 = C^2 - 4620*C + 909059944;
    m34 = 5.7585e5*C + 7.1665e7;
    m35 = 14.333e7;
    m36 = 7.1665e7 - 5.7585e5*C;
    megainmax=41.1405;
    
    a(1,1:3) =     [1   m11*m12 m11*m13];
    b(1,1:3) = m11*[m14     m15     m16];

    a(2,1:3) =     [1   m21*m22 m21*m23];
    b(2,1:3) = m21*[m24     m25     m26]; 

    a(3,1:3) = [1 m31*m32 m31*m33];
    b(3,1:3) = m31*[m34 m35 m36]/megainmax;
end

if nargout == 0
    
    if flags.do_zilany2009 || flags.do_zilany2009cat
        
        N = round(fs/2);
        
        insig = [1; zeros(N-1,1)];
        Nr_cascaded = size(b,1);
        for k = 1:Nr_cascaded
            insig = filter(b(k,:),a(k,:),insig);
        end
        
        [h,w] = freqz(insig,1,N);
        f = (w/pi)*fs/2;
        figure; 
        semilogx(f,20*log10(abs(h)));
        xlim([20 fs/2]);
        
        ylim([-50 10])
        grid on
        
        xlabel('Frequency (Hz)');
        xlabel('IIR middleearfilter');
        
    end
    
end

% if flags.do_plot
%     % Manually calculate the frequency response
%     fmid = abs(fftreal(b));
%     % Half the filter length.
%     n2=length(fmid);
%     % x-values for plotting.
%     xplot=linspace(0,fs/2,n2);
%     loglog(xplot/1000,fmid);
%     xlabel('Frequency (kHz)');
%     ylabel('FIR middleearfilter');
% end



