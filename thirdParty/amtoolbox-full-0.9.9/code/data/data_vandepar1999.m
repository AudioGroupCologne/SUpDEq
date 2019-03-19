function data = data_vandepar1999(varargin)
%DATA_VANDEPAR1999 Returns data points from the van de Par and Kohlrausch (1999) paper
%   Usage: data = data_vandepar1999(figure,nfc)
%
%   DATA_VANDEPAR1999(figure,nfc) returns data from van de Par and
%   Kohlrausch (1999)
%
%   The flags for the figure may be one of:
%
%     'fig1_N0S0'     Returns the N0S0 values for figure 1
% 
%     'fig1_N0Spi'     Returns the N0Spi values for figure 1
%
%     'fig1_NpiSo'    Returns the NpiS0 values for figure 1
%
%   The nfc (center frequency of noise) may be one of: 
%     
%     'nfc125'
%     
%     'nfc250'
%     
%     'nfc500'
%     
%     'nfc1000'
%     
%     'nfc2000'
%     
%     'nfc4000'
%
%   Examples:
%   ---------
% 
%   To get data for the fig. 1 van de Par and Kohlrausch (1999) for the 
%   N0S0 condition with 125 Hz center frequency use :
%
%     data_vandepar1999('fig1_N0S0','nfc125');
%
%   References:
%     S. van de Par and A. Kohlrausch. Dependence of binaural masking level
%     differences on center frequency, masker bandwidth, and interaural
%     parameters. J. Acoust. Soc. Am., 106(4):1940-1947, 1999.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_vandepar1999.php

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


%   AUTHOR: Martina Kreuzbichler


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type={'missingflag','fig1_N0S0','fig1_N0Spi','fig1_NpiS0'};
definput.flags.nfc = {'missingnfc','nfc125','nfc250','nfc500','nfc1000','nfc2000','nfc4000'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

if flags.do_missingnfc
  nfcnames=[sprintf('%s, ',definput.flags.nfc{2:end-2}),...
             sprintf('%s or %s',definput.flags.nfc{end-1},definput.flags.nfc{end})];
  error('%s: You must specify one of the following center frequencies: %s.',upper(mfilename),nfcnames);
end;

%% ------ Data points from the paper ------------------------------------
%
% Data for the given figure
if flags.do_fig1_N0S0
    if flags.do_nfc125
        data= [3 2 1 0.8 -0.5 -6.5];
    elseif flags.do_nfc250
        data = [2.5 2.5 -1 0 -3.5 -8.5 -12.5];
    elseif flags.do_nfc500
        data = [2.5 1.5 -0.5 -2 -3 -7 -11 -11];
    elseif flags.do_nfc1000
        data = [2.5 1.5 -1 -2 -3.5 -6 -8.5 -11.5 -15.5];
    elseif flags.do_nfc2000
        data = [3 2 0 -2 -3.5 -4 -4.5 -9 -12.5 -17];
    elseif flags.do_nfc4000
        data = [2 3 0 -1.5 -2.5 -5.5 -4.5 -5.5 -8.5 -13];
    end
    
elseif flags.do_fig1_N0Spi
    if flags.do_nfc125
        data= [-20 -19.5 -18.5 -17 -16.7 -19.7];
    elseif flags.do_nfc250
        data = [-22.5 -23 -22.5 -23.5 -23 -25 -27.5];
    elseif flags.do_nfc500
        data = [-21 -22 -23.5 -23.3 -23 -24 -25.5 -29.5];
    elseif flags.do_nfc1000
        data = [-18.5 -19.5 -19.5 -20 -21 -19 -18 -22 -22.5];
    elseif flags.do_nfc2000
        data = [-13.5  -14.5 -13.7 -16 -16.5 -13.3 -13.7 -13.5 -17 -21];
    elseif flags.do_nfc4000
        data = [-9.5 -8.5 -9.5 -8.5 -11.7 -12.5 -11.5 -10.5 -10 -15];
    end
    
elseif flags.do_fig1_NpiS0
    if flags.do_nfc125
        data=  [-14 -14.5 -12.5 -8.5 -7 -12];
    elseif flags.do_nfc250
        data = [-17 -18.5 -18 -16 -14.5 -19 -21];
    elseif flags.do_nfc500
        data = [-17.5 -20.5 -18.5 -18.5 -18 -19 -22.5 -26];
    elseif flags.do_nfc1000
        data = [-16 -17.5 -16.5 -18 -18 -17 -16 -19.5 -22.5];
    elseif flags.do_nfc2000
        amt_disp('Not available');
        data = [];
    elseif flags.do_nfc4000
        amt_disp('Not available');
        data = [];
    end
        
end
