function data = data_breebaart2001(varargin)
%DATA_BREEBAART2001 Returns data points from the Breebaart et al. (2001b) paper
%   Usage: data = data_breebaart(figure,nfc)
%
%   data_breebaart(figure,nfc) returns data from Breebaart et al. (2001b)
%
%   The flags fo the figure may be one of:
%
%     'fig3'      Returns the N0Spi values for figure 3
% 
%     'fig6'      Returns the NpiS0 values for figure 6
%
%   The nfc (center frequency of noise) may be one of: 
%
%     'nfc125','nfc250','nfc500','nfc1000': figure 3 and figure 6
%
%     'nfc2000','nfc4000':                  figure 3 only
%
%   For the data points one can choose between
%
%     'breebaartmodel'  Returns the data points derived from the model of
%                       Breebaart et al. (2001b). This is the default.
%                       
%     'vandepaar'       Returns the data points Breebaart used in Fig. 1
%                       from van de Paar and Kohlrausch (1999)
%
%   Examples:
%   ---------
% 
%   To get data for the fig. 3 Breebaart et al. (2001b) for the 
%   condition with 125 Hz center frequency use :
%
%     data_breebaart2001('fig3','nfc125');
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074--1088, August 2001.
%     
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. II. Dependence on spectral
%     parameters. J. Acoust. Soc. Am., 110:1089--1104, August 2001.
%     
%     S. van de Par and A. Kohlrausch. Dependence of binaural masking level
%     differences on center frequency, masker bandwidth, and interaural
%     parameters. J. Acoust. Soc. Am., 106(4):1940--1947, 1999.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_breebaart2001.php


%   #Author: Martina Kreuzbichler (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type={'missingflag','fig3','fig6'};
definput.flags.nfc = {'missingnfc','nfc125','nfc250','nfc500','nfc1000','nfc2000','nfc4000'};
definput.flags.datapoints = {'breebaartmodel','vandepar'};

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
if flags.do_fig3
    if flags.do_nfc125
        if flags.do_breebaartmodel  
            data= [-20.5 -21 -21.5 -22 -21 -23.5];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_N0Spi','nfc125');
        end
    elseif flags.do_nfc250
        if flags.do_breebaartmodel
            data = [-22 -22 -20.5 -21 -21 -24 -30];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_N0Spi','nfc250');
        end    
    elseif flags.do_nfc500
        if flags.do_breebaartmodel 
            data = [-20 -21 -21 -20 -21 -22 -26 -28.5];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_N0Spi','nfc500');
        end
    elseif flags.do_nfc1000
        if flags.do_breebaartmodel 
            data = [-18 -17.5 -18 -17 -18 -18 -20 -22 -25.5];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_N0Spi','nfc1000');
        end
    elseif flags.do_nfc2000
        if flags.do_breebaartmodel 
            data = [-13.5 -10.5 -11 -14 -14 -15 -15 -17.5 -18 -23];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_N0Spi','nfc2000');
        end
    elseif flags.do_nfc4000
        if flags.do_breebaartmodel 
            data = [-9.5 -10 -11.5 -14 -15.5 -17 -16.5 -16 -17.5 -21];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_N0Spi','nfc4000');
        end
    end
    
elseif flags.do_fig6
    if flags.do_nfc125
        if flags.do_breebaartmodel 
            data= [-12.5 -12 -11.5 -8.5 -6 -8.5];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_NpiS0','nfc125');
        end
    elseif flags.do_nfc250
        if flags.do_breebaartmodel 
            data = [-17 -16.5 -18.5 -19 -16 -20 -24];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_NpiS0','nfc250');
        end
    elseif flags.do_nfc500
        if flags.do_breebaartmodel 
            data = [-16.5 -19 -19 -20 -17.5 -18.5 -21.5 -25.5];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_NpiS0','nfc500');
        end
    elseif flags.do_nfc1000
        if flags.do_breebaartmodel 
            data = [-15 -17.5 -14.5 -16 -16.5 -15.5 -17 -21.5 -23.5];
        elseif flags.do_vandepar
            data = data_vandepar1999('fig1_NpiS0','nfc1000');
        end
    end
    
end
    
    


