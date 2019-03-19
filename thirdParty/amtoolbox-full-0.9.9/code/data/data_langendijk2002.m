function data = data_langendijk2002(varargin)
%DATA_LANGENDIJK2002  Data from Langendijk & Bronkhorst (2002)
%   Usage: data = data_langendijk2002(flag)
%
%   Output parameters:
%      data  : The data points from the given figure;
%
%   DATA_LANGENDIJK2002(flag) returns data points from the paper by 
%   Langendijk & Bronkhorst (2002). 
%
%   In the case of response patterns (Fig. 7 & 9) the first row of data*
%   describes target position and the second one belongs to the response
%   position.  In the case of DTF data (Fig. 11) the first dimension of the
%   data matrix describes frequency and the second one angle position, the
%   *first column* defines the actual angle positions.
%
%   The flag may be one of:
%
%     'P3_b'     Data from Fig.9; listener: P3, condition: 'baseline'.
%              
%     'P3_2o'    Data from Fig.9; listener: P3, condition: '2-oct'.
%              
%     'P3_1ol'   Data from Fig.9; listener: P3, condition: '1-oct(low)'.
%              
%     'P3_1om'   Data from Fig.9; listener: P3, condition: '1-oct(middle)'.
%              
%     'P3_1oh'   Data from Fig.9; listener: P3, condition: '1-oct(high)'.
%              
%     'P6_b'     Data from Fig.9; listener: P6, condition: 'baseline'.
%              
%     'P6_2o'    Data from Fig.9; listener: P6, condition: '2-oct'.
%              
%     'P6_1ol'   Data from Fig.9; listener: P6, condition: '1-oct(low)'.
%              
%     'P6_1om'   Data from Fig.9; listener: P6, condition: '1-oct(middle)'.
%              
%     'P6_1oh'   Data from Fig.9; listener: P6, condition: '1-oct(high)'.
%              
%     'P3_dtf'   Precalculated DTF data of P3 from Fig.11.
%              
%     'P6_dtf'   Precalculated DTF data of P6 from Fig.11.
%
%     'P3_dtf_bmp' DTFs calculated of P3 from the bitmap of the JASA paper.
%
%     'P6_dtf_bmp' DTFs calculated of P6 from the bitmap of the JASA paper.
%
%     'expdata'  Create the whole dataset required for exp_langendijk2002,
%                e.g. after adjusting response data.
%                *BE ADVISED**: The calculation takes a long time.
%
%   If no flag is given, the function will print the list of valid flags.
%
%   References:
%     E. Langendijk and A. Bronkhorst. Contribution of spectral cues to human
%     sound localization. J. Acoust. Soc. Am., 112:1583-1596, 2002.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_langendijk2002.php

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

% TODO: explain Data in description;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-19
% Bugfixes and minor adjustments: Sebastian Grill 2011-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  definput.flags.type={...
      'missingflag',...
      'P3_b','P3_2o','P3_1ol','P3_1om','P3_1oh',...
      'P6_b','P6_2o','P6_1ol','P6_1om','P6_1oh',...
      'P3_dtf','P6_dtf','P3_dtf_bmp','P6_dtf_bmp','expdata'
                      };
  % Parse input options
  [flags,keyvals]  = ltfatarghelper({},definput,varargin);

  if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
               sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
  end;

  
  target=-55:29:235;
  idb=round(0.5:0.1:11.4);
  idc=round(0.5:0.2:11.4);

  if flags.do_P3_b
    target=target(idb);
    response=zeros(1,110);
    response(1:10) =[-55,-55,-55,-55,-55,-55,-53,-54,-45,-45];
    response(11:20)=[-55,-40,-38,-34,-30,-29,-28,-22,-25,-15];
    response(21:30)=[-12,00,03,05,09,15,20,24,30,140];
    response(31:40)=[21,29,31,32,34,45,58,75,110,122];
    response(41:50)=[118,119,123,134,135,140,163,174,176,210];
    response(51:60)=[106,107,123,124,137,162,183,210,219,225];
    response(61:70)=[-55,95,113,119,120,121,126,161,181,214];
    response(71:80)=[119,130,142,143,170,184,190,199,225,235];
    response(81:90)=[146,153,154,155,157,161,180,181,181,210];
    response(91:100)=[93,94,170,180,182,183,184,191,230,231];
    response(101:110)=[230,231,235,235,235,235,235,235,235,234];
    data=[target;response];
  end;

  if flags.do_P3_2o
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-45,-12,165,181,182];
    response(6:10) =[-53,-30,-12,-11,192];
    response(11:15)=[-20,128,141,170,213];
    response(16:20)=[-55,-54,130,175,204];
    response(21:25)=[-30,72,172,180,202];
    response(26:30)=[164,181,219,235,234];
    response(31:35)=[-55,108,184,210,235];
    response(36:40)=[116,182,183,184,210];
    response(41:45)=[113,165,180,181,225];
    response(46:50)=[-54,-41,134,169,184];
    response(51:55)=[-55,234,234,235,235];
    data=[target;response];
  end;

  if flags.do_P3_1ol
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-54,-50,-49,-35,-30];
    response(6:10) =[-40,-37,-30,-23,-12];
    response(11:15)=[-25,-17,-3,5,13];
    response(16:20)=[-37,40,39,108,110];
    response(21:25)=[148,161,174,175,190];
    response(26:30)=[-50,-35,103,150,176];
    response(31:35)=[125,126,130,147,168];
    response(36:40)=[-55,110,167,168,175];
    response(41:45)=[144,152,154,180,190];
    response(46:50)=[173,179,180,215,225];
    response(51:55)=[-55,235,234,235,235];
    data=[target;response];
  end;

  if flags.do_P3_1om
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-55,-55,-54,-50,-48];
    response(6:10) =[-35,-28,-27,-25,-24];
    response(11:15)=[-51,-48,-25,-23,-18];
    response(16:20)=[-55,112,210,234,235];
    response(21:25)=[122,122,160,172,185];
    response(26:30)=[122,138,208,213,223];
    response(31:35)=[122,163,176,219,235];
    response(36:40)=[-55,-54,177,234,235];
    response(41:45)=[-55,128,163,180,224];
    response(46:50)=[-55,-54,150,151,183];
    response(51:55)=[-55,-54,228,234,235];
    data=[target;response];
  end;
  
  if flags.do_P3_1oh
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-55,-54,-50,-42,235];
    response(6:10) =[-55,-50,-40,-22,175];
    response(11:15)=[-29,3,60,145,150];
    response(16:20)=[14,130,166,168,180];
    response(21:25)=[40,122,123,148,180];
    response(26:30)=[122,123,128,208,215];
    response(31:35)=[119,133,138,145,175];
    response(36:40)=[122,175,180,181,182];
    response(41:45)=[110,175,175,176,190];
    response(46:50)=[-55,174,184,200,201];
    response(51:55)=[234,234,235,235,235];
    data=[target;response];
  end;

  if flags.do_P6_b
    target=target(idb);
    response=zeros(1,110);
    response(1:10) =[-55,-55,-55,-55,-55,-55,-55,-55,-55,-55];
    response(11:20)=[-50,-46,-41,-35,-36,-31,-32,-30,-25,-10];
    response(21:30)=[-25,-13,-8,-1,0,1,3,5,5,10];
    response(31:40)=[15,20,28,28,29,31,33,33,36,50];
    response(41:50)=[55,55,60,64,73,74,80,82,85,95];
    response(51:60)=[90,90,91,93,94,100,101,103,120,130];
    response(61:70)=[70,77,90,90,91,100,105,130,130,131];
    response(71:80)=[-50,-50,-40,0,95,100,101,105,162,167];
    response(81:90)=[145,146,150,160,165,169,170,171,180,181];
    response(91:100)=[185,186,195,196,200,209,210,211,215,220];
    response(101:110)=[210,215,220,225,230,231,235,235,234,234];
    data=[target;response];
  end;

  if flags.do_P6_2o
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-55,-54,-30,110,165];
    response(6:10) =[-55,-45,-40,-35,-20];
    response(11:15)=[-50,-35,-25,-20,125];
    response(16:20)=[-55,-50,-30,-31,180];
    response(21:25)=[-45,125,145,180,185];
    response(26:30)=[-55,-45,-30,-15,140];
    response(31:35)=[-55,-50,-45,-44,-20];
    response(36:40)=[-45,-40,-41,-35,170];
    response(41:45)=[-50,-35,100,180,181];
    response(46:50)=[-55,-54,-50,-45,140];
    response(51:55)=[-55,-55,-54,-54,-53];
    data=[target;response];
  end;
  
  if flags.do_P6_1ol
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-55,-54,-53,-49,-50];
    response(6:10) =[-35,-34,-30,-20,-19];
    response(11:15)=[-25,-10,-9,-11,5];
    response(16:20)=[25,30,31,35,40];
    response(21:25)=[65,70,95,96,100];
    response(26:30)=[60,85,90,91,110];
    response(31:35)=[75,85,90,100,105];
    response(36:40)=[-35,-25,-26,-20,80];
    response(41:45)=[0,185,184,186,190];
    response(46:50)=[175,195,205,210,215];
    response(51:55)=[225,230,233,234,235];
    data=[target;response];
  end;
  
  if flags.do_P6_1om
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-55,-54,-55,-54,-50];
    response(6:10) =[-55,-45,-44,-46,-35];
    response(11:15)=[-15,-10,5,6,10];
    response(16:20)=[5,15,20,30,31];
    response(21:25)=[45,60,61,65,66];
    response(26:30)=[-45,60,100,120,121];
    response(31:35)=[-55,80,105,175,176];
    response(36:40)=[-50,-40,-41,-10,10,];
    response(41:45)=[-55,-54,-50,-30,195];
    response(46:50)=[-50,100,105,106,135];
    response(51:55)=[-55,-55,-54,-54,235];
    data=[target;response];
  end;
  
  if flags.do_P6_1oh
    target=target(idc);
    response=zeros(1,55);
    response(1:5)  =[-55,-54,-55,-54,-35];
    response(6:10) =[-55,-54,-50,-51,-35];
    response(11:15)=[-15,-11,-10,-9,-5];
    response(16:20)=[10,25,35,55,105];
    response(21:25)=[40,45,50,51,75];
    response(26:30)=[-10,70,72,80,102];
    response(31:35)=[90,95,100,101,105];
    response(36:40)=[-55,-25,40,105,145];
    response(41:45)=[-55,-50,-45,-30,-25];
    response(46:50)=[-35,-34,-5,162,180];
    response(51:55)=[-55,-50,-54,235,234];
    data=[target;response];
  end;

  if flags.do_P3_dtf
      data=amt_load('langendijk2002','P3_dtf.mat');
  end

  if flags.do_P3_dtf_bmp
      [med,pol]=bmp2gr(amt_load('langendijk2002','P3_dtf.bmp'));
      data=[pol;med];
  end  
  
  if flags.do_P6_dtf
      data=amt_load('langendijk2002','P6_dtf.mat');
  end

  if flags.do_P6_dtf_bmp
      [med,pol]=bmp2gr(amt_load('langendijk2002','P6_dtf.bmp'));
      data=[pol;med];
  end  
  
  if flags.do_expdata
    listener='P3';
    amt_disp('Calculation for P3','progress');
    temp=data_langendijk2002([listener '_dtf']);
    fs=temp.stimPar.SamplingRate;           
    pol=temp.posLIN(:,5)';
    temp=temp.ampMdB;
    med=temp(1:end,:);
    amt_disp('  Condition: b','progress');
    temp=data_langendijk2002([listener '_b']);
    targetb=temp(1,:); responseb=temp(2,:);
    medir=gr2ir(med,'b',fs);
    amt_disp('  Condition: 2o','progress');
    temp=data_langendijk2002([listener '_2o']);
    targetc=temp(1,:); response2o=temp(2,:);
    medir2o=gr2ir(med,'2o',fs);
    amt_disp('  Condition: 1ol','progress');
    temp=data_langendijk2002([listener '_1ol']);
    response1ol=temp(2,:);
    medir1ol=gr2ir(med,'1ol',fs);
    amt_disp('  Condition: 1om','progress');
    temp=data_langendijk2002([listener '_1om']);
    response1om=temp(2,:);
    medir1om=gr2ir(med,'1om',fs);
    amt_disp('  Condition: 1oh','progress');
    temp=data_langendijk2002([listener '_1oh']);
    response1oh=temp(2,:);
    medir1oh=gr2ir(med,'1oh',fs);
    save(fullfile(amt_auxdatapath, 'langendijk2002', [listener '_data.mat']),'-regexp','^med|^pol|^response|^target');
    
    clear
    amt_disp('Calculation for P6','progress');
    listener='P6';
    temp=data_langendijk2002([listener '_dtf']);
    fs=temp.stimPar.SamplingRate;
    pol=temp.posLIN(:,5)';
    temp=temp.ampMdB;
    med=temp(1:end,:);
    amt_disp('  Condition: b','progress');
    temp=data_langendijk2002([listener '_b']);
    targetb=temp(1,:); responseb=temp(2,:);
    medir=gr2ir(med,'b',fs);
    amt_disp('  Condition: 2o','progress');
    temp=data_langendijk2002([listener '_2o']);
    targetc=temp(1,:); response2o=temp(2,:);
    medir2o=gr2ir(med,'2o',fs);
    amt_disp('  Condition: 1ol','progress');
    temp=data_langendijk2002([listener '_1ol']);
    response1ol=temp(2,:);
    medir1ol=gr2ir(med,'1ol',fs);
    amt_disp('  Condition: 1om','progress');
    temp=data_langendijk2002([listener '_1om']);
    response1om=temp(2,:);
    medir1om=gr2ir(med,'1om',fs);
    amt_disp('  Condition: 1oh','progress');
    temp=data_langendijk2002([listener '_1oh']);
    response1oh=temp(2,:);
    medir1oh=gr2ir(med,'1oh',fs);
    save(fullfile(amt_auxdatapath, 'langendijk2002', [listener '_data.mat']),'-regexp','^med|^pol|^response|^target');
  end

end

function [med,pol]=bmp2gr(bmpname)
% BMP2GR converts a bitmap to gain response data (in dB) with conditions
% Usage:  [med,pol]=bmp2gr(bmpname)
% Input arguments:
%       bmpname:    filename of bitmap
% Output arguments:
%       medir:   	gain response data (in dB) on median plane for 
%                   53 polar angles defined in pol
%       pol:     	53 equally spaced polar angles between -55° and 235°
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bmp=imread(bmpname,'bmp');
hrtf=zeros(size(bmp,2),size(bmp,1));

sq=4;
for jj=1:size(hrtf,2)
    if jj<=sq || jj>=size(hrtf,2)-sq
        av=0;
    else
        av=sq;
    end
    for ii=1:size(hrtf,1)
        temp=mean(bmp(end-jj+1-av:end-jj+1+av,ii,1));
        if temp<=25 && temp >=15
            hrtf(ii,jj)=-15;
        elseif temp<=50 && temp>40
            hrtf(ii,jj)=-10;
        elseif temp<=85 && temp>65
            hrtf(ii,jj)=-5;
        elseif temp<=150 && temp>120
            hrtf(ii,jj)=0;
        elseif temp<=185 && temp>160
            hrtf(ii,jj)=5;
        elseif temp<=210 && temp>199
            hrtf(ii,jj)=10;
        elseif temp>220
            hrtf(ii,jj)=15;
        else
            if ii==1
                hrtf(ii,jj)=0;
            else
                hrtf(ii,jj)=hrtf(ii-1,jj);
            end
        end
    end
end

sq=3;
for jj=1:size(hrtf,2)
    if jj<=sq || jj>=size(hrtf,2)-sq
        av=0;
    else
        av=sq;
    end
    for ii=1:size(hrtf,1)
        hrtf(ii,jj)=mean(hrtf(ii,jj-av:jj+av));
    end
end

% figure; pcolor(hrtf'),shading flat
posidx=5:9:475;
pol=-55:290/52:235;
med=hrtf(:,posidx);
% figure; pcolor(med'),shading flat

end

function [medir]=gr2ir(med,cond,fs)
% GR2IR converts given gain responses MED (in dB) to impulse responses
% MEDIR; furthermore several conditions according to langendijk et al.
% (2002) can be defined
% Usage:            [medir]=gr2ir(med,cond,fs)
% Input arguments:
%       med:        gain responses (in dB)
%       cond:       condition, 
%                   possibilities:  baseline    'b'
%                                   2 octaves   '2o'
%                                   1 oct (low) '1ol'
%                                   1 oct (mid) '1om'
%                                   1 oct (high)'1oh'
%       fs:      	sampling frequency
% Output arguments:
%       medir:   	impulse responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default settings
if ~exist('cond','var')
    cond='b';
end
if ~exist('fs','var')
    fs=48000;
end

imp=zeros(240,1);imp(1)=1;
n=256;
len=240;
frq=logspace(log10(2000),log10(16000),size(med,1));

% frequency indices
f1=find(frq>=4000,1);
f2=find(frq>=5700,1);
f3=find(frq>=8000,1);
f4=find(frq>=11300,1);

switch cond
    case 'b' % baseline
        medir=zeros(n,size(med,2));
        for ii=1:size(med,2)
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '2o' % 2 oct (4-16kHz)
        medir=zeros(n,size(med,2));
        med2o=med;
        for ii=1:size(med,2)
            med2o(f1:end,ii)=mean(med(f1:end,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med2o(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1ol' % 1 oct (low:4-8kHz)
        medir=zeros(n,size(med,2));
        med1ol=med;
        for ii=1:size(med,2)
            med1ol(f1:f3,ii)=mean(med(f1:f3,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1ol(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1om' % 1 oct (middle:5.7-11.3kHz)
        medir=zeros(n,size(med,2));
        med1om=med;
        for ii=1:size(med,2)
            med1om(f2:f4,ii)=mean(med(f2:f4,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1om(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1oh' % 1 oct (high:8-16kHz)
        medir=zeros(n,size(med,2));
        med1oh=med;
        for ii=1:size(med,2)
            med1oh(f3:end,ii)=mean(med(f3:end,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1oh(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end
end
end

