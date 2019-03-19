function [V,Y,E,CF]=verhulst2012(sign,fs,fc,spl,normalizeRMS,subject,irregularities)
%VERHULST2012 Process a signal with the cochlear model by Verhulst et. al. 2012
%   Usage: output = verhulst2012(insig,fs,fc,spl,normalizeRMS,subject,irregularities,sheraPo)
%
%   Input parameters:
%        sign           : the input signal to be processed. Each column is processed
%                         in parallel, so it is possible to run several
%                         simulation in parallel
%        fs             : sampling rate
%        fc             : list of frequencies specifying the probe positions on
%                         the basilar membrane, or 'all' to probe all 1000
%                         cochlear sections
%        spl            : array of the sound pressure levels that correspond to
%                         value 1 of the correspondent input signal
%        normalizeRms   : arry to control the normalization of each
%                         signal. With value 1 normalize the energy of the signal, so the
%                         relative spl value correspond to the rms of the signal (default 0)
%        subject        : the subject number controls the cochlear irregulatiries (default 1)
%        irregularities : array that enable (1) or disable (0)
%                         irregularities and nonlinearities for each simulation (default 1)                
%       
%   Output parameters:
%       V  : velocity of the basilar membrane sections V(time,section,channel) 
%       Y  : displacement of the basilar membrane sections Y(time,section,channel)
%       E  : sound pressure at the middle ear
%       CF : center frequencies of the probed basiliar membrane sections
%
%   This function computes the basilar membrane displacement and the
%   velocity of the movement at different positions employing a faster
%   implementation of the nonlinear time-domain model of cochlea by
%   Verhulsts, Dau, Shera 2012, through the method described in Alto? et
%   al. 2014
%
%   The processing is implemented as follows:
%
%   1) the input signal is resampled to the 96 kHz sampling rate
%      employed in the cochlea model
%
%   2) the list of frequencies in fc are converted in to probe
%      positions in a manner that the frequencies are divided evenly into
%      low and high frequency categories. 
%
%   3) the signals are processed in parallel
%
%   4) the values obtained are resampled back to the original sampling
%      rate
%
%   Requirements and installation: 
%   ------------------------------
%
%   1) Python >2.6 is required with numpy and scipi packages. On Linux, use sudo apt-get install python-scipy python-numpy
% 
%   2) Compiled files with a C-compiler, e.g. gcc. In amtbase/src/verhulst start make (Linux) or make.bat (Windows)
%
%   3) On linux, when problems with GFORTRAN lib appear, try sudo ln -sf /usr/lib64/libgfortran.so.3.0.0 /mymatlabroot/sys/os/glnxa64/libgfortran.so.3 (mymatlabroot is usually /usr/local/MATLAB/version
%               
%   See also: verhulst2012,
%
%   References:
%     S. Verhulst, T. Dau, and C. A. Shera. Nonlinear time-domain cochlear
%     model for transient stimulation and human otoacoustic emission. J.
%     Acoust. Soc. Am., 132(6):3842 - 3848, 2012.
%     
%     A. Altoè, S. Verhulst, and V. Pulkki. Transmission line cochlear
%     models: improved accuracy and efficiency. J. Acoust. Soc. Am.,
%     136(EL302), 2014.
%     
%
%   AUTHOR: Alessandro Altoe' 
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/verhulst2012.php

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

if nargin<5
    [channels,idx]=min(size(sign));
    normalizeRMS=zeros(channels,1);
end
if nargin < 6
    subject = 1;
end
if nargin<7
    [channels,idx]=min(size(sign));
    irregularities=ones(1,channels);

end
% sign=double(sign); %force to write signal as double array
Fs=96000;
sectionsNo=1000;
[channels,idx]=min(size(sign));
if(idx==2) %transpose it (python C-style row major order)
    sign=sign';
end
stim=zeros(channels,length(resample(sign(1,:),Fs,fs)));
for i=1:channels
    stim(i,:)=resample(sign(i,:),Fs,fs);
%     figure; plot(stim(i,:));
    if normalizeRMS(i)
        s_rms=rms(stim(i,:));
        stim(i,:)=stim(i,:)./s_rms;
    end
end
sheraPo=0.061;
if(isstr(fc) && strcmp(fc,'all')) %if probing all sections 1001 output (1000 sections plus the middle ear)
    p=sectionsNo;
else %else pass it as a column vector
    [p,idx]=max(size(fc));
    if(idx==2)
        fc=fc'; 
    end
    fc=round(fc);
end
probes=fc;
[path,name,ext]=fileparts(which('verhulst2012'));
[path,name,ext]=fileparts(path);
act_path=pwd;
cd(strcat(path,'/bin/verhulst2012/')); 
if ~exist('tridiag.so','file')
	error('tridiag.so library is missing. Goto to AMT/bin/verhulst2012 and compile by executing the makefile');
end
save('input.mat','stim','Fs','channels','spl','subject','sheraPo','irregularities','probes','-v7');
system('python run_cochlear_model.py');
l=length(stim(1,:));
rl=length(resample(stim(1,:),fs,Fs));
V=zeros(rl,p,channels);
Y=zeros(rl,p,channels);
E=zeros(rl,channels);
CF=zeros(p,1);
for i=1:channels
    fname=strcat('out/v',int2str(i),'.np');
    f=fopen(fname,'r');
    V(:,:,i)=resample(fread(f,[p,l],'double','n')',fs,Fs);
    fclose(f);
    fname=strcat('out/y',int2str(i),'.np');
    f=fopen(fname,'r');
    Y(:,:,i)=resample(fread(f,[p,l],'double','n')',fs,Fs);
    fclose(f);
    fname=strcat('out/E',int2str(i),'.np');
    f=fopen(fname,'r');
    E(:,i)=resample(fread(f,[l,1],'double','n'),fs,Fs);
    fclose(f);
    if(i==1)
    fname=strcat('out/F',int2str(i),'.np');
    f=fopen(fname,'r');
    CF=fread(f,[p,1],'double','n');
    fclose(f);
    end
end
cd(act_path);

