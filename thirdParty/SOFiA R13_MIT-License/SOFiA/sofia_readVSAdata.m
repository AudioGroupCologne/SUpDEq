% /// ASAR/MARA Research Group
%  
% Technology Arts Sciences TH Köln
% Technical University of Berlin
% Deutsche Telekom Laboratories
% University of Rostock
% WDR Westdeutscher Rundfunk
% IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis
% 
% VariSphear -> SOFiA data import R13-0215
% 
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%              
% 
% This file is part of the SOFiA toolbox under MIT-License
%
%
% [timeDataCH1, timeDataCH2] = readVSA_data([downSample], ...
%                                    [normalize], [directory])
% ------------------------------------------------------------------------
% timeDataCH1/CH2     Structs with fields:
%
%                     .impulseResponses     [Channels x Samples]
%                     .centerIR             [1 x Samples]
%                     .FS
%                     .quadratureGrid       [AZ1 EL1 W1; ...; AZn ELn Wn]
%                     .metaData             Cell Array, VSA Metadata
%                     .downSample
%                     .averageAirTemp       Temperature in DEG
%                     .irOverlay            Plot this for a good total
%                                           overview of the dataset.
% ------------------------------------------------------------------------
% downSample         Downsampling factor  [default = 1: No downsampling]
%                    Downsampling is done using DECIMATE and a FIR low
%                    pass filter of order 30. See MATLAB documentation
%                    for more information.
%                    !!! MATLAB Signal Processing Library required
%
% normalize          Normalize flag 1:on, 0:off         [default = 1: on]
%                    Normalizes the impulse responses with respect to the
%                    absolute maximum value within the complete dataset.
%
% directory          VSA dataset directory. If not defined a user dialog
%                    opens to pick a directory.
%

% CONTACT AND LICENSE INFORMATION:
% 
% /// ASAR/MARA Research Group 
%  
%     [1] Technology Arts Sciences TH Köln
%     [2] Technical University of Berlin 
%     [3] Deutsche Telekom Laboratories 
%     [4] University of Rostock
%     [5] WDR Westdeutscher Rundfunk 
%     [6] IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis toolbox
% 
% Copyright 2011-2017 Benjamin Bernschütz et al.(§)  
% 
% Contact ------------------------------------
% Technology Arts Sciences TH Köln 
% Institute of Communications Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
% 
% phone       +49 221 8275 -2496 
% cell phone  +49 171 4176069 
% mail        rockzentrale 'at' me.com 
% --------------------------------------------
% 
% This file is part of the SOFiA sound field analysis toolbox
%
% Licence Type: MIT License
%
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
%
% (§) Christoph Pörschmann [1]     christoph.poerschmann 'at' th-koeln.de
%     Sascha Spors         [2,3,4] sascha.spors 'at' uni-rostock.de  
%     Stefan Weinzierl     [2]     stefan.weinzierl 'at' tu-berlin.de

function [timeDataCH1, timeDataCH2] = sofia_readVSAdata(downSample, normalize, directory)

disp('SOFiA VariSphear -> SOFiA data import R13-0215');

air_temperature = 20; %DEFAULT if no temperatures found in vsd-dataset

if nargin == 0
    downSample = 1;
end

if nargin < 2
    normalize = 1;
end

if nargin < 3
    directory = uigetdir();
    if directory == 0
        error('No directory picked');
    end
end

files=dir(fullfile(directory, '*.mat'));
arraycounterCH1 = 0;
arraycounterCH2 = 0;
indexListCH1    = 0;
indexListCH2    = 0;

Temp1=[];
Temp2=[];

for filecounter=1:size(files,1)
    
    imported = files(filecounter).name;
    vsd=[];
    
    if ~isempty(imported)
        if mod(filecounter,5) == 0
            fprintf('|');
        end
        if mod(filecounter, 200) == 0
            fprintf('\n');
        end
               
        load(fullfile(directory, files(filecounter).name), 'vsd')
        
        if ~isempty(vsd)
            
            if vsd.CH == 1 || vsd.HalfSphere == 1% -------------------------------------------- CH1
                
                if strcmp(vsd.StorageType,'MAT')
                    if downSample ~= 1
                        if arraycounterCH1 == 0
                            ir_length = length(decimate(cast(vsd.ImpulseResponse,'double'),downSample ,'FIR'));
                            timeDataCH1.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length);      
                            indexListCH1 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH1.impulseResponses(vsd.SamplePosNo,:)=decimate(cast(vsd.ImpulseResponse,'double'),downSample ,'FIR');
                        timeDataCH1.FS=vsd.FS/downSample;
                    else
                        if arraycounterCH1 == 0
                            ir_length = length(vsd.ImpulseResponse);
                            timeDataCH1.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length); 
                            indexListCH1 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH1.impulseResponses(vsd.SamplePosNo,:)=vsd.ImpulseResponse;
                        timeDataCH1.FS=vsd.FS;
                    end
                else
                    try                 
                        irdata = wavread(fullfile(directory,vsd.ImpulseResponse(5:end)))';
                    catch
                        error(['ERROR - Data file missing: ', vsd.ImpulseResponse(5:end)]);
                    end
                    
                    if downSample ~= 1
                        if arraycounterCH1 == 0
                            ir_length = length(decimate(cast(irdata,'double'),downSample ,'FIR'));
                            timeDataCH1.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length); 
                            indexListCH1 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH1.impulseResponses(vsd.SamplePosNo,:) = decimate(cast(irdata,'double'),downSample ,'FIR');
                        timeDataCH1.FS=vsd.FS/downSample;
                    else
                        if arraycounterCH1 == 0
                            ir_length = length(irdata);
                            timeDataCH1.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length); 
                            indexListCH1 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH1.impulseResponses(vsd.SamplePosNo,:) = cast(irdata,'double');
                        timeDataCH1.FS=vsd.FS;
                    end
                end
                
                timeDataCH1.metaData{vsd.SamplePosNo}=rmfield(vsd,'ImpulseResponse');
                timeDataCH1.radius = vsd.Radius;
                timeDataCH1.quadratureGrid(vsd.SamplePosNo,1)=vsd.Azimuth*pi/180;
                timeDataCH1.quadratureGrid(vsd.SamplePosNo,2)=vsd.Elevation*pi/180;
                timeDataCH1.quadratureGrid(vsd.SamplePosNo,3)=vsd.GridWeight;
                
                if isfield(vsd,'Temperature')
                    Temp1(arraycounterCH1+1)=vsd.Temperature;
                else
                    Temp1(arraycounterCH1+1)=air_temperature;
                end
                indexListCH1(vsd.SamplePosNo) = 0;   

                arraycounterCH1=arraycounterCH1+1;
                
            else   % -------------------------------------------- CH2
                
                if strcmp(vsd.StorageType,'MAT')
                    if downSample ~= 1
                        if arraycounterCH2 == 0
                            ir_length = length(decimate(cast(vsd.ImpulseResponse,'double'),downSample ,'FIR'));
                            timeDataCH2.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length);
                            indexListCH2 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH2.impulseResponses(vsd.SamplePosNo,:)=decimate(cast(vsd.ImpulseResponse,'double'),downSample ,'FIR');
                        timeDataCH2.FS=vsd.FS/downSample;
                    else
                        if arraycounterCH2 == 0
                            ir_length = length(vsd.ImpulseResponse);
                            timeDataCH2.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length);
                            indexListCH2 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH2.impulseResponses(vsd.SamplePosNo,:)=vsd.ImpulseResponse;
                        timeDataCH2.FS=vsd.FS;
                    end
                else
                    try
                        irdata = wavread(fullfile(directory,vsd.ImpulseResponse(5:end)))';                        
                    catch
                        error(['ERROR - Data file missing: ', vsd.ImpulseResponse(5:end)]);
                    end
                    
                    if downSample ~= 1
                        if arraycounterCH2 == 0
                            ir_length = length(decimate(cast(irdata,'double'),downSample ,'FIR'));
                            timeDataCH2.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length);
                            indexListCH2 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH2.impulseResponses(vsd.SamplePosNo,:) = decimate(cast(irdata,'double'),downSample ,'FIR');
                        timeDataCH2.FS=vsd.FS/downSample;
                    else
                        if arraycounterCH2 == 0
                            ir_length = length(irdata);
                            timeDataCH2.impulseResponses = zeros(vsd.TotalSamplepoints, ir_length);
                            indexListCH2 = ones(vsd.TotalSamplepoints,1);
                            clear ir_length
                        end
                        timeDataCH2.impulseResponses(vsd.SamplePosNo,:) = cast(irdata,'double');
                        timeDataCH2.FS=vsd.FS;
                    end
                end
                
                timeDataCH2.metaData{vsd.SamplePosNo}=rmfield(vsd,'ImpulseResponse');
                timeDataCH2.radius = vsd.Radius;
                timeDataCH2.quadratureGrid(vsd.SamplePosNo,1)=vsd.Azimuth*pi/180;
                timeDataCH2.quadratureGrid(vsd.SamplePosNo,2)=vsd.Elevation*pi/180;
                timeDataCH2.quadratureGrid(vsd.SamplePosNo,3)=vsd.GridWeight;
                
                if isfield(vsd,'Temperature')
                    Temp2(arraycounterCH2+1)=vsd.Temperature;
                else
                    Temp2(arraycounterCH2+1)=air_temperature;
                end        
                indexListCH2(vsd.SamplePosNo) = 0;
                arraycounterCH2=arraycounterCH2+1;
            end
        end
    end
end

% CenterIR
centerPath     = fullfile(directory,'CENTER');
centerLocated  = false;
centerOkay     = false;
centerIsInline = false;

if exist(centerPath,'dir')
    centerFiles     = dir(fullfile(centerPath, '*.mat'));
    averagedLocated = false;
    for ctFileCounter = 1:size(centerFiles,1)
        importedCT = centerFiles(ctFileCounter).name;
        if strfind(importedCT,'Averaged') & strfind(importedCT,'CTIR'); 
            averagedLocated = true;
            break
        end
    end
    if averagedLocated
        try
            load(fullfile(centerPath,importedCT),'vsd')
            centerVsd = vsd;
            centerLocated = true;
            centerIsInline = true;
        catch
        end
    end
else
    centerPath  = directory;
end

if ~centerLocated
    [importedCT, centerPath] = uigetfile(fullfile(centerPath,'*.mat; *.wav'),'Pick the center impulse response file');
    if importedCT
        try
            if strfind(importedCT,'.wav')
                centerIsInline = false;
                centerLocated  = true;
            else
                centerIsInline = true;
                load(fullfile(centerPath,importedCT),'vsd')
                centerVsd      = vsd;
                centerLocated  = true;
            end
        catch
        end
    else
        centerLocated  = false;
        centerOkay     = false;
    end
end

if centerLocated    
    if centerIsInline
        if strcmp(centerVsd.StorageType,'MAT')
            waveFile = [];
            if downSample ~= 1
                timeDataCH1.centerIR = decimate(cast(centerVsd.ImpulseResponse,'double'),downSample ,'FIR');
            else
                timeDataCH1.centerIR = centerVsd.ImpulseResponse;
            end
            centerOkay = true;
        else            
            waveFile = strrep(fullfile(centerPath, importedCT),'.mat','.wav');
        end
    else
        waveFile = fullfile(centerPath,importedCT);
    end
    
    if ~isempty(waveFile)
        if exist(waveFile,'file')
            irdata = wavread(waveFile)';
            centerOkay = true;
        else
            error(['Inline reference not found: ', waveFile]);
        end
        
        if downSample ~= 1
            timeDataCH1.centerIR = decimate(cast(irdata,'double'),downSample ,'FIR');
        else
            timeDataCH1.centerIR = cast(irdata,'double');
        end
    end
    
    if normalize == 1
        timeDataCH1.centerIR = timeDataCH1.centerIR./max(max(abs(timeDataCH1.centerIR)));
    end
    
    if arraycounterCH2>0
        timeDataCH2.centerIR = timeDataCH1.centerIR;
    end
end

if centerOkay
    fprintf(['\n\nCenterIR: ',importedCT]);
    if size(timeDataCH1.centerIR,2) ~= size(timeDataCH1.impulseResponses,2)
        fprintf(['\nWARNING: Impulse responses and CenterIR do not have the same length!']);
    end
    
else
    fprintf(['\n\nWARNING: No valid center file found > timeData.centerIR field is not written!']);    
end


if sum(indexListCH1) > 0
   fprintf(['\n\nWARNING: Folowing samples missing for CH1:']);
   for i = 1:length(indexListCH1)
       if indexListCH1(i) == 1
           fprintf(' %d',i);
       end
   end
end

if sum(indexListCH2) > 0
   fprintf(['\nWARNING: Folowing samples missing for CH2:']);
   for i = 1:length(indexListCH2)
       if indexListCH2(i) == 1
           fprintf(' %d',i);
       end
   end
end

if arraycounterCH1>0
    fprintf('\n\n');
    
    if normalize == 1
        timeDataCH1.impulseResponses = timeDataCH1.impulseResponses./max(max(abs(timeDataCH1.impulseResponses)));
    end
    
    timeDataCH1.downSample=downSample;
    timeDataCH1.averageAirTemp = mean(Temp1);
    timeDataCH1.irOverlay = sum(timeDataCH1.impulseResponses,1);
    timeDataCH1.irOverlay = abs(timeDataCH1.irOverlay/max(abs(timeDataCH1.irOverlay)));
    
    if arraycounterCH2>0
        disp(['CH1: ',num2str(arraycounterCH1),' spatial sampling points imported.']);
        disp(['CH2: ',num2str(arraycounterCH2),' spatial sampling points imported.']);
        if normalize == 1
            timeDataCH2.impulseResponses=timeDataCH2.impulseResponses./max(max(abs(timeDataCH2.impulseResponses)));
        end
        timeDataCH2.downSample=downSample;
        timeDataCH2.averageAirTemp = mean(Temp2);
        timeDataCH2.irOverlay = sum(timeDataCH2.impulseResponses,1);
        timeDataCH2.irOverlay = abs(timeDataCH2.irOverlay/max(abs(timeDataCH2.irOverlay)));
    else
        disp([num2str(arraycounterCH1),' spatial sampling points imported.']);
        timeDataCH2=[];
    end    
else
    disp(['\n\nNothing found to import.']);
    timeDataCH1=[];
    timeDataCH2=[];
end

fprintf('\n');
