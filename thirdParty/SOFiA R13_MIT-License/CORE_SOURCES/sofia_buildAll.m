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
% sofia_buildAll.m (Automatic Builder)
% 
% Copyright 2011-2017 Benjamin Bernschütz, rockzentrale 'AT' me.com  
%                        and Nils Peters, nils 'AT' icsi.berkeley.edu 
% 
% This file is part of the SOFiA toolbox under MIT-License
% 
% void = sofia_buildAll(configuration)
% 
% configuration: [] or 'DEBUG' 
%
% This file compiles the SOFiA C++ sources. Run the file inside the
% \CORE_SOURCES folder.  
%
% Remember: You need a C++ Compiler and need to configure MATLAB to 
%           use this compiler. Further the BOOST C++ ist required.
%           Read "howtocompile.txt" for more information.
%
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
%     Nils Peters                  nils 'at' icsi.berkeley.edu


function sofia_buildAll(configuration)
clc

if nargin == 0
    configuration = 'RELEASE';
end

disp(' ');
disp('*** building SOFiA ***');
disp(' ');

sourcefiles = dir('*.cpp');
sources=[];
for i=1:size(sourcefiles,1)
    sources{i}  = sourcefiles(i).name;
    target_m{i} = strrep(sourcefiles(i).name,'.cpp','.m');
end

cHeaderSources = dir(['HEADER',filesep(),'*.cpp']);
cHeaders=[];

for i=1:size(cHeaderSources,1)
    cHeaders = [cHeaders, 'HEADER',filesep(),cHeaderSources(i).name,' '];
end

cHeaderPath = [pwd(),filesep(),'HEADER'];


for k = 1:length(sources);
    eval(sprintf('disp(''building %s'');',char(sources(k))));
    % mex compiling
    sprintf('mex %s %s -I''%s'' -D%s -outdir ../;',char(sources(k)),char(cHeaders), char(cHeaderPath),char(configuration))
    eval(sprintf('mex %s %s -I''%s'' -D%s -outdir ../;',char(sources(k)),char(cHeaders), char(cHeaderPath),char(configuration)));  
        
    % header file generation: writing the header into an extra .m file 
    eval(sprintf('fir=fopen(''%s'' ,''r'',''n'',''ISO-8859-1'');',char(sources(k)))); 
    eval(sprintf('fiw=fopen(''../%s'',''w'',''n'',''ISO-8859-1'');',char(target_m(k))));
    tline = fgetl(fir); tline = fgetl(fir);
    while 1
        tline = fgetl(fir);
        if isempty(tline)
            fprintf(fiw,' \n');
            break;
        elseif char(tline(1)) == '@'
            fprintf(fiw, ' \n');        
        elseif char(tline(1)) == '%'
            fprintf(fiw, '%s\n',char(tline));
        else                
            break;
        end;
    end;            
 fclose(fir);
 fclose(fiw);
 end;

disp(' ');    
disp('***     done       ***');
disp(' ');

