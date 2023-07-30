function amt_mex(varargin)
%AMT_MEX   Compile Mex/Oct interfaces
%   Usage:  amt_mex;
%           amt_mex(flags);
%
%   AMT_MEX compiles AMT-related binaries on your system. 
%   
%   Requirements
%   ------------
%   
%   1) Mex compiler working in your Matlab/Octave environment. Type help mex to check whether
%      and which compiler is available in your environment.
%
%   2) GCC installed on your OS. 
%      On Windows, GCC can be installed from <https://gcc.gnu.org>.
%      On Linux, GCC is most probably available within your distribution. 
%
%   Details
%   -------
%
%   The action of AMT_MEX is determined by one of the following flags:
%
%     'compile'  Compile stuff. This is the default.
%                In Matlab, all comp_.c and comp_.cpp files from theare
%                mex directory are compiled to system-dependent mex files. 
%                In Octave, all comp_.cc files from the oct directory
%                are compiled to oct files. Then, the remaining files from 
%                the mex directory are compiled to mex and moved to oct.
%                In both environments, other binaries are handled then. 
%                On Windows, make.bat from the bin directory is executed.
%                On Linux, make is executed using the makefile file from the bin directory.
%
%     'clean'    Removes the compiled functions.
%                In Matlab, all system-dependent mex files from the mex
%                directory are removed.
%                In Octave, all oct, o, and mex files from the oct
%                directory are removed.
%                In both environments, other binaries are cleared by calling 
%                clean and make clean on Windows and other systems, respectively, 
%                in the bin directory. 
%
%   See also: amt_start
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_mex.php


%   #Author : Peter Søndergaard (2015)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

[~, kv] = amt_configuration;
bp = [kv.path, '/'];

defnopos.flags.command={'compile','clean'};
[flags,kv]=ltfatarghelper({},defnopos,varargin);

% Remember the current directory.
curdir=pwd;

if isoctave
  extname='oct';
else
  extname='mex';
end;

% -------------- Handle cleaning --------------------------------
  
if flags.do_clean
    
  if ~isoctave
    % Delete files permanently (bypass trashcan on Windows/Mac
    % but remember the old state
    oldstate = recycle('off');
  end;
  
  amt_disp(['========= Cleaning ',extname,' interfaces ==========']);
  if isoctave
    deletefiles([bp,'oct'],'*.oct');
    deletefiles([bp,'oct'],'*.o');
    deletefiles([bp,'oct'],'*.mex');
  else
    deletefiles([bp,'mex'],['*.',mexext]);
  end;

  amt_disp('========= Cleaning environment interfaces ==========');
  cd([bp,'environments']); 
  if ispc, 
    [~, output]=system('clean'); 
  else
    [~, output]=system('make clean');
  end
  amt_disp(output);
  
  if ~isoctave
    recycle(oldstate);
  end;
  % Jump back to the original directory.
  cd(curdir); 
end;

% -------------- Handle compiling  --------------------------------

if flags.do_compile

  amt_disp(['========= Compiling ',extname,' interfaces ==========']);
  % Get the list of files.
  if isoctave
    ext='oct';
    L=dir([bp,filesep,'oct',filesep,'comp_*.cc']);
    L=[L; dir([bp,filesep,'oct',filesep,'comp_*.cpp'])];
  else
    ext=mexext;
    L=dir([bp,filesep,'mex',filesep,'comp_*.c']);
    L=[L; dir([bp,filesep,'mex',filesep,'comp_*.cpp'])];
  end;
  filenames = arrayfun(@(lEl) lEl.name,L,'UniformOutput',0);
  
  if compile_amt(bp,ext,filenames)>1
    amt_disp(['Error: The ',extname,' interfaces was not built.']);
  else
    amt_disp('Done.');
  end;
  
  if isoctave
    % Compile MEXs instead of missing OCTs
    Lmex=dir([bp,filesep,'mex',filesep,'comp_*.c']); 
    mexnamesstrip = arrayfun(@(lEl) lEl.name(1:end-2),Lmex,'UniformOutput',0);
    octnamesstrip = cellfun(@(lEl) lEl(1:end-3),filenames,'UniformOutput',0);
    
    mexdiffstrip = setdiff(mexnamesstrip,octnamesstrip);
    
    mexdiff = cellfun(@(lEl) [lEl,'.c'],mexdiffstrip,'UniformOutput',0);
    if ~isempty(mexdiff)
        amt_disp('========= Compiling MEX interfaces ==========')
        if compile_amt(bp,'mex',mexdiff)>1;                
            amt_disp(['Error: The ',extname,' interfaces was not built.']);
        else
            if movefile([bp,filesep,'mex',filesep,'*.mex'],...
                        [bp,filesep,'oct'],'f');
               amt_disp('Done.');
            else
               error(['Error: Compilation sucessful, but MEX files were not '...
               'moved from mex to oct directory. Check your write permissions.\n']); 
            end
        end;
    end
 end;

  amt_disp('========= Compiling environments interfaces ==========');
  cd([bp,'environments']); 
  [~, output]=system('make');  
  amt_disp(output);
  amt_disp('Done.');
 
% Jump back to the original directory.
cd(curdir);
end


function deletefiles(base,files)

L=dir([base,filesep,files]);
for ii=1:numel(L)
    s=[base,filesep,L(ii).name];
    delete(s);
end;



function status=compile_amt(bp,ext,filenames)

% If we exit early, it is because of an error, so set status=1
status=1;

    if strcmpi(ext(1:3),'oct')
        cd([bp,'oct']); 
    else
        cd([bp,'mex']);
    end

for ii=1:numel(filenames)
    filename = filenames{ii};
    dotPos = strfind(filename,'.');
    objname  = [filename(1:dotPos(end)),ext];
    objdirinfo = dir(objname);
    
    % Make-like behaviour: build only the files where the src file is
    % newer than the object file, or the object file is missing.
    L = dir(filename);
    if isempty(objdirinfo) || (objdirinfo.datenum<L(1).datenum)
        
        amt_disp(['Compiling ', filename]);
        if isoctave
          if ~strcmpi(ext(1:3),'oct')
              mkoctfile('-mex','-I.','-I../src',filename);
          else
              mkoctfile('-I.','-I../src',filename);
          end
        else
          mex('-I.','-I../src',filename);
        end;                
        
    end;        
    status=0;
end;





                    

    
    
    
    




