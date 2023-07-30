function [out,status] = amt_extern(environment, directory, module, input, outstruct)
%AMT_EXTERN Process a function from an external environment
%   Usage: 
%     output = amt_extern(environment, directory, module, input, outstruct);
%
%   Input parameters:
%
%     environment  : String selecting the environment to be used.
%                    Currently must be Python.
%     directory    : Diretory inside the environment. 
%                    For example, verhulst2018 for the model from 
%                    Verhulst et al. (2018).
%     module       : Module name to be called within the environment.
%                    For example, run_cochlear_model.py for the model
%                    from Verhulst et al. (2018).
%     input        : Structure with the input parameters. All fields 
%                    will be saved in the directory out within the
%                    directory in the file input.mat.
%     outstruct    : Structure defining the output parameters. Each field's
%                    name defines the variable name to be read from a file
%                    named by the field with the ending .np in the 
%                    out directory. Each field must contain a vector 
%                    defining the size of the variable to be read. Up to 
%                    three dimensions are handled. 
%
%   Output parameters:
%
%     out          : Structure with the output defined by outstruct. The
%                    fields of out will have names as in outstruct. Each
%                    field will have size as given by the corresponding 
%                    vector in outstruct.
%     status       : Structure with the status and results, see system.
%
%   See also: verhulst2012 verhulst2018
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_extern.php


%   #Author: Piotr Majdak (2021): programmed for the AMT 1.0
%   #Author: Alejandro Osses (2023): compatible with AMT paths containing spaces (on Unix)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



switch environment
  case 'Python'
      % change to directory
    act_path=pwd;
    
    cd(fullfile(amt_basepath,'environments', directory));
      % input parameters
    if ~isempty(input)      
      save(fullfile('out','input.mat'), '-struct', 'input', '-v7');
    end
      % call
    dir_AMT_nospaces = amt_basepath;
    try
      if isunix % For Mac and Linux, spaces will be replaced by '\ ':          
        dir_AMT_nospaces = strrep(dir_AMT_nospaces, ' ' , '\ ' );
        dir_prefix = '';
        dir_suffix = '';
      end
      if ~isunix % i.e., is windows
        dir_prefix = '"'; % If there are spaces, the arguments will be between " ":
        dir_suffix = '"';
      end
    end
    argument2evaluate = fullfile(dir_AMT_nospaces,'environments', directory, module);
    argument2evaluate = [dir_prefix argument2evaluate dir_suffix];
    
    [stat,res]=system(['python ' argument2evaluate]);
    
    status.status=stat;
    status.res=res;
    if stat ~= 0
        amt_disp();
        amt_disp(res);
        error('AMT_EXTERN: Something went wrong calling Python (see message above)');
    end    
      % collect the output
    fn=fieldnames(outstruct);
    for ii=1:length(fn)
        varname=fn{ii};
        var=outstruct.(varname);
        dim1=var(1); 
        if length(var)<2, dim2=1; else dim2=var(2); end
        if length(var)<3, dim3=1; else dim3=var(3); end
        for jj=1:dim3
          filename=fullfile(amt_basepath,'environments', directory, 'out',[varname int2str(jj) '.np']);
          f=fopen(filename,'r');
          out.(varname)(:,:,jj)=fread(f,[dim1 dim2],'double','n');
          fclose(f);
          delete(filename);
        end
    end
      % clean up     
    if ~isempty(input), delete(fullfile('out','input.mat')); end
    cd(act_path);  
  otherwise
    error(['Environment ' environment ' is not supported']);
end


