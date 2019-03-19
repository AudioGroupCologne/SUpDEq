%demo_modelinitiative Demo for the interface to the model initiative using the Breebaart et al. (2001) model
%
%   demo_modelinitiative starts a server waiting for requests to process binaural
%   signals with the breebaart2001 model. See the BInit interface for more
%   details. 
%
%   demo_modelinitiative processes request placed in the temporary system directory.
%   Define the variable directory in order to use an other directory. 
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/demo_modelinitiative_breebaart2001.php

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

% PM 15.4.2017

if ~exist('directory','var'), directory=tempdir; end;


trials=1;
disp(['Waiting for requests in ' directory]);
delete(fullfile(directory,'interval_*.wav'));
while 1

	% wait for the signal to process the intervals
    while ~exist(fullfile(directory,'interval_1.wav'),'file');
        pause(.1);
    end

	% process all interval_*.wav files by the pathway model
	l=dir(fullfile(directory,'interval_*.wav'));
    pathway_out_1 = cell(length(l),1);
    pathway_out_2 = cell(length(l),1);
    pathway_out_3 = cell(length(l),1);
    pathway_out_4 = cell(length(l),1);
    % read parameters if provided
    par_path=[32000 0 0];
    if exist('a_priori.csv','file'), par_path = csvread('a_priori.csv'); end
    % load and process wave files
    for ii=length(l):-1:1
        % check if we can access the file
      fid=-1;
      while fid==-1
        fid=fopen(fullfile(directory,['interval_'  num2str(ii) '.wav']),'r');
      end
      fclose(fid);
        % load it
      [wave,fs] = audioread(fullfile(directory,['interval_'  num2str(ii) '.wav']));
        % delete it
      delete(fullfile(directory,['interval_'  num2str(ii) '.wav']));
        % call the model 
      disp(['Model started with interval #' num2str(ii)]);
      [pathway_out_1{ii}, pathway_out_2{ii}, pathway_out_3{ii}, pathway_out_4{ii}] = ...
         breebaart2001_preproc(wave,par_path(1),par_path(2),par_path(3));
    end

    % process detector model
    if exist(fullfile(directory,'decision_parameters.mat'),'file'), 
      par = load(fullfile(directory,'decision_parameters.mat')); 
      response = breebaart2001_centralproc(pathway_out_1, pathway_out_3, pathway_out_4 , par.decision.input4);
    else
      response = breebaart2001_centralproc(pathway_out_1, pathway_out_3, pathway_out_4 , 'lbr');
    end
    
    
    % write model response in a file.
    csvwrite(fullfile(directory,'choice.dat'),response);
	
    disp(['Decision: ' num2str(response) '. ' num2str(trials) ' trials handled.']);
    trials=trials+1;
end

