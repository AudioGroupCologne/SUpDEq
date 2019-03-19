function [out, par] = emuexp(command,par,varargin)
%emuexp Emulates psychoacoustic experiments
%   Usage:   par = emuexp(command,par);
%     [out,par] = emuexp('run',par);
%     [out, par] = emuexp('run',par,'plot');
%
%   Input parameters:
%     command  : One of the following commands. 'expinit' intializes the
%                general experiment parameters. 'signalinit' intializes
%                the signal generator creating model inputs. 'modelinit'
%                intializes model parameters. 'decisioninit' intializes
%                the parameters of the decision stage in the experiment.
%                Finally, 'run' runs the experiment and lets the model decide.
%     par      : Structure of the experimental parameters used by EMUEXP.
%                Set to [] on the first call (when par is not set up yet).
%
%   Output parameters:
%     par   : Structure containing all parameters
%     out   : Vector with the experiment output. out(:,1) is the average threshold of the
%             experimental variable. out(:,2) is the standard deviation of the variable
%             across all runs. out(:,3:end) provides the individual experimental variables
%             used in each trial.
%
%   par = EMUEXP(init_command,par) initializes the various parts of the
%   psychoacoustic experiment to be emulated depending on init_command.
%
%   out = EMUEXP('run',par) runs the experiment defined
%   by the structure par and outputs the experimental result in out.
%
%   [out, par] = EMUEXP('run',par) runs the experiment and outputs
%   more details on the experiment parameters in par.
%
%   [out, par] = EMUEXP('run',par,'plot') runs the experiment and
%   plots the experiment progress.
%
%   Initialization
%   --------------
%
%   Experiment
%   *********
%
%   par=EMUEXP('expinit',[],exp) initilizes the experiment wiht key-value pairs provided
%   in a cell array exp. The following pairs are required:
%
%     'intnum',intnum                number of intervals in a trial,
%                                    e.g. 3 sets up a 3-afc experiment.
%
%     'rule',down_up                 vector with down-up-rule
%                                    e.g. [2 1] sets up a 2-down, 1-up experiment.
%
%     'expvarstart',expvarstart      step size of the experimental variable at the
%                                    beginning of the experiment
%
%     'expvarsteprule',factor_turns  vector with a factor and number of turn arounds.
%                                    The factor affects the step size of the experimental
%                                    variable after the number of turn arounds, e.g. [0.5 2]
%                                    multiplies the stepsize by 0.5 after two turn arounds.
%
%     'stepmin',min_threshturn       vector with minimal step size and number of turn arounds
%                                    after reaching that minimal step size for the threshold
%                                    calculation. E.g. [1 8] means that after reaching the
%                                    step size 1, the experiment will continue for
%                                    another 8 reversals before terminating.
%
%   Signal generator
%   ***************
%
%   par=EMUEXP('signalinit',par,sig) intializes the signal generator creating
%   signals for the model with key-value pairs provided in the cell array sig. The signal
%   generator is called with those parameters in each trial of the experiment.
%   Up to 15 input parameters are supported. One of inputs must be 'inttyp': In each
%   experimental interval, this input will be replaced by 'target' or 'reference'
%   depending on the interval type. One of the inputs must be 'expvar': In each trial,
%   this input will be replaced by the value of the experimental variable. The
%   following pairs are required:
%
%     'name',name         string which defines the name of the signal
%                         generation
%
%     'inputX',inputX    input parameter X needed for the signal generator
%
%   Model called in each interval
%   ****************************
%
%   par=EMUEXP('modelinit',par,mod) initializes the model called in each interval with
%   the key-value pairs provided in mod. Up to 10 input parameters are supported.
%   One of the inputs must contain the keyword 'expsignal'.
%   This keyword is replaced in the 'run' routine with the output of
%   the signal generation function:
%
%     'name',name         string which defines the name of the model
%                         function
%
%     'inputX',intputX    input parameter X needed by the model
%
%     'outputs',outputs   indicies of used model outputs for the decision
%                         e.g. [1 2 6]: output 1,2 and 6 used
%
%   Decision stage called in each trial
%   **********************************
%
%   par=EMUEXP('decisioninit',par,dec) initializes the decision stage of the experiment
%   with key-value pairs provided in dec. Up to 10 input parameters are supported.
%   All inputs containing the keyword 'modelout' are
%   replaced with the outputs of the model function during an
%   experimental run. Therefore the number of inputs with the keyword
%   'modelout' must be equal to number of 'outputs' defined in
%   'modelinit'. An output of the modelfunction contains a cell with an
%   entry for each interval. E.g. param1{1} contains the first output of
%   the model function of the first interval and param3{2} contains the
%   third output of the modelfunction of the second interval. Therefore the
%   decision function must be implemented so that the inputs of the decision
%   function are cells with entries for each interval.
%   Following parameters are required:
%
%     'name',name         name of the decision fuction
%     'inputX',intputX    input parameter X needed by the decision function
%
%   Running the experiment
%   ----------------------
%
%   After the initialization, the experiment can be started by
%   out = EMUEXP('run',par);. The threshold will be in out.
%
%   See also: exp_breebaart2001, demo_breebaart2001
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/emuexp.php

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


% AUTHOR: Martina Kreuzbichler

% turn warning off
warning('off','MATLAB:nargchk:deprecated')

switch command
    case 'expinit'
        definput.keyvals.intnum = [];
        definput.keyvals.rule = [];
        definput.keyvals.expvarstepstart = [];
        definput.keyvals.expvarsteprule = [];
        definput.keyvals.stepmin = [];
        definput.keyvals.expvarstart = [];
        definput.keyvals.interface = 'AMT';
        definput.keyvals.fs = [];
        definput.keyvals.directory = [];

        if iscell(varargin{1}) && nargin == 3
            [~,kvexp]=ltfatarghelper({},definput,varargin{:});
        else
            [~,kvexp]=ltfatarghelper({},definput,varargin);
        end

       % out = setstructfields(kvexp, par);
       out = par;
       out.exp = kvexp;

    case 'modelinit'
        definput.keyvals.name=[];
        definput.keyvals.input1 = [];
        definput.keyvals.input2 = [];
        definput.keyvals.input3 = [];
        definput.keyvals.input4 = [];
        definput.keyvals.input5 = [];
        definput.keyvals.input6 = [];
        definput.keyvals.input7 = [];
        definput.keyvals.input8 = [];
        definput.keyvals.input9 = [];
        definput.keyvals.input10 = [];
        definput.keyvals.outputs = [];

        if iscell(varargin{1}) && nargin == 3
            [~,kvmodel]=ltfatarghelper({},definput,varargin{:});
        else
            [~,kvmodel]=ltfatarghelper({},definput,varargin);
        end

        % check what outputs of model are needed
        outputnumber = nargout(kvmodel.name);

        callmodelstring = [];

        for outputcounter = 1:outputnumber
            if any(kvmodel.outputs==outputcounter)
                modelstring = sprintf('modelout.par%i',outputcounter);
                if isempty(callmodelstring)
                    callmodelstring = ['[' modelstring '{interval_num}'];
                else
                    callmodelstring =  [callmodelstring ',' modelstring '{interval_num}'];
                end
            else
                callmodelstring = [callmodelstring ',~'];
            end
        end

        modelinputs = struct2cell(kvmodel);

        % delete name of model function & outputs
        modelinputs = modelinputs(2:end-1);

        % delete empty cells
        modelinputs = modelinputs(~cellfun('isempty',modelinputs));

        callmodelstring = [callmodelstring ']=' kvmodel.name '('];

        for modelinputscounter = 1:length(modelinputs)
            callmodelstring = [callmodelstring num2str(modelinputs{modelinputscounter}) ','];
        end
        callmodelstring(end:end+1) = ');';

        out = par;
        out.model = kvmodel;
        out.callstrings.model = callmodelstring;

        %TODO CALLMODELSTRING

    case 'signalinit'
        definput.keyvals.name= [];
        definput.keyvals.input1 = [];
        definput.keyvals.input2 = [];
        definput.keyvals.input3 = [];
        definput.keyvals.input4 = [];
        definput.keyvals.input5 = [];
        definput.keyvals.input6 = [];
        definput.keyvals.input7 = [];
        definput.keyvals.input8 = [];
        definput.keyvals.input9 = [];
        definput.keyvals.input10 = [];
        definput.keyvals.input11 = [];
        definput.keyvals.input12 = [];
        definput.keyvals.input13 = [];
        definput.keyvals.input14 = [];
        definput.keyvals.input15 = [];


        if iscell(varargin{1}) && nargin == 3
            [~,kvsignal]=ltfatarghelper({},definput,varargin{:});
        else
            [~,kvsignal]=ltfatarghelper({},definput,varargin);
        end

        out = par;
        out.signal = kvsignal;

    case 'decisioninit'
        definput.keyvals.name= [];
        definput.keyvals.input1 = [];
        definput.keyvals.input2 = [];
        definput.keyvals.input3 = [];
        definput.keyvals.input4 = [];
        definput.keyvals.input5 = [];
        definput.keyvals.input6 = [];
        definput.keyvals.input7 = [];
        definput.keyvals.input8 = [];
        definput.keyvals.input9 = [];
        definput.keyvals.input10 = [];

        definput.commands.plot = {'noplot','plot'};

        if iscell(varargin{1}) && nargin == 3
            [~,kvdecision]=ltfatarghelper({},definput,varargin{:});
        else
            [~,kvdecision]=ltfatarghelper({},definput,varargin);
        end

        out = par;
        out.decision = kvdecision;

    case 'run'

        definput.flags.plot = {'noplot','plot'};

        [commands,~]=ltfatarghelper({},definput,varargin);
        if strcmp(par.exp.interface,'ModelInitiative'),          
          csvwrite(fullfile(par.exp.directory,'a_priori.csv'),...
            [par.model.input2, par.model.input3, par.model.input4]);
          decision=par.decision;
          save(fullfile(par.exp.directory,'decision_parameters.mat'),'decision');
          delete(fullfile(par.exp.directory,'detector_out.csv'));
        end

        % find experimental variable and inttyp variable
        sigparnames = fieldnames(par.signal);
        for count = 1: length(sigparnames)
            name = sigparnames(count);
            if strcmp(getfield(par.signal, name{:}),'expvar')
                experimentvar = name{:};
            elseif strcmp(getfield(par.signal, name{:}),'inttyp')
                inttypvar = name{:};
            end
        end

        stepsize = par.exp.expvarstepstart;
        truecounter = 0;
        par.signal.(experimentvar) = par.exp.expvarstart;
        expparvalue = [];
        downturn = 0;
        upturn = 1;
        turncounter = 0;
        lastturn = [];
        checkmodelout = 0;
        condition = 1;
        wrongcounter = -1;
        trialcounter = 1;

        while condition

            for interval_num=1:par.exp.intnum

                if interval_num == 1
                    par.signal.(inttypvar) = 'target';
                else
                    par.signal.(inttypvar) = 'reference';
                end

                signalinputs = struct2cell(par.signal);

                % delete name of signal function
                signalinputs = signalinputs(2:end);

                % delete empty cells
                signalinputs = signalinputs(~cellfun('isempty',signalinputs));

                % call signalfunction
                testsignal = feval(par.signal.name,signalinputs{:});

                % call model or transfer files to the model
                switch par.exp.interface
                  case 'AMT'
                    % find experimental signal variable
                    modelparnames = fieldnames(par.model);
                    for count = 1: length(modelparnames)
                        name = modelparnames(count);
                        if strcmp(getfield(par.model, name{:}),'expsignal')
                            experimentsignal = name{:};
                            break
                        end
                    end
                    par.model.(experimentsignal) = testsignal;
                    % call model
                    par.callstrings.model = strrep(par.callstrings.model,'expsignal', 'testsignal');
                    eval(par.callstrings.model);
                  case 'ModelInitiative'
                    audiowrite(fullfile(par.exp.directory,['interval_' num2str(interval_num) '.wav']), testsignal, par.exp.fs);
                end
            end % for each interval

            switch par.exp.interface
              case 'AMT'
                % find model output variable, only at the first time
                if checkmodelout == 0
                    modeloutvarcount = 1;
                    decisionparnames = fieldnames(par.decision);
                    for count = 1:length(decisionparnames)
                        name = decisionparnames(count);
                        if strcmp(getfield(par.decision, name{:}),'modelout')
                            modeloutvar{modeloutvarcount} = name{:};
                            modeloutvarcount = modeloutvarcount+1;
                        end
                    end
                    checkmodelout = 1;
                end

                %set model outputs to decsion inputs
                for count = 1:modeloutvarcount-1
                    modeloutname = sprintf('par%i',par.model.outputs(count));
                    par.decision.(modeloutvar{count}) = modelout.(modeloutname);
                end

                decisioninputs = struct2cell(par.decision);

                % delete name of model function & outputs
                decisioninputs = decisioninputs(2:end);

                % delete empty cells
                decisioninputs = decisioninputs(~cellfun('isempty',decisioninputs));

                % call decision
                decision = feval(par.decision.name,decisioninputs{:});
              case 'ModelInitiative'
                amt_disp(['Trial #' num2str(trialcounter) ', ModelInitiative on ' par.exp.directory],'volatile');
                while ~exist(fullfile(par.exp.directory,'detector_out.csv'),'file');
                  pause(.1);
                end
                fid=-1;
                while fid==-1
                  fid=fopen(fullfile(par.exp.directory,'detector_out.csv'),'r');
                end
                fclose(fid);                
                decision=csvread(fullfile(par.exp.directory,'detector_out.csv'));
                delete(fullfile(par.exp.directory,'detector_out.csv'));
            end

            % store expparvalue
            expparvalue = [expparvalue par.signal.(experimentvar)];

            % count reversals
            % wrong answers are par.exp.rule(2)and no low point reversal
            if decision ~= 1 && wrongcounter == par.exp.rule(2)-1 && downturn == 0
                turncounter = turncounter + 1;
                downturn = 1;
                upturn = 0;
                if stepsize == par.exp.stepmin(1)
                    lastturn = [lastturn par.signal.(experimentvar)];
                end

             % right answers are par.exp.rule(1)
             % and no high point reversal
            elseif decision == 1 && truecounter == par.exp.rule(1)-1 && upturn == 0
                turncounter = turncounter + 1;
                downturn = 0;
                upturn = 1;
                if stepsize == par.exp.stepmin(1)
                    lastturn = [lastturn par.signal.(experimentvar)];
                end
            end

            % change stepsize after par.exp.expvarsteprule(2) reversals, if stepsize
            % is not already par.exp.stepmin(1) dB
            if turncounter == par.exp.expvarsteprule(2) && truecounter == 1 && ...
                stepsize ~= par.exp.stepmin(1)
                stepsize = stepsize * par.exp.expvarsteprule(1);
                turncounter = 0;
            end


            if decision == 1
                wrongcounter = 0;
                truecounter = truecounter + 1;
            else
                truecounter = 0;
                wrongcounter = wrongcounter +1;
            end

            % amount of right answers is par.exp.rule(1)
            if truecounter == par.exp.rule(1)
                par.signal.(experimentvar) = par.signal.(experimentvar) - stepsize;
                truecounter = 0;
            end

            % amount of wrong answers is par.exp.rule(2)
            if wrongcounter == par.exp.rule(2)
                par.signal.(experimentvar) = par.signal.(experimentvar) + stepsize;
                wrongcounter = 0;
            end

            % amount of reversals is par.stepmin(2)
            if size(lastturn,2) == par.exp.stepmin(2)
                threshold = median(lastturn);
                threshstd = std(lastturn,1);
                condition = 0;
                out = [threshold threshstd expparvalue];
            end
            
            trialcounter = trialcounter+1;

        end
        %clear persistent variables
        clear (par.decision.name);

        if commands.do_plot
            figure
            for plotcounter = 1:(size(expparvalue,2)-1)
                if expparvalue(plotcounter) < expparvalue(plotcounter+1)
                    style = 'or';
                    stem(plotcounter,(expparvalue(plotcounter)),'or',...
                        'LineStyle','none')
                    hold on
                else
                    style = '+g';
                    stem(plotcounter,(expparvalue(plotcounter)),'+g',...
                        'LineStyle','none')
                    hold on
                end
            end

            % special case last entry
            if decision == 1
                stem(plotcounter+1,(expparvalue(plotcounter+1)),'+g',...
                        'LineStyle','none')
            else
                stem(plotcounter+1,(expparvalue(plotcounter+1)),'or',...
                        'LineStyle','none')
            end
            title(['Threshold (Median): ' num2str(threshold) ...
                'dB, Std: ' num2str(threshstd) 'dB'])

        end

end

