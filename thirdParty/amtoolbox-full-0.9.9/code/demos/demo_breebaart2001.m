function demo_breebaart2001(varargin)
%DEMO_BREEBAART2001 Demo for binaural processingmodel from Breebaart et al. (2001)
%
%   DEMO_BREEBAART2001(flags) demonstrates how to compute N0S0, N0Spi and
%   NpiS0 thresholds for different combinations of the monaural and
%   binaural central processor for a center frequency of 1000 or 4000 Hz 
%   using the model from Breebaart et al. (2014).
%
%   
%   The N0Spi condition is default. The following combinations for 4000 Hz 
%   are calculated:
%
%     lbr  Compute decision from combination of left, binaural and
%          right result. If the binaural result is zero, don't use it. 
%          In this condition lbr is equal to lBr, because the binaural
%          result is not zero.
%
%     b    Compute decision from binaural result only.
%
%     lr   Compute decision from left and right result.
%
%   
%   Set flag to the following flags to shows other conditions:
%
%     'N0S0'  N0S0 thresholds for differnt combinations of the monaural and
%             binaural processor for a center frequency of 4000 Hz.
%             The following combinations are calculated:
%
%              lbr  Compute decision from combination of left, binaural and
%                    right result. If the binaural result is zero, don't use
%                    it. In this condition lbr is equal to lr.
%
%              lBr  Compute decision from combination of left, binaural and
%                    right result even if the binaural result is zero. 
%           
%     'NpiS0' NpiS0 thresholds for differnt combinations of the monaural and
%             binaural processor for a center frequency of 1000 Hz.
%             The following combinations are calculated:
%
%             'lbr'  Compute decision from combination of left, binaural and
%                    right result. If the binaural result is zero, don't use it. 
%                    In this condition lbr is equal to lBr, because the binaural
%                    result is not zero.
%
%             'b'    Compute decision from binaural result only.
%
%             'lr'   Compute decision from left and right result.
%
%     'exact' The threshold is computed as a mean value of 6 repetitions per
%             bandwidth. To get results faster the computational results are
%             stored in cache and can be downloaded.
%
%     'fast' The threshold as the result of one experimental run per bandwidth.
%            The number of bandwidth values is half the number of values at the
%            'exact' computation.  The fast computation is default.
%
% 
%   demo_breebaart2001('N0Spi') shows :
%
%     demo_breebaart2001('N0Spi');
%
%   N0Spi thresholds  as a function of the masker bandwidth for a center
%   frequency of 4000 Hz and lbr, b and lr conditions in comparison with
%   the model results of Breebaart et al. (2001b).
%
%
%   demo_breebaart2001('N0S0') shows :
%
%     demo_breebaart2001('N0S0');
%
%   N0S0 thresholds  as a function of the masker bandwidth for a center
%   frequency of 4000 Hz and lbr and lBr conditions in comparison with
%   experimental results of van de Par and Kohlrausch (1999).
%
%
%   demo_breebaart2001('NpiS0') shows :
%
%     demo_breebaart2001('NpiS0');
%
%   NpiS0 thresholds  as a function of the masker bandwidth for a center
%   frequency of 4000 Hz and lbr, b and lr conditions in comparison with
%   the model results of Breebaart et al. (2001b).
%
%   See also: breebaart2001_centralproc breebaart2001_preproc
%   sig_breebaart2001
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. II. Dependence on spectral
%     parameters. J. Acoust. Soc. Am., 110:1089-1104, August 2001.
%     
%     S. van de Par and A. Kohlrausch. Dependence of binaural masking level
%     differences on center frequency, masker bandwidth, and interaural
%     parameters. J. Acoust. Soc. Am., 106(4):1940-1947, 1999.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_breebaart2001.php

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

% AUTHOR : Martina Kreuzbichler

definput.import={'amt_cache'};
definput.flags.type = {'N0Spi','N0S0','NpiS0'};
definput.flags.speed = {'fast','exact'};
definput.flags.interface = {'AMT', 'BInit'};
definput.keyvals.directory=tempdir;

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);


if flags.do_fast
    
    if flags.do_N0Spi
        % set experimental paramters
        parout = [];
        switch flags.interface
          case 'AMT'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',65};
          case 'BInit'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',65, ...
                'interface','BInit','directory',kv.directory,'fs',32000};
        end
        parout = emuexp('expinit',parout,expset);

        % set model parameters
        % input2 = fs; input3 = tau; input4 = ild; 
        modelset = {'name','breebaart2001_preproc','input1',...
            'expsignal', 'input2',32000,'input3',0,'input4',0,...
            'outputs',[1 3 4]};
        parout = emuexp('modelinit',parout,modelset);

        % run model
        centerfreq = 4000;
        bw = [5 25 100 500 2000];
        centralprocstring = {'lbr','b','lr'};
        nl = 65;

    amt_disp('demo_breebaart2001 will calculate 15 thresholds','progress');
		output = amt_cache('get','N0Spi_fast',flags.cachemode);		
		if isempty(output)
      runcounter=1;
			% loop for all centralprocstring conditions.
			for stringcount = 1:length(centralprocstring)

				% set decision parameters
				decisionset = {'name','breebaart2001_centralproc',...
					'input1','modelout','input2','modelout','input3',...
					'modelout','input4',centralprocstring{stringcount}};
				parout = emuexp('decisioninit',parout,decisionset);
				
				resultbwvec = zeros(1,1);
				resultvec = zeros(5,1);
				
				%loop for all bandwidths
				for bwcount = 1:length(bw)

					% set signal parameters
					%input3 = signallevel; input4 = signalduration;
					%input5 = signalphase; input7 = noiselevel;
					%input8 = noiseduration; input9 = noisephase;
					%input10 = hanning ramp duration; input11 = fs;
					signalset = {'name','sig_breebaart2001','input1',...
						'inttyp', 'input2',centerfreq,'input3',...
						'expvar','input4',0.15,'input5',pi,'input6',...
						bw(bwcount),'input7',nl,'input8',0.2,'input9',0,...
						'input10',0.025,'input11', 32000};
					parout = emuexp('signalinit',parout,signalset);

          amt_disp(['Threshold #' num2str(runcounter) ': for bw=' num2str(bw(bwcount)) 'Hz and decision=' centralprocstring{stringcount} '.'],'progress');
					result = emuexp('run',parout);
					resultbwvec(1) = result(1)-nl;

					resultvec(bwcount) = resultbwvec;
					resultbwvec = zeros(1,1);
					amt_disp(sprintf(['Progress for central processor condition ' ...
						centralprocstring{stringcount} ':' num2str(round(bwcount/length(bw)*100)) ...
						'%% calculated']),'progress');
          runcounter=runcounter+1;
					result_temp =sprintf(['N0Spi4000' centralprocstring{stringcount}]);
					output.(result_temp) = resultvec;

				end
			end
			amt_cache('set','N0Spi_fast',output);
		end

        N0Spi4000lbr = output.N0Spi4000lbr;
        N0Spi4000b = output.N0Spi4000b;
        N0Spi4000lr = output.N0Spi4000lr;
        
        % plot

        %get model data
        N0Spimodeldata4000 = data_breebaart2001('fig3','nfc4000');

        bw_model = [5 10 25 50 100 250 500 1000 2000 4000];
        bw_demo = [5 25 100 500 2000];

        figure
        plot(bw_demo,N0Spi4000lbr,'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 20])
        hold on;
        plot(bw_demo,N0Spi4000b,'-o','MarkerSize',10,'Color',[0,0.7,1]);
        plot(bw_demo,N0Spi4000lr,'-dc','MarkerSize',10);
        plot(bw_model,N0Spimodeldata4000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        legend('N0Spi condition lbr','N0Spi condition b',...
            'N0Spi condition lr','model data from Breebaart (2001)')
        title('N0Spi 4000 Hz thresholds: monaural factor = 0.0003 - fast computation')
    
    elseif flags.do_N0S0
        % set experimental paramters
        parout = [];
        expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
            'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',90};
        parout = emuexp('expinit',parout,expset);

        % set model parameters
        % input2 = fs; input3 = tau; input4 = ild; 
        modelset = {'name','breebaart2001_preproc','input1',...
            'expsignal', 'input2',32000,'input3',0,'input4',0,...
            'outputs',[1 3 4]};
        parout = emuexp('modelinit',parout,modelset);

        % run model
        centerfreq = 4000;
        bw = [5 25 100 500 2000];
        centralprocstring = {'lbr', 'lBr'};
        nl = 70;

		output = amt_cache('get','N0S0_fast',flags.cachemode);		
		if isempty(output)
			% loop for all centralprocstring conditions.
			for stringcount = 1:length(centralprocstring)

				% set decision parameters
				decisionset = {'name','breebaart2001_centralproc',...
					'input1','modelout','input2','modelout','input3',...
					'modelout','input4',centralprocstring{stringcount}};
				parout = emuexp('decisioninit',parout,decisionset);
				
				resultbwvec = zeros(1,1);
				resultvec = zeros(5,1);

				%loop for all bandwidths
				for bwcount = 1:length(bw)

					% set signal parameters
					%input3 = signallevel; input4 = signalduration;
					%input5 = signalphase; input7 = noiselevel;
					%input8 = noiseduration; input9 = noisephase;
					%input10 = hanning ramp duration; input11 = fs;
					signalset = {'name','sig_breebaart2001','input1',...
						'inttyp', 'input2',centerfreq,'input3',...
						'expvar','input4',0.15,'input5',0,'input6',...
						bw(bwcount),'input7',nl,'input8',0.2,'input9',0,...
						'input10',0.025,'input11', 32000};
					parout = emuexp('signalinit',parout,signalset);

					result = emuexp('run',parout);
					resultbwvec(1) = result(1)-nl;

					resultvec(bwcount) = resultbwvec;
					resultbwvec = zeros(1,1);


					amt_disp(sprintf(['Progress for central processor condition ' ...
						centralprocstring{stringcount} ': ' ...
						num2str(round(bwcount/length(bw)*100)) ...
						'%% calculated']),'progress');

					result_temp =sprintf(['N0S04000' centralprocstring{stringcount}]);
					output.(result_temp) = resultvec;

				end
			end
			amt_cache('set','N0S0_fast',output);
		end
        N0S04000lbr = output.N0S04000lbr;
        N0S04000lBr = output.N0S04000lBr;
        
        %plot
        N0S0expdata4000 = data_vandepar1999('fig1_N0S0','nfc4000');

        bw_model = [5 10 25 50 100 250 500 1000 2000 4000];
        bw_demo = [5 25 100 500 2000];

        figure
        plot(bw_demo,N0S04000lbr,'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -15 10])
        hold on;
        plot(bw_demo,N0S04000lBr,'-o','MarkerSize',10,'Color',[0,0.7,1]);
        plot(bw_model,N0S0expdata4000,'-sg','MarkerSize',10,'MarkerFaceColor','g');
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        legend('N0S0 condition lbr','N0S0 condition lBr','experimental data from van de Par and Kohlrausch (1999)')
        title('N0S0 4000 Hz thresholds: monaural factor = 0.0003 - fast computation')
        
        
    elseif flags.do_NpiS0
        % set experimental paramters
        parout = [];
        expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
            'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',85};
        parout = emuexp('expinit',parout,expset);

        % set model parameters
        % input2 = fs; input3 = tau; input4 = ild; 
        modelset = {'name','breebaart2001_preproc','input1',...
            'expsignal', 'input2',32000,'input3',0.0005,'input4',0,...
            'outputs',[1 3 4]};
        parout = emuexp('modelinit',parout,modelset);

        % run model
        centerfreq = 1000;
        bw = [5 25 100 500 2000];
        centralprocstring = {'lbr','b','lr'};
        nl = 70;

		output = amt_cache('get','NpiS0_fast',flags.cachemode);		
		if isempty(output)
			% loop for all centralprocstring conditions.
			for stringcount = 1:length(centralprocstring)

				% set decision parameters
				decisionset = {'name','breebaart2001_centralproc',...
					'input1','modelout','input2','modelout','input3',...
					'modelout','input4',centralprocstring{stringcount}};
				parout = emuexp('decisioninit',parout,decisionset);
				
				resultbwvec = zeros(1,1);
				resultvec = zeros(5,1);

				%loop for all bandwidths
				for bwcount = 1:length(bw)

					% set signal parameters
					%input3 = signallevel; input4 = signalduration;
					%input5 = signalphase; input7 = noiselevel;
					%input8 = noiseduration; input9 = noisephase;
					%input10 = hanning ramp duration; input11 = fs;
					signalset = {'name','sig_breebaart2001','input1',...
						'inttyp', 'input2',centerfreq,'input3',...
						'expvar','input4',0.15,'input5',0,'input6',...
						bw(bwcount),'input7',nl,'input8',0.2,'input9',pi,...
						'input10',0.025,'input11', 32000};
					parout = emuexp('signalinit',parout,signalset);
					
					result = emuexp('run',parout);
					resultbwvec(1) = result(1)-nl;

					resultvec(bwcount) = resultbwvec;
					resultbwvec = zeros(1,1);


					amt_disp(sprintf(['Progress for central processor condition ' ...
						centralprocstring{stringcount} ': ' ...
						num2str(round(bwcount/length(bw)*100)) ...
						'%% calculated']),'progress');

					result_temp =sprintf(['NpiS01000' centralprocstring{stringcount}]);
					output.(result_temp) = resultvec;

				end
			end
			amt_cache('set','NpiS0_fast',output);
		end

        NpiS01000lbr = output.NpiS01000lbr;
        NpiS01000b = output.NpiS01000b;
        NpiS01000lr = output.NpiS01000lr;

        % plot

        %get model data
        NpiS0modeldata1000 = data_breebaart2001('fig6','nfc1000');

        bw_model = [5 10 25 50 100 250 500 1000 2000];
        bw_demo = [5 25 100 500 2000];

        figure
        plot(bw_demo,NpiS01000lbr,'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -30 10])
        hold on;
        plot(bw_demo,NpiS01000b,'-o','MarkerSize',10,'Color',[0,0.7,1]);
        plot(bw_demo,NpiS01000lr,'-dc','MarkerSize',10);
        plot(bw_model,NpiS0modeldata1000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        legend({'NpiS0 condition lbr','NpiS0 condition b',...
            'NpiS0 condition lr','model data from Breebaart (2001)'},...
            [290,170,0.1,0.05]);
        title('NpiS0 1000 Hz thresholds: monaural factor = 0.0003 - fast computation')
    end


elseif flags.do_exact

    if flags.do_N0Spi
        [N0Spi4000lbr,N0Spi4000b,N0Spi4000lr] = amt_cache('get','N0Spi',flags.cachemode);
        if isempty(N0Spi4000lbr)
            % do computation

            % set experimental paramters
            parout = [];
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',65};
            parout = emuexp('expinit',parout,expset);

            % set model parameters
            % input2 = fs; input3 = tau; input4 = ild; 
            modelset = {'name','breebaart2001_preproc','input1',...
                'expsignal', 'input2',32000,'input3',0,'input4',0,...
                'outputs',[1 3 4]};
            parout = emuexp('modelinit',parout,modelset);

            % run model
            centerfreq = 4000;
            bw = [5 10 25 50 100 250 500 1000 2000 4000];
            centralprocstring = {'lbr','b','lr'};
            nl = 65;

            % loop for all centralprocstring conditions.
            for stringcount = 1:length(centralprocstring)

                % set decision parameters
                decisionset = {'name','breebaart2001_centralproc',...
                    'input1','modelout','input2','modelout','input3',...
                    'modelout','input4',centralprocstring{stringcount}};
                parout = emuexp('decisioninit',parout,decisionset);

                %loop for all bandwidths
                for bwcount = 1:length(bw)

                    % set signal parameters
                    %input3 = signallevel; input4 = signalduration;
                    %input5 = signalphase; input7 = noiselevel;
                    %input8 = noiseduration; input9 = noisephase;
                    %input10 = hanning ramp duration; input11 = fs;
                    signalset = {'name','sig_breebaart2001','input1',...
                        'inttyp', 'input2',centerfreq,'input3',...
                        'expvar','input4',0.3,'input5',pi,'input6',...
                        bw(bwcount),'input7',nl,'input8',0.4,'input9',0,...
                        'input10',0.05,'input11', 32000};
                    parout = emuexp('signalinit',parout,signalset);

                    % loop for experimental runs
                    for runcounter = 1:6
                        result = emuexp('run',parout);
                        resultbwvec(runcounter) = result(1)-nl;
                    end

                    resultvec(bwcount) = mean(resultbwvec);
                    resultvecstd(bwcount) = std(resultbwvec,1);
                    resultbwvec = zeros(6,1);
                    amt_disp(sprintf(['Progress for central processor condition ' ...
                        centralprocstring{stringcount} ':' num2str(round(bwcount/length(bw)*100)) ...
                        '%% calculated']),'progress');

                    result_temp =sprintf(['N0Spi4000' centralprocstring{stringcount}]);
                    output.(result_temp) = [resultvec; resultvecstd]';

                end
            end

            N0Spi4000lbr = output.N0Spi4000lbr;
            N0Spi4000b = output.N0Spi4000b;
            N0Spi4000lr = output.N0Spi4000lr;

            amt_cache('set','N0Spi',N0Spi4000lbr,N0Spi4000b,N0Spi4000lr);
        end

        % plot

        %get model data
        N0Spimodeldata4000 = data_breebaart2001('fig3','nfc4000');

        bw = [5 10 25 50 100 250 500 1000 2000 4000];

        figure
        errorbar(bw,N0Spi4000lbr(:,1),N0Spi4000lbr(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 20])
        hold on;
        errorbar(bw,N0Spi4000b(:,1),N0Spi4000b(:,2),'-o','MarkerSize',10,'Color',[0,0.7,1]);
        errorbar(bw,N0Spi4000lr(:,1),N0Spi4000lr(:,2),'-dc','MarkerSize',10);
        plot(bw,N0Spimodeldata4000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        legend('N0Spi condition lbr','N0Spi condition b',...
            'N0Spi condition lr','model data from Breebaart (2001)')
        title('N0Spi 4000 Hz thresholds: monaural factor = 0.0003')

    elseif flags.do_N0S0
        [N0S04000lbr,N0S04000lBr] = amt_cache('get','N0S0',flags.cachemode);
        if isempty(N0S04000lbr)
            % do computation

            % set experimental paramters
            parout = [];
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',90};
            parout = emuexp('expinit',parout,expset);

            % set model parameters
            % input2 = fs; input3 = tau; input4 = ild; 
            modelset = {'name','breebaart2001_preproc','input1',...
                'expsignal', 'input2',32000,'input3',0,'input4',0,...
                'outputs',[1 3 4]};
            parout = emuexp('modelinit',parout,modelset);

            % run model
            centerfreq = 4000;
            bw = [5 10 25 50 100 250 500 1000 2000 4000];
            centralprocstring = {'lbr', 'lBr'};
            nl = 70;

            % loop for all centralprocstring conditions.
            for stringcount = 1:length(centralprocstring)

                % set decision parameters
                decisionset = {'name','breebaart2001_centralproc',...
                    'input1','modelout','input2','modelout','input3',...
                    'modelout','input4',centralprocstring{stringcount}};
                parout = emuexp('decisioninit',parout,decisionset);

                %loop for all bandwidths
                for bwcount = 1:length(bw)

                    % set signal parameters
                    %input3 = signallevel; input4 = signalduration;
                    %input5 = signalphase; input7 = noiselevel;
                    %input8 = noiseduration; input9 = noisephase;
                    %input10 = hanning ramp duration; input11 = fs;
                    signalset = {'name','sig_breebaart2001','input1',...
                        'inttyp', 'input2',centerfreq,'input3',...
                        'expvar','input4',0.3,'input5',0,'input6',...
                        bw(bwcount),'input7',nl,'input8',0.4,'input9',0,...
                        'input10',0.05,'input11', 32000};
                    parout = emuexp('signalinit',parout,signalset);

                    % loop for experimental runs
                    for runcounter = 1:6
                        result = emuexp('run',parout);
                        resultbwvec(runcounter) = result(1)-nl;
                    end

                    resultvec(bwcount) = mean(resultbwvec);
                    resultvecstd(bwcount) = std(resultbwvec,1);
                    resultbwvec = zeros(6,1);
                    amt_disp(sprintf(['Progress for central processor condition ' ...
                        centralprocstring{stringcount} ': ' ...
                        num2str(round(bwcount/length(bw)*100)) ...
                        '%% calculated']),'progress');

                    result_temp =sprintf(['N0S04000' centralprocstring{stringcount}]);
                    output.(result_temp) = [resultvec; resultvecstd]';

                end
            end
            N0S04000lbr = output.N0S04000lbr;
            N0S04000lBr = output.N0S04000lBr;


            amt_cache('set','N0S0',N0S04000lbr,N0S04000lBr);
        end

        % plot

        % get experimental data
        N0S0expdata4000 = data_vandepar1999('fig1_N0S0','nfc4000');

        bw = [5 10 25 50 100 250 500 1000 2000 4000];

        figure
        errorbar(bw,N0S04000lbr(:,1),N0S04000lbr(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -15 10])
        hold on;
        errorbar(bw,N0S04000lBr(:,1),N0S04000lBr(:,2),'-o','MarkerSize',10,'Color',[0,0.7,1]);
        plot(bw,N0S0expdata4000,'-sg','MarkerSize',10,'MarkerFaceColor','g');
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        legend('N0S0 condition lbr','N0S0 condition lBr','experimental data from van de Par and Kohlrausch (1999)')
        title('N0S0 4000 Hz thresholds: monaural factor = 0.0003')


    elseif flags.do_NpiS0
        [NpiS01000lbr,NpiS01000b,NpiS01000lr] = amt_cache('get','NpiS0',flags.cachemode);
        if isempty(NpiS01000lbr)
            % do computation

            % set experimental paramters
            parout = [];
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',85};
            parout = emuexp('expinit',parout,expset);

            % set model parameters
            % input2 = fs; input3 = tau; input4 = ild; 
            modelset = {'name','breebaart2001_preproc','input1',...
                'expsignal', 'input2',32000,'input3',0.0005,'input4',0,...
                'outputs',[1 3 4]};
            parout = emuexp('modelinit',parout,modelset);

            % run model
            centerfreq = 1000;
            bw = [5 10 25 50 100 250 500 1000 2000];
            centralprocstring = {'lbr','b','lr'};
            nl = 70;

            % loop for all centralprocstring conditions.
            for stringcount = 1:length(centralprocstring)

                % set decision parameters
                decisionset = {'name','breebaart2001_centralproc',...
                    'input1','modelout','input2','modelout','input3',...
                    'modelout','input4',centralprocstring{stringcount}};
                parout = emuexp('decisioninit',parout,decisionset);

                %loop for all bandwidths
                for bwcount = 1:length(bw)

                    % set signal parameters
                    %input3 = signallevel; input4 = signalduration;
                    %input5 = signalphase; input7 = noiselevel;
                    %input8 = noiseduration; input9 = noisephase;
                    %input10 = hanning ramp duration; input11 = fs;
                    signalset = {'name','sig_breebaart2001','input1',...
                        'inttyp', 'input2',centerfreq,'input3',...
                        'expvar','input4',0.3,'input5',0,'input6',...
                        bw(bwcount),'input7',nl,'input8',0.4,'input9',pi,...
                        'input10',0.05,'input11', 32000};
                    parout = emuexp('signalinit',parout,signalset);

                    % loop for experimental runs
                    for runcounter = 1:6
                        result = emuexp('run',parout);
                        resultbwvec(runcounter) = result(1)-nl;
                    end

                    resultvec(bwcount) = mean(resultbwvec);
                    resultvecstd(bwcount) = std(resultbwvec,1);
                    resultbwvec = zeros(6,1);
                    amt_disp(sprintf(['Progress for central processor condition ' ...
                        centralprocstring{stringcount} ':' num2str(round(bwcount/length(bw)*100)) ...
                        '%% calculated']),'progress');

                    result_temp =sprintf(['NpiS01000' centralprocstring{stringcount}]);
                    output.(result_temp) = [resultvec; resultvecstd]';

                end
            end

            NpiS01000lbr = output.NpiS01000lbr;
            NpiS01000b = output.NpiS01000b;
            NpiS01000lr = output.NpiS01000lr;

            amt_cache('set','NpiS0',NpiS01000lbr,NpiS01000b,NpiS01000lr);
        end

        % plot

        %get model data
        NpiS0modeldata1000 = data_breebaart2001('fig6','nfc1000');

        bw = [5 10 25 50 100 250 500 1000 2000];

        figure
        errorbar(bw,NpiS01000lbr(:,1),NpiS01000lbr(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -30 10])
        hold on;
        errorbar(bw,NpiS01000b(:,1),NpiS01000b(:,2),'-o','MarkerSize',10,'Color',[0,0.7,1]);
        errorbar(bw,NpiS01000lr(:,1),NpiS01000lr(:,2),'-dc','MarkerSize',10);
        plot(bw,NpiS0modeldata1000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        legend({'NpiS0 condition lbr','NpiS0 condition b',...
            'NpiS0 condition lr','model data from Breebaart (2001)'},...
            [290,170,0.1,0.05]);
        title('NpiS0 1000 Hz thresholds: monaural factor = 0.0003')
    end

end












