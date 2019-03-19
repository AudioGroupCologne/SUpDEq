function output=exp_breebaart2001(varargin)
%EXP_BREEBAART2001   Figures from Breebaart et al.(2001a and 2001b)
%   Usage: output = exp_breebaart2001(flags)
%
%   EXP_BREEBAART2001(flags) reproduces experiments from the 
%   Breebaart et al (2001a and 2001b) papers.
%
%   The following flags can be specified;
%
%     'a_fig2'   Reproduce Fig. 2 from Breebaart et al. (2001a):
%                Output of the peripheral preprocessor for a 500-Hz tone
%                (left panel) and a 4000-Hz tone(right panel) of 100-ms
%                duration. The output is caluculated for a filter which
%                is nearest to the frequency of the tone.
%
%     'a_fig6'   Reproduce Fig. 6 from Breebaart et al. (2001a):
%                Left panel: Idealized EI-activity for a wideband diotic
%                noise (0-4000 Hz) with an overall level of 70 dB SPL
%                (offset = default = 100 dB) for an auditory filter
%                which is nearest to the frequency of the tone.
%                Right panel: change in the activity pattern of the left
%                pannel if a 500-Hz interaurally ou-of-phase signal
%                (Spi) is added with a level of 50 dB (offset = default
%                = 100 dB). Both signal and master have a duration of 1
%                second and the figures show the mean output of the
%                first 100 ms.
%
%     'b_fig3'   Reproduce Fig. 3 from Breebaart et al. (2001b):
%                N0Spi thresholds as a function of the masker bandwidth
%                for a constant overall level of the masker. The six
%                panels represent center frequencies of 125, 250, 500,
%                1000, 2000 and 4000 Hz, respectivly. The filled squares
%                are model predictions from Breebaart et al. (2001b).
%                The stars are model predictions calculated with our
%                implementation of the Breebaart model for a combination
%                of binaural and monaural decisions. The open squares
%                are data adapted from van de Par and Kohlrausch (1999).
%
%     'b_fig6'   Reproduce Fig. 6 from Breebaart et al. (2001b):
%                NpiSo thresholds as a function of masker bandwidth for
%                125-Hz (upper-left panel), 250-Hz (upper-right-panel),
%                500-Hz (lower left panel), and 1000-Hz center frequency
%                (lower-right panel). The open sybols are data adapted
%                from van de Par and Kohlrausch (1999), the filled
%                symbols are model preictions from Breebaart et al.
%                (2001b) and the stars are model predicitons calculated
%                with out implementation of the Breebaart model for a
%                combination of binaural and monaural decisions.
%
%     'fig1_vandepar1999'         Reproduce Fig.1 from van de Par and
%                                 Kohlrausch for the N0S0 condition:
%                                 N0S0 thresholds as a fuction of the
%                                 masker bandwidth expressed in
%                                 signal-to-overall-noise power ratio.
%                                 The six panels show data at various
%                                 center frequencies. The circels are
%                                 experimental data and the stars are
%                                 model predicitions calculated with our
%                                 implemntation ot the Breebaart model
%                                 using the monaural decision.
%
%
%   Further, cache flags (see amt_cache) and plot flags can be specified:
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%   Figure 3 of Breebaart et al. (2001b) can also be calculated using the
%   interface of the model initiative. Assuming a model interface
%   waiting for binaural signals in a directory DIR, use
%   EXP_BREEBAART2001('b_fig3','redo','ModelInitiative','directory',DIR);. The AMT
%   will start the experiment, collect the responses from the model server,
%   and plot the figure. 
%
%   See also: data_breebaart2001 data_vandepar1999 breebaart2001_centralproc breebaart2001_preproc sig_breebaart2001 
%
%   Example:
%   ---------
%
%   To display Figure 2 of Breebaart et al. (2001a) use :
%
%     out = exp_breebaart2001('a_fig2');
%
%   To display Figure 6 of Breebaart et al. (2001a) use :
%
%     out = exp_breebaart2001('a_fig6');
%
%   To display Figure 3 of Breebaart et al. (2001b) use :
%
%     out = exp_breebaart2001('b_fig3');
%
%   To display Figure 6 of Breebaart et al. (2001b) use :
%
%     out = exp_breebaart2001('b_fig6');
%
%   To display Figure 1 of van de Par and Kohlrausch (1999) use :
%
%     out = exp_breebaart2001('fig1_vandepar1999');
%
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074-1088, August 2001.
%     
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. II. Dependence on spectral
%     parameters. J. Acoust. Soc. Am., 110:1089-1104, August 2001.
%     
%     S. van de Par and A. Kohlrausch. Dependence of binaural masking level
%     differences on center frequency, masker bandwidth, and interaural
%     parameters. J. Acoust. Soc. Am., 106(4):1940-1947, 1999.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_breebaart2001.php

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
  
%  AUTHOR: Martina Kreuzbichler


warning off;
%% ------ Check input options --------------------------------------------

  definput.import={'amt_cache'};
  definput.flags.type = {'missingflag','a_fig2','a_fig6','b_fig3','b_fig6',...
      'fig1_vandepar1999'};
  definput.flags.plot = {'plot','noplot'};
  definput.flags.interface = {'AMT', 'ModelInitiative'};
  definput.keyvals.directory=tempdir;

  % Parse input options
  [flags,kv]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

% Breebaart et al. (2001a) fig. 2
if flags.do_a_fig2
    
    testsigadlo = amt_cache('get','a_fig2',flags.cachemode);
    
    if isempty(testsigadlo)
        definput.import = {'auditoryfilterbank','ihcenvelope','adaptloop','breebaart2001_eicell'};
        definput.importdefaults={'fhigh',8000,'ihc_breebaart','adt_breebaart'};

        [flags,keyvals,~,~,~]  = ltfatarghelper({'flow', 'fhigh', ...
                            'basef'},definput,varargin);


        testsig(:,1) = sin(2*pi*(0:(0.1*48000)-1)'*500/48000);
        testsig(:,2) = sin(2*pi*(0:(0.1*48000)-1)'*4000/48000);

        testsig = setdbspl(testsig,70);
        testsig = [testsig; zeros(0.1*48000,2)];

        for counter = 1:2

            testsigoutmiddle(:,counter) = ...
                breebaart2001_outmiddlefilter(testsig(:,counter),48000);

            [testsigaudfilt(:,:,counter), ~] = ...
                auditoryfilterbank(testsigoutmiddle(:,counter),48000,'argimport',flags,keyvals);

            testsigihc(:,:,counter) = ihcenvelope(testsigaudfilt(:,:,counter),...
                48000,'argimport',flags,keyvals);

            testsigadlo(:,:,counter) = adaptloop(testsigihc(:,:,counter),48000,...
                'argimport',flags,keyvals);
        end
        
        amt_cache('set','a_fig2',testsigadlo)
        
    end
    
    output = testsigadlo;
    
% Breebaart et al. (2001a) fig. 6   
elseif flags.do_a_fig6
    
    [ei_map_n_mean,ei_map_diff_mean] = ...
        amt_cache('get','a_fig6',flags.cachemode);
    
    if isempty(ei_map_n_mean)
    % do computation
    
        definput.import = {'auditoryfilterbank','ihcenvelope','adaptloop','breebaart2001_eicell'};
        definput.importdefaults={'fhigh',8000,'ihc_breebaart','adt_breebaart'};

        [flags,keyvals,flow,fhigh,basef]  = ltfatarghelper({'flow', 'fhigh', ...
                            'basef'},definput,varargin);

        fs = 32000;

        insig(:,1) = sig_bandpassnoise(2000,fs,1,70,4000);
        insig(:,2) = insig(:,1);
        insig = setdbspl(insig,70);
        insig(:,3) = sin(2*pi*(0:(1*fs)-1)'*500/fs);
        insig(:,4) = sin(2*pi*(0:(1*fs)-1)'*500/fs-pi);
        insig(:,3:4) = setdbspl(insig(:,3:4),50);
        insig(:,3:4) = insig(:,3:4) + insig(:,1:2);
        signalcounter = 1;

        for counter = 1:2

            sig = insig(:,signalcounter:signalcounter+1);

            for earnumber = 1:2
                outsig(:,earnumber) = breebaart2001_outmiddlefilter(sig(:,earnumber),fs);
            end

            [outsig, ~] = auditoryfilterbank(outsig,fs,'argimport',flags,keyvals);

            outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

            outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);
            outsig = squeeze(outsig(:,9,:)); % 9 is filter at ~ 520 Hz

            ei_map{counter} = zeros(41,41,length(outsig));

            taucount = 1;
            ildcount = 1;
            for tau = -0.002:0.0001:0.002
                for ild = -10:0.5:10
                    n = round( abs(tau) * fs );
                    l=outsig(:,1);
                    r=outsig(:,2);
                    if tau > 0,
                        l = [zeros(n,1) ; l(1:end-n)];
                    else
                        r = [zeros(n,1) ; r(1:end-n)];
                    end
                    % apply characteristic ILD:
                    l=gaindb(l, ild/2);
                    r=gaindb(r,-ild/2);
                    % compute instanteneous EI output:
                    x = (l - r).^2;
                    % temporal smoothing:
                    A=[1 -exp(-1/(fs*30e-3))];
                    B=[1-exp(-1/(fs*30e-3))];
                    y= filtfilt(B,A,x);

                    z =0.1*log(0.00002*y+1);
                    ei_map{counter}(taucount,ildcount,:) = z; 
                    ildcount = ildcount + 1;
                end
                ildcount = 1;
                taucount = taucount + 1;
            end
            signalcounter = signalcounter + 2;
        end
        ei_map_n_mean = mean(ei_map{1}(:,:,1:3200,1),3); % mean of first 100 ms
        ei_map_sn_mean = mean(ei_map{2}(:,:,1:3200,1),3); % mean of first 100 ms
        ei_map_diff_mean = ei_map_sn_mean-ei_map_n_mean;

        amt_cache('set','a_fig6',ei_map_n_mean,ei_map_diff_mean);
    end
    
    output = struct('ei_map_noise',ei_map_n_mean,'ei_map_diff',ei_map_diff_mean);



% Breebaart et al. (2001b) fig. 3
elseif flags.do_b_fig3
    
    [N0Spi125,N0Spi250,N0Spi500,N0Spi1000,N0Spi2000,N0Spi4000] = ...
        amt_cache('get','b_fig3',flags.cachemode);
    
    if isempty(N0Spi125)
    % do computation
        
        parout = [];
        switch flags.interface
          case 'AMT'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
              'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',65};
          case 'ModelInitiative'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',65, ...
                'interface','ModelInitiative','directory',kv.directory,'fs',32000};
        end
        parout = emuexp('expinit',parout,expset);
        
        % input2 = fs; input3 = tau; input4 = ild; 
        modelset = {'name','breebaart2001_preproc','input1',...
            'expsignal', 'input2',32000,'input3',0,'input4',0,...
            'outputs',[1 3 4]};
        parout = emuexp('modelinit',parout,modelset);
        
        decisionset = {'name','breebaart2001_centralproc','input1',...
            'modelout','input2','modelout', 'input3','modelout',...
            'input4','lbr'};
        parout = emuexp('decisioninit',parout,decisionset);

        fc = [125 250 500 1000 2000 4000];
        bw = [5 10 25 50 100 250 500 1000 2000 4000];
        %bwadd = [500 1000 2000 4000];
        nl = 65;

        % loop for all center freq.
        for fccount = 1:length(fc)

            %loop for all bandwidths
            for bwcount = 1:length(bw)
                  % calculate only for bandwidths max. twice as large as fc
                if bw(bwcount) <= 2*fc(fccount)
                  %input3 = signallevel; input4 = signalduration;
                  %input5 = signalphase; input7 = noiselevel;
                  %input8 = noiseduration; input9 = noisephase;
                  %input10 = hanning ramp duration; input11 = fs;
                  signalset = {'name','sig_breebaart2001','input1',...
                      'inttyp', 'input2',fc(fccount),'input3',...
                      'expvar','input4',0.3,'input5',pi,'input6',...
                      bw(bwcount),'input7',nl,'input8',0.4,'input9',0,...
                      'input10',0.05,'input11', 32000};
                  parout = emuexp('signalinit',parout,signalset);

                  % loop for six experimental runs
                  resultbwvec = zeros(6,1);
                  for runcounter = 1:6
                      amt_disp(['Calculating: center frequency = ', num2str(fc(fccount)), ...
                          ' Hz, bandwidth = ' num2str(bw(bwcount)) ' Hz, run #' num2str(runcounter)],'progress');

                      result = emuexp('run',parout);
                      resultbwvec(runcounter) = result(1)-nl;
                  end
                  resultvec(bwcount) = mean(resultbwvec);
                  resultvecstd(bwcount) = std(resultbwvec,1);                  
                end
            end

            N0Spi_temp =sprintf('N0Spi%i',fc(fccount));
            output.(N0Spi_temp) = [resultvec; resultvecstd]';

            %if fccount <= length(bwadd)
            %    bw(end+1) = bwadd(fccount);
            %end

        end
        
        N0Spi125 = output.N0Spi125;
        N0Spi250 = output.N0Spi250;
        N0Spi500 = output.N0Spi500;
        N0Spi1000 = output.N0Spi1000;
        N0Spi2000 = output.N0Spi2000;
        N0Spi4000 = output.N0Spi4000;
        
        switch flags.interface
          case 'AMT'
            amt_cache('set','b_fig3',N0Spi125,N0Spi250,N0Spi500,N0Spi1000,...
                N0Spi2000,N0Spi4000);
        end
   
    else
        output = struct('N0Spi125',N0Spi125,'N0Spi250',N0Spi250,...
            'N0Spi500',N0Spi500,'N0Spi1000',N0Spi1000,'N0Spi2000',...
            N0Spi2000,'N0Spi40000',N0Spi4000); 
    end

elseif flags.do_b_fig6
    
    [NpiS0125,NpiS0250,NpiS0500,NpiS01000] = ...
        amt_cache('get','b_fig6',flags.cachemode);
    
    if isempty(NpiS0125)
    % do computation
        
        parout = [];
        switch flags.interface
          case 'AMT'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
              'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',65};
          case 'ModelInitiative'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',65, ...
                'interface','ModelInitiative','directory',kv.directory,'fs',32000};
        end
        parout = emuexp('expinit',parout,expset);

        decisionset = {'name','breebaart2001_centralproc','input1',...
            'modelout','input2','modelout', 'input3','modelout',...
            'input4','lbr'};
        parout = emuexp('decisioninit',parout,decisionset);

        fc = [125 250 500 1000];
        bw = [5 10 25 50 100 250];
        bwadd = [500 1000 2000];
        tau = [0.0039 0.002 0.001 0.0005];
        nl = 70;

        % loop for all center freq.
        for fccount = 1:length(fc)
            
            % input2 = fs; input3 = tau; input4 = ild; 
            modelset = {'name','breebaart2001_preproc','input1',...
                'expsignal','input2',32000,'input3',tau(fccount),...
                'input4',0,'outputs',[1 3 4]};
            parout = emuexp('modelinit',parout,modelset);

            %loop for all bandwidths
            for bwcount = 1:length(bw)
                
                %input3 = signallevel; input4 = signalduration;
                %input5 = signalphase; input7 = noiselevel;
                %input8 = noiseduration; input9 = noisephase;
                %input10 = hanning ramp duration; input11 = fs;
                signalset = {'name','sig_breebaart2001','input1',...
                    'inttyp', 'input2',fc(fccount),'input3',...
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
                amt_disp(sprintf(['Progress for %i Hz center frequency: ' ...
                    num2str(round(bwcount/length(bw)*100)) ...
                    '%% calculated'],fc(fccount)),'progress'); 

            end

            NpiS0_temp =sprintf('NpiS0%i',fc(fccount));
            output.(NpiS0_temp) = [resultvec; resultvecstd]';

            if fccount <= length(bwadd)
                bw(end+1) = bwadd(fccount);
            end

        end
        
        NpiS0125 = output.NpiS0125;
        NpiS0250 = output.NpiS0250;
        NpiS0500 = output.NpiS0500;
        NpiS01000 = output.NpiS01000;

        switch flags.interface
          case 'AMT'
            amt_cache('set','b_fig6',NpiS0125,NpiS0250,NpiS0500,NpiS01000);
        end    
    
    else
    output = struct('NpiS0125',NpiS0125,'NpiS0250',NpiS0250,'NpiS0500',...
        NpiS0500,'NpiS01000',NpiS01000); 
    end 
    
elseif flags.do_fig1_vandepar1999
    
    [N0S0125,N0S0250,N0S0500,N0S01000,N0S02000,N0S04000] = ...
        amt_cache('get','fig1_vandepar1999',flags.cachemode);
    
    if isempty(N0S0125)
        parout = [];
        switch flags.interface
          case 'AMT'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',90};
          case 'ModelInitiative'
            expset = {'intnum',3,'rule',[2 1],'expvarstepstart',8,...
                'expvarsteprule',[0.5 2],'stepmin',[1 8],'expvarstart',90, ...
                'interface','ModelInitiative','directory',kv.directory,'fs',32000};
        end
        parout = emuexp('expinit',parout,expset);
        
        % input2 = fs; input3 = tau; input4 = ild; 
        modelset = {'name','breebaart2001_preproc','input1',...
            'expsignal', 'input2',32000,'input3',0,'input4',0,...
            'outputs',[1 3 4]};
        parout = emuexp('modelinit',parout,modelset);
        
        decisionset = {'name','breebaart2001_centralproc','input1',...
            'modelout','input2','modelout', 'input3','modelout',...
            'input4','lbr'};
        parout = emuexp('decisioninit',parout,decisionset);

        fc = [125 250 500 1000 2000 4000];
        bw = [5 10 25 50 100 250];
        bwadd = [500 1000 2000 4000];
        nl = 70;

        % loop for all center freq.
        for fccount = 1:length(fc)

            %loop for all bandwidths
            for bwcount = 1:length(bw)

                %input3 = signallevel; input4 = signalduration;
                %input5 = signalphase; input7 = noiselevel;
                %input8 = noiseduration; input9 = noisephase;
                %input10 = hanning ramp duration; input11 = fs;
                signalset = {'name','sig_breebaart2001','input1',...
                    'inttyp', 'input2',fc(fccount),'input3',...
                    'expvar','input4',0.3,'input5',0,'input6',...
                    bw(bwcount),'input7',nl,'input8',0.4,'input9',0,...
                    'input10',0.05,'input11', 32000};
                parout = emuexp('signalinit',parout,signalset);

                % loop for experimental runs
                for runcounter = 1:6
                    amt_disp(['Calculating: center frequency = ', num2str(fc(fccount)), ...
                        ' Hz, bandwidth = ' num2str(bw(bwcount)) ' Hz, run #' num2str(runcounter)],'progress');                  
                    result = emuexp('run',parout);
                    resultbwvec(runcounter) = result(1)-nl;
                end

                resultvec(bwcount) = mean(resultbwvec);
                resultvecstd(bwcount) = std(resultbwvec,1);
                resultbwvec = zeros(6,1);
            end

            N0S0_temp =sprintf('N0S0%i',fc(fccount));
            output.(N0S0_temp) = [resultvec; resultvecstd]';

            if fccount <= length(bwadd)
                bw(end+1) = bwadd(fccount);
            end

        end
        
        N0S0125 = output.N0S0125;
        N0S0250 = output.N0S0250;
        N0S0500 = output.N0S0500;
        N0S01000 = output.N0S01000;
        N0S02000 = output.N0S02000;
        N0S04000 = output.N0S04000;

        switch flags.interface
          case 'AMT'
            amt_cache('set','fig1_vandepar1999',N0S0125,N0S0250,N0S0500,...
                N0S01000,N0S02000,N0S04000);
        end
        
    else
        output = struct('N0S0125',N0S0125,'N0S0250',N0S0250,...
            'N0S0500',N0S0500,'N0S01000',N0S01000,'N0S02000',N0S02000,...
            'N0S040000',N0S04000);  
    end

    
end


if flags.do_plot
    if flags.do_a_fig2
        fax = 0:0.2/9599:0.2;
        hafig2 = figure;
        set(hafig2, 'Position', [400 400 1300 450])
        subplot(1,2,1)
        plot(fax,output(:,9,1))
        axis([0 0.2 -1000 8500])
        set(gca,'XTick',[0 0.05 0.1 0.15 0.2]);
        xlabel('time [s]')
        ylabel('Output [MU')
        text('Position',[0.1 7000],'string','500 Hz tone')
        subplot(1,2,2)
        plot(fax,output(:,25,2))
        axis([0 0.2 -1000 8500])
        set(gca,'XTick',[0 0.05 0.1 0.15 0.2]);
        xlabel('time [s]')
        ylabel('Output [MU]')
        text('Position',[0.1 7000],'string','4000 Hz tone')
        
        set(gcf,'NextPlot','add');
        axes;
        h = title('Fig. 2 Output of the peripheral preprocessor');
        set(gca,'Visible','off');
        set(h,'Visible','on'); 

        
    elseif flags.do_a_fig6
        
        [x1,y1] = meshgrid(-10:0.5:10, -2:0.1:2);
        
        hfig1=figure;
        set(hfig1, 'Position', [100 100 1400 500])
        subplot(1,2,1)
        surf(x1,y1,ei_map_n_mean)
        xlabel('\alpha [dB]')
        ylabel('\tau [ms]')
        zlabel('Activity [MU]')
        axis([-10 10 -2 2 0 0.4])
        set(gca,'XTick',[-10 0 10]);
        set(gca,'YDir','reverse')
        set(gca,'YTick',[-2 0 2]);
        set(gca,'ZTick',[0 0.2 0.4]);
        hcb1=colorbar;
        caxis([0 0.4]);
        set(hcb1,'YTick',[0:0.1:0.4]);
        grid off
        view(-53,28)
        
        subplot(1,2,2)
        surf(x1,y1,ei_map_diff_mean)
        xlabel('\alpha [dB]')
        ylabel('\tau [ms]')
        zlabel('Activity [MU]')
        axis([-10 10 -2 2 -0.05 0.211])
        set(gca,'XTick',[-10 0 10]);
        set(gca,'YDir','reverse')
        set(gca,'YTick',[-2 0 2]);
        set(gca,'ZTick',[0 0.1 0.2]);
        hcb2 = colorbar;
        caxis([0 0.2]);
        set(hcb2,'YTick',[0:0.05:0.2]);
        grid off
        view(-53,28)
        set(gca,'LooseInset',get(gca,'TightInset'))
        
        set(gcf,'NextPlot','add');
        axes;
        h = title('Fig. 6');
        set(gca,'Visible','off');
        set(h,'Visible','on'); 


        
        
    elseif flags.do_b_fig3
        
        % get experimental data
        N0Spiexpdata125 = data_vandepar1999('fig1_N0Spi','nfc125');
        N0Spiexpdata250 = data_vandepar1999('fig1_N0Spi','nfc250');
        N0Spiexpdata500 = data_vandepar1999('fig1_N0Spi','nfc500');
        N0Spiexpdata1000 = data_vandepar1999('fig1_N0Spi','nfc1000');
        N0Spiexpdata2000 = data_vandepar1999('fig1_N0Spi','nfc2000');
        N0Spiexpdata4000 = data_vandepar1999('fig1_N0Spi','nfc4000');
        
        % get model data
        N0Spimodeldata125 = data_breebaart2001('fig3','nfc125');
        N0Spimodeldata250 = data_breebaart2001('fig3','nfc250');
        N0Spimodeldata500 = data_breebaart2001('fig3','nfc500');
        N0Spimodeldata1000 = data_breebaart2001('fig3','nfc1000');
        N0Spimodeldata2000 = data_breebaart2001('fig3','nfc2000');
        N0Spimodeldata4000 = data_breebaart2001('fig3','nfc4000');

        hfig2=figure;
        set(hfig2, 'Position', [100 100 1100 700])
        
        bw = [5 10 25 50 100 250];
        subplot(3,2,1);
        errorbar(bw,N0Spi125(:,1),N0Spi125(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata125,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,N0Spiexpdata125,'-sg','MarkerSize',10)
        text('Position',[1000 -10],'string','125 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500];
        subplot(3,2,2);
        errorbar(bw,N0Spi250(:,1),N0Spi250(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata250,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,N0Spiexpdata250,'-sg','MarkerSize',10)
        text('Position',[1000 -10],'string','250 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500 1000];
        subplot(3,2,3);
        errorbar(bw,N0Spi500(:,1),N0Spi500(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata500,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,N0Spiexpdata500,'-sg','MarkerSize',10)
        text('Position',[1000 -10],'string','500 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500 1000 2000];
        subplot(3,2,4);
        errorbar(bw,N0Spi1000(:,1),N0Spi1000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata1000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,N0Spiexpdata1000,'-sg','MarkerSize',10)
        text('Position',[1000 -10],'string','1000 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        
        bw = [5 10 25 50 100 250 500 1000 2000 4000];
        subplot(3,2,5);
        errorbar(bw,N0Spi2000(:,1),N0Spi2000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata2000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,N0Spiexpdata2000,'-sg','MarkerSize',10)
        text('Position',[1000 -10],'string','2000 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        
        subplot(3,2,6);
        errorbar(bw,N0Spi4000(:,1),N0Spi4000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 5])
        hold on;
        plot(bw,N0Spimodeldata4000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,N0Spiexpdata4000,'-sg','MarkerSize',10)
        text('Position',[1000 0],'string','4000 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        legend({'modeled data SSL = 65 dB, NL 65 dB  - monaural factor = 0.0003',...
            'model data from Breebaart (2001)',...
            'experimental data from van de Par and Kohlrausch (1999)'},...
            'Position' ,[0.75,0.943,0.1,0.04]);
        set(gcf,'NextPlot','add');
        axes;
        h = title('Fig. 3 N0Spi Thresholds');
        set(gca,'Visible','off');
        set(h,'Visible','on'); 
        
        
    elseif flags.do_b_fig6
        
        % get experimental data
        NpiS0expdata125 = data_vandepar1999('fig1_NpiS0','nfc125');
        NpiS0expdata250 = data_vandepar1999('fig1_NpiS0','nfc250');
        NpiS0expdata500 = data_vandepar1999('fig1_NpiS0','nfc500');
        NpiS0expdata1000 = data_vandepar1999('fig1_NpiS0','nfc1000');

        
        % get model data
        NpiS0modeldata125 = data_breebaart2001('fig6','nfc125');
        NpiS0modeldata250 = data_breebaart2001('fig6','nfc250');
        NpiS0modeldata500 = data_breebaart2001('fig6','nfc500');
        NpiS0modeldata1000 = data_breebaart2001('fig6','nfc1000');

        hfig3=figure;
        set(hfig3, 'Position', [100 100 1100 700])
        
        bw = [5 10 25 50 100 250];
        subplot(2,2,1);
        errorbar(bw,NpiS0125(:,1),NpiS0125(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 0])
        hold on;
        plot(bw,NpiS0modeldata125,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,NpiS0expdata125,'-sg','MarkerSize',10)
        text('Position',[5 -30],'string','125 Hz; 3.9 ms')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500];
        subplot(2,2,2);
        errorbar(bw,NpiS0250(:,1),NpiS0250(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,NpiS0modeldata250,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,NpiS0expdata250,'-sg','MarkerSize',10)
        text('Position',[5 -30],'string','250 Hz; 2 ms')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500 1000];
        subplot(2,2,3);
        errorbar(bw,NpiS0500(:,1),NpiS0500(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,NpiS0modeldata500,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,NpiS0expdata500,'-sg','MarkerSize',10)
        text('Position',[5 -30],'string','500 Hz; 1 ms')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N');

        bw = [5 10 25 50 100 250 500 1000 2000];
        subplot(2,2,4);
        errorbar(bw,NpiS01000(:,1),NpiS01000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,NpiS0modeldata1000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,NpiS0expdata1000,'-sg','MarkerSize',10)
        text('Position',[5 -30],'string','1000 Hz; 0.5 ms')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N');

        legend({'modeled data SSL = 85 dB, NL 70 dB  - monaural factor = 0.0003',...
            'model data from Breebaart (2001)',...
            'experimental data from van de Par and Kohlrausch (1999)'},...
            'Position' ,[0.75,0.943,0.1,0.04]);
        set(gcf,'NextPlot','add');
        axes;
        h = title('Fig. 6 NpiS0 Thresholds');
        set(gca,'Visible','off');
        set(h,'Visible','on'); 
        
        
        
    elseif flags.do_fig1_vandepar1999
        
        %get experimental data
        N0S0expdata125 = data_vandepar1999('fig1_N0S0','nfc125');
        N0S0expdata250 = data_vandepar1999('fig1_N0S0','nfc250');
        N0S0expdata500 = data_vandepar1999('fig1_N0S0','nfc500');
        N0S0expdata1000 = data_vandepar1999('fig1_N0S0','nfc1000');
        N0S0expdata2000 = data_vandepar1999('fig1_N0S0','nfc2000');
        N0S0expdata4000 = data_vandepar1999('fig1_N0S0','nfc4000');
        
        hfig1=figure;
        set(hfig1, 'Position', [100 100 700 700])
        
        bw = [5 10 25 50 100 250];
        s1 = subplot(3,2,1);
        p1 = get(s1,'pos');
        errorbar(bw,N0S0125(:,1),N0S0125(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([4 5000 -30 10])
        hold on;
        plot(bw,N0S0expdata125,'-og')
        text('Position',[1000 5],'string','125 Hz')
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500];
        s2 = subplot(3,2,2);
        p2 = get(s2,'pos');
        p2(1) = p1(1) + p1(3);
        set(s2, 'pos', p2);
        errorbar(bw,N0S0250(:,1),N0S0250(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([4 5000 -30 10])
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        hold on;
        plot(bw,N0S0expdata250,'-og')
        text('Position',[1000 5],'string','250 Hz')

        bw = [5 10 25 50 100 250 500 1000];
        s3 = subplot(3,2,3);
        p3 = get(s3,'pos');
        p3(2) = p1(2) - p1(4) + 0.005;
        set(s3, 'pos', p3);
        errorbar(bw,N0S0500(:,1),N0S0500(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([1 10000 -30 10])
        set(gca,'XTick',[]);
        set(gca,'YTick',[-30 -20 -10 0]);
        hold on;
        plot(bw,N0S0expdata500,'-og')
        text('Position',[1000 5],'string','500 Hz')
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500 1000 2000];
        s4 = subplot(3,2,4);
        p4 = get(s4,'pos');
        p4(1) = p3(1) + p3(3);
        p4(2) = p3(2);
        set(s4, 'pos', p4);
        errorbar(bw,N0S01000(:,1),N0S01000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([1 10000 -30 10])
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        hold on;
        plot(bw,N0S0expdata1000,'-og')
        text('Position',[1000 5],'string','1 kHz')
        xlabel('Bandwidth [Hz]');

        bw = [5 10 25 50 100 250 500 1000 2000 4000];
        s5 = subplot(3,2,5);
        p5 = get(s5,'pos');
        p5(2) = p3(2) - p3(4) + 0.005;
        set(s5, 'pos', p5);
        errorbar(bw,N0S02000(:,1),N0S02000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[1 10 100 1000 10000]);
        axis([1 10000 -30 10])
        set(gca,'YTick',[-30 -20 -10 0]);
        hold on;
        plot(bw,N0S0expdata2000,'-og')
        text('Position',[1000 5],'string','2 kHz')
        xlabel('Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        s6 = subplot(3,2,6);
        p6 = get(s6,'pos');
        p6(1) = p5(1) + p5(3);
        p6(2) = p5(2);
        set(s6, 'pos', p6);
        errorbar(bw,N0S04000(:,1),N0S04000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[ 10 100 1000 10000]);
        axis([1 10000 -30 10])
        set(gca,'YTick',[]);
        hold on;
        plot(bw,N0S0expdata4000,'-og')
        text('Position',[1000 5],'string','4 kHz')
        xlabel('Bandwidth [Hz]');
        
        legend({'modeled data SSL = 90 dB, NL 70 dB  - monaural factor = 0.0003',...
            'experimental data from van de Par and Kohlrausch (1999)'},...
            'Position' ,[p2(1)-p2(3)*0.2,0.15,0.1,0.05]);
        set(gcf,'NextPlot','add');
        axes;
        h = title('Fig. 1 N_0S_0 Thresholds');
        set(gca,'Visible','off');
        set(h,'Visible','on'); 
           
    end
           
end
    
    
