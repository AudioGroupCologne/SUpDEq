function varargout = exp_laback2023(varargin)
%EXP_LABACK2023 - Figures from Laback (2023)
%
%   Usage: exp_laback2023('fig5')
%
%   EXP_LABACK2023 predicts the effect of different types of precursors 
%   on ILD-based perceived lateralization. It is based on AN model by 
%   Smalt, Heinz and Strickland (2014), which is a binaural version of the  
%   Auditory Nerve Model by Zilany and Bruce (JASA 2006, 2007), that 
%   involves ipsi- and contralateral cochlear compression control via 
%   efferent feedback.
%
%
%   Examples:
%   ---------
%
%   To display Fig.5 use :
%
%     exp_laback2023('fig5');
%
%   To run 'experiment 2' (with model zilany2014), condition 1, use
%     exp_laback2023('fig5','i_cond',1);
%
%   See also: laback2023 smalt2014 zilany2014
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_laback2023.php


% #Author: Bernhard Laback, 2022
% #Author: Clara Hollomey, adaptations for AMT, 2023
% #Author: Alejandro Osses, further adaptations for AMT 1.4, 2023

definput.flags.type = {'missingflag','fig5'};
definput.keyvals.i_cond = [1 2 3 4];
[flags,kv]  = ltfatarghelper({},definput,varargin);

idxs_cond = kv.i_cond;

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

ILD_vector = 0:2:10; % vector of ILD's
num_ILD = length(ILD_vector); % remove kv.NrILD
% remove: kv.step
modelreps = 3;

if flags.do_fig5
    idxs_model = [1 2]; % runs both models
end
if flags.do_fig5 || flags.do_fig5ab || flags.do_fig5cd
    for i_model = idxs_model
        definput.import = {'laback2023'};
        if i_model == 1
            definput.importdefaults={'exp1'};
            SubExp = 0; % used in sig_laback2023.m for choosing the param set
        else
            definput.importdefaults={'exp2'};
            SubExp = 0; % used in sig_laback2023.m for choosing the param set
        end

        [fg,kv]  = ltfatarghelper({},definput,varargin); % obtained the keyvalues

        % Taking some variables out of the keyvalues (easier for debugging)
        model = kv.model; 
        fs = kv.fs;
        % end Taking

        if kv.PrecILD < 0
            error('Precursor ILD has to be >= 0');
        end
        if kv.PrecITD < 0
            error('Precursor ITD has to be >= 0!');
        end

        switch model
            case 1
                model_str = 'smalt2014';
                Modeldelay = 240;   %Delay between stimulus and model response (in samples) [default: 240]
            case 2
                model_str = 'zilany2014';
                Modeldelay = 400;   %Delay between stimulus and model response (in samples) [default: 240]
        end

        % Define time vectors for P and T
        nvec = 0: 1/fs : kv.DurN-1/fs; 
        tvec = 0: 1/fs : kv.DurT-1/fs;
        
        % Memory allocation for results matrices
        ResultMatrixvecLOW  = zeros(modelreps,8,num_ILD);
        RespLMatrixLOW      = zeros(modelreps,8,num_ILD);
        RespRMatrixLOW      = zeros(modelreps,8,num_ILD);

        ResultMatrixvecMID  = zeros(modelreps,8,num_ILD);
        RespLMatrixMID      = zeros(modelreps,8,num_ILD);
        RespRMatrixMID      = zeros(modelreps,8,num_ILD);

        ResultMatrixvecHIGH = zeros(modelreps,8,num_ILD);
        RespLMatrixHIGH     = zeros(modelreps,8,num_ILD);
        RespRMatrixHIGH     = zeros(modelreps,8,num_ILD);

        NrFiberTypes = length(kv.spontvec);
        for i_fibtype=1:NrFiberTypes %three different fiber types (based on SR)
            % spontaneous rate (in spikes/s) of the fiber BEFORE 
            % refractory effects stimulus parameters (default: 30)
            switch i_fibtype
                case 1
                    i_fibtype_str = 'low-SR';
                case 2
                    i_fibtype_str = 'mid-SR';
                case 3
                    i_fibtype_str = 'high-SR';
            end
            amt_disp(sprintf('Model=%s; Fibertype=%s', model_str, i_fibtype_str));
            for i_modelreps=1:modelreps %loop for repetition
                amt_disp(sprintf('\tRepetition %.0f of %.0f', i_modelreps,modelreps));
                
                % Initialisation of the variables containing the results per 
                %   condition and ILD.
                ResultMatrix = nan(8,num_ILD); % 8=2*4 conditions
                RespLMatrix  = nan(size(ResultMatrix));
                RespRMatrix  = nan(size(ResultMatrix));
                  
                for i_cond = idxs_cond %loop for four different stimulus conditions (1: no T, 2:diotic P, 3: ipsi P, 4: contra P)

                    RespL = zeros(2, num_ILD);
                    RespR = zeros(2, num_ILD);
                    diffvec = zeros(2, num_ILD);
                    
                    %Define precursor conditions
                    switch i_cond
                        case 1 % no precursor
                            precmode = 0; %precursor mode: 0 = without precursor, 1 = with precursor 
                            ILDNoise = 0; %ILD of precursor noise, positive value favoring left ear
                            ILDNoiseLin = 10^(ILDNoise/20);
                            PrecITD_side = [];
                            PrecITD_t= [];
                        case 2 % diotic precursor
                            precmode = 1;  
                            ILDNoise = 0; 
                            ILDNoiseLin = 10^(ILDNoise/20);
                            PrecITD_side = 0;
                            PrecITD_t=0;
                        case 3 % precursor favoring LEFT ear (ipsilateral, because target is always left leading)
                            precmode = 1; 
                            ILDNoise = kv.PrecILD;
                            ILDNoiseLin = 10^(ILDNoise/20);
                            PrecITD_side = -1;
                            PrecITD_t=kv.PrecITD;
                        case 4 % precursor favoring right ear (contralater, because target is always left leading)
                            precmode = 1; 
                            ILDNoise = kv.PrecILD*-1; 
                            ILDNoiseLin = 10^(ILDNoise/20);
                            PrecITD_side = 1;
                            PrecITD_t=kv.PrecITD;
                    end
                    % clear stim RespL RespR diffvec timeout meout mocr c1filterout c2filterout c1vihc c2vihc vihc psthbins psthtime pr psth wstart wend;

                    stim = []; 
                    
                    %TODO:get stimuli sorted
                    %[stimTarg, stimPrec, noiseILDvec, TotDur, B_P, A_P, PrecNoise] = sig_laback2023(precmode,...
                    %    kv.frozennoise, kv.Experiment, tvec, nvec, kv.CF, kv.Fs, kv.F0cue,...
                    %    kv.RampDur, kv.stimdb, kv.DurT, kv.DurN, kv.DurG, kv.factor, ILDNoiseLin,...
                    %    PrecITD_side, PrecITD_t);
                    [stimTarg, ~, noiseILDvec, TotDur, B_P, A_P, PrecNoise] = sig_laback2023(precmode,...
                        tvec, nvec, ILDNoiseLin, PrecITD_side, PrecITD_t, model, 'SubExp',SubExp);
                    
                    if kv.AdaptMode == 1
                        %Generate Pre-precursor
                        %PrecPrecNoise=randn(length(kv.preprecvec),1);
                        PrecPrec = filter(B_P, A_P, PrecNoise);
                        PrecPrec = local_raisedcosine(PrecPrec,fs,kv.RampDur);


                        %Add ILD to PrecPrec
                        PrecPrecILDLin = 10^(kv.PrecPrecILD/20);
                        PrecPrec = PrecPrec * rms(stimTarg)/rms(PrecPrec);
                        PrecPrecVec(1,:) = PrecPrec * sqrt(PrecPrecILDLin);
                        PrecPrecVec(2,:) = PrecPrec / sqrt(PrecPrecILDLin);

                        %TotDur=DurT+DurN+DurG+DurN+DurG; %Duration of stimulus
                    end

                    %% Run the models
                    for MOCstat=1:2 %Moc off vs. on modes
                        if MOCstat==1
                            flags_extra = {'no_MOC'};
                            % if model == 2
                            %     flags_extra(end+1:end+2) = {'ANmode',3}; % default back-end
                            % end
                        elseif MOCstat==2
                            flags_extra = {'MOC'};
                            % if model == 2
                            %     flags_extra(end+1:end+2) = {'ANmode',1}; % LSO from Ashida2016
                            % end
                        end

                        for i_ILD=1:num_ILD
                            ILD = ILD_vector(i_ILD); % dB, for calculating all ILDs around zero: add -12
                            ILD_linear = 10^(ILD/20);
                            
                            if  precmode == 0 % used in cond=1, zilany2014
                                %Apply ILD

                                switch kv.AdaptMode
                                    case 0
                                        stim(:,1) = stimTarg*sqrt(ILD_linear); %positive ILD value enhances left ear
                                        stim(:,2) = stimTarg/sqrt(ILD_linear);
                                    case 1
                                        error('Not validated yet by AO') % convert first to column vectors
                                        stim(1,:) = [PrecPrecVec(1,:) gapvecLong stimTarg'*sqrt(ILD_linear)]; %stimnoise
                                        stim(2,:) = [PrecPrecVec(2,:) gapvecLong stimTarg'/sqrt(ILD_linear)];
                                end

                            elseif  precmode == 1 % used in cond=2, 3, 4, zilany2014

                                %Apply ILD
                                targL = stimTarg*sqrt(ILD_linear); %positive ILD value enhances left ear
                                targR = stimTarg/sqrt(ILD_linear);

                                if kv.AdaptMode == 0
                                    stim(:,1) = [noiseILDvec(:,1); kv.gapvec; targL]; %stimnoise
                                    stim(:,2) = [noiseILDvec(:,2); kv.gapvec; targR];
                                    
                                elseif kv.AdaptMode ==1
                                    error('Not validated yet by AO') 
                                    stim(:,1) = [PrecPrecVec(1,:) kv.gapvecLong noiseILDvec(1,:) kv.gapvec targL']; %stimnoise
                                    stim(:,2) = [PrecPrecVec(2,:) kv.gapvecLong noiseILDvec(2,:) kv.gapvec targR'];
                                end                  
                            end

                            output = laback2023(stim,fs,'argimport',fg,kv,'fiberType',i_fibtype,'precmode',precmode, ...
                                'Modeldelay',Modeldelay,'Onset',kv.Onset, ...
                                'noiseILDvec',noiseILDvec,'length_stimTarg',length(stimTarg), ... % flags from the stimuli
                                flags_extra{:});
                            
                            % Nothing to do...already implemented in laback2023.m
                            RespL(MOCstat, i_ILD)   = output.RespL;
                            RespR(MOCstat, i_ILD)   = output.RespR;
                            diffvec(MOCstat, i_ILD) = output.diffvec;
                            
                            disp('')
                        end % for i_ILD
                        disp('')
                    end % for MOCstat
                    disp('')

                    switch i_cond
                        case 1
                            idx2assign = [1 2];
                        case 2
                            idx2assign = [3 4];
                        case 3
                            idx2assign = [5 6];
                        case 4
                            idx2assign = [7 8];
                    end
                    ResultMatrix(idx2assign,:) = diffvec(1:2,:);
                    RespLMatrix(idx2assign,:) = RespL(1:2,:);
                    RespRMatrix(idx2assign,:) = RespR(1:2,:);

                end % i_cond

                %different P conditions
                switch i_fibtype
                    case 1
                      ResultMatrixvecLOW(i_modelreps,:,:)=ResultMatrix;
                      if kv.model == 2
                          ResultMatrixvecLOW_LSO(i_modelreps,:,:)=ResultMatrix;
                      end
                      RespLMatrixLOW(i_modelreps,:,:)=RespLMatrix;
                      RespRMatrixLOW(i_modelreps,:,:)=RespRMatrix;
                    case 2
                      ResultMatrixvecMID(i_modelreps,:,:)=ResultMatrix;
                      if kv.model == 2
                          ResultMatrixvecMID_LSO(i_modelreps,:,:)=ResultMatrix;
                      end
                      RespLMatrixMID(i_modelreps,:,:)=RespLMatrix;
                      RespRMatrixMID(i_modelreps,:,:)=RespRMatrix;
                    case 3
                      ResultMatrixvecHIGH(i_modelreps,:,:)=ResultMatrix;
                      if kv.model == 2
                          ResultMatrixvecHIGH_LSO(i_modelreps,:,:)=ResultMatrix;
                      end
                      RespLMatrixHIGH(i_modelreps,:,:)=RespLMatrix;
                      RespRMatrixHIGH(i_modelreps,:,:)=RespRMatrix;
                end

            end %repetition
        end %different fiber types

        %% Calculate Mean and SD across model repetitions
        Results_Ave_LOW = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_LOW = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_Ave_MID = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_MID = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_Ave_HIGH = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_HIGH = zeros(i_ILD, length(ResultMatrix(:,1)));

        Results_Ave_LOW_L = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_LOW_L = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_Ave_MID_L = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_MID_L = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_Ave_HIGH_L = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_HIGH_L = zeros(i_ILD, length(ResultMatrix(:,1)));

        Results_Ave_LOW_R = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_LOW_R = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_Ave_MID_R = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_MID_R = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_Ave_HIGH_R = zeros(i_ILD, length(ResultMatrix(:,1)));
        Results_STD_HIGH_R = zeros(i_ILD, length(ResultMatrix(:,1)));

        Results_Ave_LOW_LSO = [];
        Results_STD_LOW_LSO = [];
        Results_Ave_MID_LSO = [];
        Results_STD_MID_LSO = [];
        Results_Ave_HIGH_LSO = [];
        Results_STD_HIGH_LSO = [];

        if modelreps ==1 error('Attention: Error in calculating mean'); end
        %for ILDvec=1:ii
        for ILDidx=1:i_ILD
            for condvec=1:length(ResultMatrix(:,1))

                Results_Ave_LOW(ILDidx,condvec) = mean(ResultMatrixvecLOW(:,condvec,ILDidx));
                Results_STD_LOW(ILDidx,condvec) = std(ResultMatrixvecLOW(:,condvec,ILDidx));
                Results_Ave_MID(ILDidx,condvec) = mean(ResultMatrixvecMID(:,condvec,ILDidx));
                Results_STD_MID(ILDidx,condvec) = std(ResultMatrixvecMID(:,condvec,ILDidx));
                Results_Ave_HIGH(ILDidx,condvec) = mean(ResultMatrixvecHIGH(:,condvec,ILDidx));
                Results_STD_HIGH(ILDidx,condvec) = std(ResultMatrixvecHIGH(:,condvec,ILDidx));

                Results_Ave_LOW_L(ILDidx,condvec) = mean(RespLMatrixLOW(:,condvec,ILDidx));
                Results_STD_LOW_L(ILDidx,condvec) = std(RespLMatrixLOW(:,condvec,ILDidx));
                Results_Ave_MID_L(ILDidx,condvec) = mean(RespLMatrixMID(:,condvec,ILDidx));
                Results_STD_MID_L(ILDidx,condvec) = std(RespLMatrixMID(:,condvec,ILDidx));
                Results_Ave_HIGH_L(ILDidx,condvec) = mean(RespLMatrixHIGH(:,condvec,ILDidx));
                Results_STD_HIGH_L(ILDidx,condvec) = std(RespLMatrixHIGH(:,condvec,ILDidx));

                Results_Ave_LOW_R(ILDidx,condvec) = mean(RespRMatrixLOW(:,condvec,ILDidx));
                Results_STD_LOW_R(ILDidx,condvec) = std(RespRMatrixLOW(:,condvec,ILDidx));
                Results_Ave_MID_R(ILDidx,condvec) = mean(RespRMatrixMID(:,condvec,ILDidx));
                Results_STD_MID_R(ILDidx,condvec) = std(RespRMatrixMID(:,condvec,ILDidx));
                Results_Ave_HIGH_R(ILDidx,condvec) = mean(RespRMatrixHIGH(:,condvec,ILDidx));
                Results_STD_HIGH_R(ILDidx,condvec) = std(RespRMatrixHIGH(:,condvec,ILDidx));

                if kv.model == 2
                    %LSO Rate Difference
                    Results_Ave_LOW_LSO(ILDidx,condvec) = mean(ResultMatrixvecLOW_LSO(:,condvec,ILDidx));
                    Results_STD_LOW_LSO(ILDidx,condvec) = std(ResultMatrixvecLOW_LSO(:,condvec,ILDidx));
                    Results_Ave_MID_LSO(ILDidx,condvec) = mean(ResultMatrixvecMID_LSO(:,condvec,ILDidx));
                    Results_STD_MID_LSO(ILDidx,condvec) = std(ResultMatrixvecMID_LSO(:,condvec,ILDidx));
                    Results_Ave_HIGH_LSO(ILDidx,condvec) = mean(ResultMatrixvecHIGH_LSO(:,condvec,ILDidx));
                    Results_STD_HIGH_LSO(ILDidx,condvec) = std(ResultMatrixvecHIGH_LSO(:,condvec,ILDidx));           
                end

            end
        end

        %ResultMatrix= ResultMatrix';

        %Plotting (as in Laback, 2023).
        if i_model == 1
            LSOmode = 0;
        else
            LSOmode = 1;
        end
        %Plotting_MOC_LSO
        local_plot(LSOmode, Results_Ave_LOW, Results_Ave_MID, Results_Ave_HIGH,...
            Results_STD_LOW, Results_STD_MID, Results_STD_HIGH, ILD_vector, ILD,...
            Results_Ave_LOW_LSO, Results_STD_LOW_LSO, Results_Ave_MID_LSO, Results_STD_MID_LSO,...
            Results_Ave_HIGH_LSO, Results_STD_HIGH_LSO,i_model,flags);

    end % i_exp (1=smalt2014, 2=zilany2014)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% local functions
function local_plot(LSOmode, Results_Ave_LOW, Results_Ave_MID, Results_Ave_HIGH,...
    Results_STD_LOW, Results_STD_MID, Results_STD_HIGH, ILDvec, ILD,...
    Results_Ave_LOW_LSO, Results_STD_LOW_LSO, Results_Ave_MID_LSO, Results_STD_MID_LSO,...
    Results_Ave_HIGH_LSO, Results_STD_HIGH_LSO,jj,flags)

if LSOmode == 0
    W=[83 17 0]; %Weighting Coefficients for Low-, Mid-, and High-SR fibers of AN model: these can be optimized using the script  
    W_LSO=[80 20 0]; %Weighting Coefficients for Low, Mid, and High-SR fibers of LSO model; %[70 20 10]

    %Calculate Weighted Mean
    Results_Ave_WMean = ( Results_Ave_LOW*W(1) + Results_Ave_MID*W(2) + Results_Ave_HIGH*W(3) ) / sum(W);
    Results_STD_WMean = ( Results_STD_LOW*W(1) + Results_STD_MID*W(2) + Results_STD_HIGH*W(3) ) / sum(W);
    NrILD=length(Results_Ave_WMean(:,1));
    %ILDvec = 0:10/(NrILD-1):10;
    %ILDvec=0:2:10;
    %NrILD=6;

    %MOC on=====================================================
    figure(1); %Low-SR fibers
    subplot(4,4,1);
    hold on;
    errorbar(ILDvec,Results_Ave_LOW(:,2)', Results_STD_LOW(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_LOW(:,4)', Results_STD_LOW(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_LOW(:,6)', Results_STD_LOW(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_LOW(:,8)', Results_STD_LOW(:,8)', 'c'); %Contra
    axis([0 ILD -10 70]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,60,'LOW-SR fibers')
    legend('NoP','Diotic','Ipsi','Contra');
    hold off;

    %Mid-SR fibers
    subplot(4,4,2);
    hold on;
    errorbar(ILDvec,Results_Ave_MID(:,2)', Results_STD_MID(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_MID(:,4)', Results_STD_MID(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_MID(:,6)', Results_STD_MID(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_MID(:,8)', Results_STD_MID(:,8)', 'c'); %Contra
    axis([0 ILD -10 160]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,140,'MID-SR fibers')
    hold off;

    %High-SR fibers
    subplot(4,4,3);
    hold on;
    errorbar(ILDvec,Results_Ave_HIGH(:,2)', Results_STD_HIGH(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_HIGH(:,4)', Results_STD_HIGH(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_HIGH(:,6)', Results_STD_HIGH(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_HIGH(:,8)', Results_STD_HIGH(:,8)', 'c'); %Contra
    axis([0 ILD 0 160]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,140,'HIGH-SR fibers')
    hold off;

    %Weighted Mean
    subplot(4,4,4);
    hold on;
    errorbar(ILDvec,Results_Ave_WMean(:,2)', Results_STD_WMean(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_WMean(:,4)', Results_STD_WMean(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_WMean(:,6)', Results_STD_WMean(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_WMean(:,8)', Results_STD_WMean(:,8)', 'c'); %Contra
    axis([0 ILD 0 70]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,60,'Weighted Mean')
    hold off;

    %MOC off=====================================================
    figure(1); %Low-SR fibers
    subplot(4,4,5);

    hold on;
    errorbar(ILDvec,Results_Ave_LOW(:,1)', Results_STD_LOW(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_LOW(:,3)', Results_STD_LOW(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_LOW(:,5)', Results_STD_LOW(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_LOW(:,7)', Results_STD_LOW(:,7)', 'c'); %Contra
    axis([0 ILD -10 70]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,60,'LOW-SR fibers')
    legend('NoP','Diotic','Ipsi','Contra');
    hold off;

    %Mid-SR fibers
    subplot(4,4,6);

    hold on;
    errorbar(ILDvec,Results_Ave_MID(:,1)', Results_STD_MID(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_MID(:,3)', Results_STD_MID(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_MID(:,5)', Results_STD_MID(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_MID(:,7)', Results_STD_MID(:,7)', 'c'); %Contra
    axis([0 ILD -10 160]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,140,'MID-SR fibers')
    hold off;

    %High-SR fibers
    subplot(4,4,7);

    hold on;
    errorbar(ILDvec,Results_Ave_HIGH(:,1)', Results_STD_HIGH(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_HIGH(:,3)', Results_STD_HIGH(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_HIGH(:,5)', Results_STD_HIGH(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_HIGH(:,7)', Results_STD_HIGH(:,7)', 'c'); %Contra
    axis([0 ILD 0 160]); xlabel('Target ILD (dB)'); ylabel('Response Azimuth (deg)'); text(2,140,'HIGH-SR fibers')
    hold off;

    %Weighted Mean
    subplot(4,4,8);

    hold on;
    errorbar(ILDvec,Results_Ave_WMean(:,1)', Results_STD_WMean(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_WMean(:,3)', Results_STD_WMean(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_WMean(:,5)', Results_STD_WMean(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_WMean(:,7)', Results_STD_WMean(:,7)', 'c'); %Contra
    axis([0 ILD 0 70]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(1,60,'Weighted Mean')
    hold off;
end

if LSOmode == 1
    %LSO MODEL
    %Calculate Weighted Mean
    %W=[33.3 33.3 33.3]; %Libermann
    W_LSO=[80 20 0]; %Weighting Coefficients for Low, Mid, and High-SR fibers of LSO model; %[70 20 10]

    
    Results_Ave_WMean_LSO = ( Results_Ave_LOW_LSO*W_LSO(1) + Results_Ave_MID_LSO*W_LSO(2) + Results_Ave_HIGH_LSO*W_LSO(3) ) / sum(W_LSO);
    Results_STD_WMean_LSO = ( Results_STD_LOW_LSO*W_LSO(1) + Results_STD_MID_LSO*W_LSO(2) + Results_STD_HIGH_LSO*W_LSO(3) ) / sum(W_LSO);
    %ILDvec=0:5:10;
    %ILDvec=0:2:10;

    %MOC on=====================================================
    figure(1); %Low-SR fibers
    subplot(4,4,9);
    hold on;
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,2)', Results_STD_LOW_LSO(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,4)', Results_STD_LOW_LSO(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,6)', Results_STD_LOW_LSO(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,8)', Results_STD_LOW_LSO(:,8)', 'c'); %Contra
    axis([0 ILD -30 30]); xlabel('Target ILD (dB)'); ylabel('ILD_LSO (spike rate difference in Hz)'); text(2,20,'LOW-SR fibers')
    legend('NoP','Diotic','Ipsi','Contra');
    hold off;

    %Mid-SR fibers
    subplot(4,4,10);
    hold on;
    errorbar(ILDvec,Results_Ave_MID_LSO(:,2)', Results_STD_MID_LSO(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_MID_LSO(:,4)', Results_STD_MID_LSO(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_MID_LSO(:,6)', Results_STD_MID_LSO(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_MID_LSO(:,8)', Results_STD_MID_LSO(:,8)', 'c'); %Contra
    axis([0 ILD -60 160]); xlabel('Target ILD (dB)'); ylabel('ILD_LSO (spike rate difference in Hz)'); text(2,100,'MID-SR fibers')
    hold off;

    %High-SR fibers
    subplot(4,4,11)
    hold on;
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,2)', Results_STD_HIGH_LSO(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,4)', Results_STD_HIGH_LSO(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,6)', Results_STD_HIGH_LSO(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,8)', Results_STD_HIGH_LSO(:,8)', 'c'); %Contra
    axis([0 ILD -30 80]); xlabel('Target ILD (dB)'); ylabel('ILD_LSO (spike rate difference in Hz)'); text(2,35,'HIGH-SR fibers')
    hold off;

    %Weighted Mean
    subplot(4,4,12)
    hold on;
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,2)', Results_STD_WMean_LSO(:,2)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,4)', Results_STD_WMean_LSO(:,4)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,6)', Results_STD_WMean_LSO(:,6)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,8)', Results_STD_WMean_LSO(:,8)', 'c'); %Contra
    axis([0 ILD -30 80]); xlabel('Target ILD (dB)'); ylabel('ILD_LSO (spike rate difference in Hz)'); text(2,35,'Weighted Mean')
    hold off;

    %===%MOC off=====================================================
    figure(1); %Low-SR fibers
    subplot(4,4,13)
    hold on;
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,1)', Results_STD_LOW_LSO(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,3)', Results_STD_LOW_LSO(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,5)', Results_STD_LOW_LSO(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_LOW_LSO(:,7)', Results_STD_LOW_LSO(:,7)', 'c'); %Contra
    axis([0 ILD -10 40]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,35,'LOW-SR fibers')
    legend('NoP','Diotic','Ipsi','Contra');
    hold off;

    %Mid-SR fibers
    subplot(4,4,14)
    hold on;
    errorbar(ILDvec,Results_Ave_MID_LSO(:,1)', Results_STD_MID_LSO(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_MID_LSO(:,3)', Results_STD_MID_LSO(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_MID_LSO(:,5)', Results_STD_MID_LSO(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_MID_LSO(:,7)', Results_STD_MID_LSO(:,7)', 'c'); %Contra
    axis([0 ILD -10 40]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,35,'MID-SR fibers')
    hold off;

    %High-SR fibers
    subplot(4,4,15)
    hold on;
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,1)', Results_STD_HIGH_LSO(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,3)', Results_STD_HIGH_LSO(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,5)', Results_STD_HIGH_LSO(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_HIGH_LSO(:,7)', Results_STD_HIGH_LSO(:,7)', 'c'); %Contra
    axis([0 ILD 0 40]); xlabel('Target ILD (dB)'); ylabel('Response Azimuth (deg)'); text(2,35,'HIGH-SR fibers')
    hold off;

    %Weighted Mean
    subplot(4,4,16)
    hold on;
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,1)', Results_STD_WMean_LSO(:,1)', 'b'); %NoPrec
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,3)', Results_STD_WMean_LSO(:,3)', 'r'); %Diotic
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,5)', Results_STD_WMean_LSO(:,5)', 'g'); %Ipsi
    errorbar(ILDvec,Results_Ave_WMean_LSO(:,7)', Results_STD_WMean_LSO(:,7)', 'c'); %Contra
    axis([0 ILD 0 40]); xlabel('Target ILD (dB)'); ylabel('ILD_INT (spike rate difference in Hz)'); text(2,35,'Weighted Mean')
    hold off;
end

