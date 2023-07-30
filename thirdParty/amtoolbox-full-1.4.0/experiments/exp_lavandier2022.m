function varargout = exp_lavandier2022(varargin)
%EXP_LAVANDIER2022 Experiments of Lavandier (2022)
%
%   Usage: [] = exp_lavandier2022(flag) 
%
%
%   exp_lavandier(flag) reproduces figures of the study from 
%   Lavandier et al.(2022).
%
%   The following flags can be specified
%
%     'fig1' calculate prediction of exp 1 from Lavandier et al 2012, 
%            using the model lavandier2022 and ear signals (not BRIRs)
%
%     'fig2' calculate prediction of exp 1 from Lavandier et al. 2012, 
%            using model jelfs2011 and BRIRs (not ear signals)
%
%     'fig3' calculate prediction of exp 4 from Collin & Lavandier 2013 
%            with model vicente2020nh
%
%     'fig4' calculate predictions of experiment 3 of Lavandier & Culling 2008, 
%            with model leclere2015
%
%     'fig5' calculate predictions of the broadband conditions measured by 
%            Rana & Buchholz(2018) with the model vicente2020
%
%     'fig6' calculate predictions of exp 1 from Deroche et al (2014) with 
%            model prudhomme2020
%
%   Examples:
%   ---------
%
%   To display results for Fig.1 from Lavandier et al. (2022) use :
%
%     exp_lavandier2022('fig1');
%
%   To display results for Fig.2 from Lavandier et al. (2022) use :
%
%     exp_lavandier2022('fig2');
%
%   To display results for Fig.3 from Lavandier et al. (2022) use :
%
%     exp_lavandier2022('fig3');
%
%   To display results for Fig.4 from Lavandier et al. (2022) use :
%
%     exp_lavandier2022('fig4');
%
%   To display results for Fig.5 from Lavandier et al. (2022) use :
%
%     exp_lavandier2022('fig5');
%
%   To display results for Fig.6 from Lavandier et al. (2022) use :
%
%     exp_lavandier2022('fig6');
%
%   See also: lavandier2022 vicente2020nh vicente2020 prudhomme2020 leclere2015
%   jelfs2011
%
%   References:
%     Rana and Buchholz. Effect of audibility on better-ear glimpsing as a
%     function of frequency in normal-hearing and hearing-impaired listeners.
%     J. Acoust. Soc. Am., 143(4):2195--2206, 2018.
%     
%     M. Lavandier, T. Vicente, and L. Prud'homme. A series of snr-based
%     speech intelligibility models in the auditory modeling toolbox. Acta
%     Acustica, 2022.
%     
%     M. Lavandier, S. Jelfs, J. Culling, A. Watkins, A. Raimond, and
%     S. Makin. Binaural prediction of speech intelligibility in reverberant
%     rooms with multiple noise sources. J. Acoust. Soc. Am.,
%     131(1):218--231, 2012.
%     
%     M. Lavandier and J. Culling. Speech segregation in rooms: Monaural,
%     binaural and interacting effects of reverberation on target and
%     interferer. J. Acoust. Soc. Am., 123(4):2237--2248, 2008.
%     
%     T. Leclère, M. Lavandier, and J. Culling. Speech intelligibility
%     prediction in reverberation: Towards an integrated model of speech
%     transmission, spatial unmasking and binaural de-reverberation. J.
%     Acoust. Soc. Am., 137(6):3335--3345, 2015.
%     
%     M. Deroche, J. Culling, M. Chatterjee, and C. Limb. Speech recognition
%     against harmonic and inharmonic complexes: Spectral dips and
%     periodicity. J. Acoust. Soc. Am., 135(5):2873--2884, 2014.
%     
%     L. Prud'homme, M. Lavandier, and V. Best. A harmonic-cancellation-based
%     model to predict speech intelligibility against a harmonic masker. J.
%     Acoust. Soc. Am., 148(5):3246--3254, 2020.
%     
%     B. Collin and M. Lavandier. Binaural speech intelligibility in rooms
%     with variations in spatial location of sources and modulation depth of
%     noise interferers. J. Acoust. Soc. Am., 134(2):1146--1159, 2013.
%     
%     T. Vicente, M. Lavandier, and J. Buchholz. A binaural model
%     implementing an internal noise to predict the effect of hearing
%     impairment on speech intelligibility in non-stationary noises. J.
%     Acoust. Soc. Am., 148(5):3305--3317, 2020.
%     
%     T. Vicente and M. Lavandier. Further validation of a binaural model
%     predicting speech intelligibility against envelope-modulated noises.
%     Hearing Research, 390(107937), 2020.
%     
%     S. Jelfs, J. Culling, and M. Lavandier. Revision and validation of a
%     binaural model for speech intelligibility in noise. Hearing Research,
%     2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_lavandier2022.php


%   #Author: Matthieu Lavandier (2022)
%   #Author: Clara Hollomey (2022): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.import={'amt_cache'};
definput.flags.type = {'missingflag', 'fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'fig6'};


[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},...
             definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.', ...
      upper(mfilename),flagnames);
end


if flags.do_fig1
    %Prediction of exp 1 from Lavandier et al 2012, using model lavandier2022 and ear signals (not BRIRs)
    
    %%Tested conditions:
    %1=sl (short distance, left), 2=sf (short distance, front), 3=sr- (short distance, right),
    %4=ll (long distance, left), 5=lf, 6=lr, 
    %7=slp (short distance, left, processed=SEIR instead of BRIR), 8=sfp, 9=srp, 10=llp, 11=lfp, 12=lrp

    %predictions with all target sentences AVERAGED as the target input but using just one SSN for
    %the interferer (the stats change very little by concatenating more noises)

    %load data
    filepath1 = amt_load('lavandier2012', 'SRT.txt');
    SRT=dlmread(filepath1);
    filepath2 = amt_load('lavandier2012', 'StdError.txt');
    StdError=dlmread(filepath2);

    %load masker signals
    %keep only 4.2 second of the maskers, to avoid the non-stationary end of
    %the BRIR stimuli (reverberant decrease, non representative of the long-term
    %stimuli). Not really required for the SEIR that are very short (negligeble effect
    %on the long term spectrum)
    Fs=zeros(12,1);

    [Nsl, Fs(1)]= amt_load('lavandier2012', 'masker_n1sl.wav');
    [Nsf, Fs(2)]= amt_load('lavandier2012', 'masker_n1sf.wav');
    [Nsr, Fs(3)]= amt_load('lavandier2012', 'masker_n1sr.wav');
    [Nll, Fs(4)]= amt_load('lavandier2012', 'masker_n1ll.wav');
    [Nlf, Fs(5)]= amt_load('lavandier2012', 'masker_n1lf.wav');
    [Nlr, Fs(6)]= amt_load('lavandier2012', 'masker_n1lr.wav');
    [Nslp, Fs(7)]= amt_load('lavandier2012', 'masker_n1slp.wav');
    [Nsfp, Fs(8)]= amt_load('lavandier2012', 'masker_n1sfp.wav');
    [Nsrp, Fs(9)]= amt_load('lavandier2012', 'masker_n1srp.wav');
    [Nllp, Fs(10)]= amt_load('lavandier2012', 'masker_n1llp.wav');
    [Nlfp, Fs(11)]= amt_load('lavandier2012', 'masker_n1lfp.wav');
    [Nlrp, Fs(12)]= amt_load('lavandier2012', 'masker_n1lrp.wav');

    x = amt_load('lavandier2012', 'target_sr.mat');
    sig_sideclose = x.sig_sideclose;
    Fs_sideclose = x.Fs_sideclose;
    duration_sideclose = x.duration_sideclose;

    y = amt_load('lavandier2012', 'target_srp.mat');
    sig_sidecloseP = y.sig_sidecloseP;
    Fs_sidecloseP = y.Fs_sidecloseP;
    duration_sidecloseP = y.duration_sidecloseP;

    %averaged target signals used as model inputs
    %(all sentences limited to the duration of the shortest sentences)
    duration=min([min(duration_sideclose), min(duration_sidecloseP)]); 
    left_Tsr=zeros(duration,length(sig_sideclose)); right_Tsr=zeros(duration,length(sig_sideclose));
    for i=1:length(sig_sideclose)
        sig=cell2mat(sig_sideclose(i));
        left_Tsr(:,i)=sig(1:duration,1)*sqrt(length(sig_sideclose));
        right_Tsr(:,i)=sig(1:duration,2)*sqrt(length(sig_sideclose));
        clear sig
    end
    
    %average of all sentences limited to the duration of the shortest sentence
    Tsr=[mean(left_Tsr,2), mean(right_Tsr,2)];   
    left_Tsrp=zeros(duration,length(sig_sidecloseP)); right_Tsrp=zeros(duration,length(sig_sidecloseP));
    for i=1:length(sig_sidecloseP)
        sig=cell2mat(sig_sidecloseP(i));
        left_Tsrp(:,i)=sig(1:duration,1)*sqrt(length(sig_sidecloseP));
        right_Tsrp(:,i)=sig(1:duration,2)*sqrt(length(sig_sidecloseP));
        clear sig
    end
    %average of all sentences limited to the duration of the shortest sentence
    Tsrp=[mean(left_Tsrp,2), mean(right_Tsrp,2)];
    
    %RMS equalization. In Lavandier et al 2012, the mean left-right rms was equalized across condition, target/masker
    %equalize processed target with processed noise (no ITD) and unprocessed target with unprocessed noise
    noise_levels=zeros(1,6);
    noise_levels(1)=local_meanrms(Nsl); noise_levels(2)=local_meanrms(Nsf); noise_levels(3)=local_meanrms(Nsr); noise_levels(4)=local_meanrms(Nll);
    noise_levels(5)=local_meanrms(Nlf); noise_levels(6)=local_meanrms(Nlr);
    noisep_levels=zeros(1,6);
    noisep_levels(1)=local_meanrms(Nslp); noisep_levels(2)=local_meanrms(Nsfp); noisep_levels(3)=local_meanrms(Nsrp); noisep_levels(4)=local_meanrms(Nllp);
    noisep_levels(5)=local_meanrms(Nlfp); noisep_levels(6)=local_meanrms(Nlrp); 
    Tsr=Tsr*10^((mean(noise_levels)-local_meanrms(Tsr))/20);
    Tsrp=Tsrp*10^((mean(noisep_levels)-local_meanrms(Tsrp))/20);

    %compute prediction
    binauralSNR(1)=lavandier2022(Tsr,Nsl,Fs(1));
    binauralSNR(2)=lavandier2022(Tsr,Nsf,Fs(1));
    binauralSNR(3)=lavandier2022(Tsr,Nsr,Fs(1));
    binauralSNR(4)=lavandier2022(Tsr,Nll,Fs(1));
    binauralSNR(5)=lavandier2022(Tsr,Nlf,Fs(1));
    binauralSNR(6)=lavandier2022(Tsr,Nlr,Fs(1));
    binauralSNR(7)=lavandier2022(Tsrp,Nslp,Fs(1));
    binauralSNR(8)=lavandier2022(Tsrp,Nsfp,Fs(1));
    binauralSNR(9)=lavandier2022(Tsrp,Nsrp,Fs(1));
    binauralSNR(10)=lavandier2022(Tsrp,Nllp,Fs(1));
    binauralSNR(11)=lavandier2022(Tsrp,Nlfp,Fs(1));
    binauralSNR(12)=lavandier2022(Tsrp,Nlrp,Fs(1));

    %transformbinauralSNR into predicted SRT
    %to compare data and prediction, mean of predictions is set to mean of SRT
    predictSRT= mean(SRT)-(binauralSNR-mean(binauralSNR));       
    %compare data and prediction
    Corr_Pearson = corr(SRT, predictSRT','type','Pearson');
    Corr_Spearman = corr(SRT, predictSRT','type','Spearman');
    MeanError = mean(abs(SRT-predictSRT'));
    LargestError = max(abs(SRT-predictSRT'));
    RMSError = sqrt(mean((SRT-predictSRT').^2));

    %Plots
    fig1=figure; errorbar((1:1:6),SRT(1:6),StdError(1:6), 'bo', 'LineWidth', 2.5),
    xlabel('SSN position (target @ right/near)'), ylabel('SRT (dB)'), set(gca, 'XTick', 1:6), set(gca, 'XTickLabel', {'left/near'; 'front/near' ; 'right/near' ; 'left/far'; 'front/far' ; 'right/far'}),
    ylim([-11 -4]); xlim([.5 6.5]); grid, hold on,
    text(0.7, -4.2, 'model: lavandier2022')
    errorbar((1:1:6),SRT(7:12),StdError(7:12), 'ro', 'LineWidth', 2.5)
    plot((1:1:6),predictSRT(1:6), 'b-', 'LineWidth', 2.5),
    plot((1:1:6),predictSRT(7:12), 'r-', 'LineWidth', 2.5),
    legend('data (BRIR)','data (SEIR)', 'model (BRIR)', 'model (SEIR)','Location', 'NorthEast')
    temp=['r= ' num2str(round(Corr_Pearson*100)/100)]; text(4.5, -9.7, temp ), temp=['MeanErr= ' num2str(round(MeanError*10)/10) ' dB']; text(4.5, -10.2, temp ), temp=['MaxErr= ' num2str(round(LargestError*10)/10) ' dB']; text(4.5, -10.7, temp ), clear temp
    
end

if flags.do_fig2
    %Prediction of exp 1 from Lavandier et al. 2012, using model jelfs2011 and BRIRs (not ear signals)
    
    %%Tested conditions:
    %1=sl (short distance, left), 2=sf (short distance, front), 3=sr- (short distance, right),
    %4=ll (long distance, left), 5=lf, 6=lr, 
    %7=slp (short distance, left, processed=SEIR instead of BRIR), 8=sfp, 9=srp, 10=llp, 11=lfp, 12=lrp

    %load data
    filepath1 = amt_load('lavandier2012', 'SRT.txt');
    SRT=dlmread(filepath1);
    filepath2 = amt_load('lavandier2012', 'StdError.txt');
    StdError=dlmread(filepath2);

    %load masker BRIRs
    Fs=zeros(12,1);

    [Nsl, Fs(1)]= amt_load('lavandier2012', 'brirs_sl.wav');
    [Nsf, Fs(2)]= amt_load('lavandier2012', 'brirs_sf.wav');
    [Nsr, Fs(3)]= amt_load('lavandier2012', 'brirs_sr.wav');
    [Nll, Fs(4)]= amt_load('lavandier2012', 'brirs_ll.wav');
    [Nlf, Fs(5)]= amt_load('lavandier2012', 'brirs_lf.wav');
    [Nlr, Fs(6)]= amt_load('lavandier2012', 'brirs_lr.wav');
    [Nslp, Fs(7)]= amt_load('lavandier2012', 'brirs_slp.wav');
    [Nsfp, Fs(8)]= amt_load('lavandier2012', 'brirs_sfp.wav');
    [Nsrp, Fs(9)]= amt_load('lavandier2012', 'brirs_srp.wav');
    [Nllp, Fs(10)]= amt_load('lavandier2012', 'brirs_llp.wav');
    [Nlfp, Fs(11)]= amt_load('lavandier2012', 'brirs_lfp.wav');
    [Nlrp, Fs(12)]= amt_load('lavandier2012', 'brirs_lrp.wav');

    % filtering of all BRIRs by the speech-spectrum filter used to create the noise stimuli 
    speech_spectrum_filter = [-4.97213795824791e-06;-5.91990499287931e-07;-8.96838230346475e-07;-7.53207928028132e-07;-1.10881467207946e-06;-9.42195356401498e-07;-1.34701951992611e-06;-1.14838780973514e-06;-1.60782531111181e-06;-1.37454082960176e-06;-1.87913667559769e-06;-1.59106241426343e-06;-2.14470992432325e-06;-1.79389826371335e-06;-2.37868403019093e-06;-1.95255347534840e-06;-2.56730595538102e-06;-2.05027185984363e-06;-2.69159659183060e-06;-2.03905187845521e-06;-2.69845577349770e-06;-1.91947469829756e-06;-2.59222156273609e-06;-1.66116649324977e-06;-2.35923926084070e-06;-1.30694149902411e-06;-2.06227764465439e-06;-8.98607140698005e-07;-1.74156264165504e-06;-4.96577740705106e-07;-1.48076890127413e-06;-1.89131796446418e-07;-1.34562446874043e-06;-2.66167110574145e-09;-1.37933045607497e-06;-4.08994687006725e-08;-1.63528102348209e-06;-2.23946727828661e-07;-2.07402172236471e-06;-6.86418673012668e-07;-2.76812897936907e-06;-1.28170449897880e-06;-3.66121844308509e-06;-2.21534014599456e-06;-4.83521307614865e-06;-3.27740440297930e-06;-6.23310143055278e-06;-4.68824191557360e-06;-7.86714645073516e-06;-6.19943784840871e-06;-9.74719387158984e-06;-8.04422688815976e-06;-1.18510579341091e-05;-1.00863753687008e-05;-1.43436000143993e-05;-1.25738288261346e-05;-1.72637646755902e-05;-1.55583547893912e-05;-2.07951434276765e-05;-1.91569379239809e-05;-2.49658369284589e-05;-2.33804348681588e-05;-2.97294609481469e-05;-2.81163574982202e-05;-3.49014189851005e-05;-3.30817747453693e-05;-4.01357101509348e-05;-3.78563345293514e-05;-4.49967046733946e-05;-4.21609947807156e-05;-4.92710241815075e-05;-4.56896414107177e-05;-5.27869597135577e-05;-4.84385309391655e-05;-5.53727368242107e-05;-5.02208677062299e-05;-5.70885822526179e-05;-5.09671081090346e-05;-5.77728533244226e-05;-5.09027122461703e-05;-5.77224382141139e-05;-5.03129376738798e-05;-5.77343744225800e-05;-4.99990164826158e-05;-5.82771608605981e-05;-5.08340817759745e-05;-6.01769424974918e-05;-5.28606069565285e-05;-6.36259210295975e-05;-5.65567206649575e-05;-6.83813268551603e-05;-6.14576274529100e-05;-7.42190532037057e-05;-6.66903506498784e-05;-7.98882028902881e-05;-7.12400506017730e-05;-8.39765125419945e-05;-7.36052024876699e-05;-8.55386242619716e-05;-7.29246894479729e-05;-8.41807996039279e-05;-6.96814458933659e-05;-8.06597963673994e-05;-6.51030204608105e-05;-7.66037774155848e-05;-6.03502958256286e-05;-7.27825463400222e-05;-5.58257561351638e-05;-6.89247826812789e-05;-5.08896991959773e-05;-6.41348160570487e-05;-4.43620847363491e-05;-5.71364471397828e-05;-3.50936279573943e-05;-4.73243890155572e-05;-2.32034999498865e-05;-3.57107892341446e-05;-1.07119876702200e-05;-2.50278062594589e-05;-1.05431081465213e-06;-1.92453007912263e-05;1.91679600902717e-06;-2.17281467485009e-05;-4.23987921749358e-06;-3.33562711603008e-05;-1.83406900760019e-05;-5.04889976582490e-05;-3.42213788826484e-05;-6.51490190648474e-05;-4.28442799602635e-05;-6.79179793223739e-05;-3.52301030943636e-05;-5.05654606968164e-05;-4.72333249490475e-06;-8.73047702043550e-06;5.02412294736132e-05;5.60838416276965e-05;0.000125900493003428;0.000137767725391313;0.000214448809856549;0.000227590193389915;0.000307041103951633;0.000318128237267956;0.000397741649067029;0.000404926977353171;0.000483527459437028;0.000486156553961337;0.000563361681997776;0.000561159511562437;0.000635892385616899;0.000627973699010909;0.000699920696206391;0.000686733983457089;0.000757268513552845;0.000739923561923206;0.000810209428891540;0.000789527839515358;0.000861107720993459;0.000839823391288519;0.000917050696443766;0.000899786071386188;0.000986958388239145;0.000976325361989439;0.00107458082493395;0.00106938590761274;0.00117518915794790;0.00116881751455367;0.00127227860502899;0.00125107879284769;0.00133613229263574;0.00128354446496815;0.00133683520834893;0.00124500622041523;0.00126522965729237;0.00113707792479545;0.00112940196413547;0.000970508321188390;0.000941256934311241;0.000757536618039012;0.000711402331944555;0.000504163792356849;0.000438275397755206;0.000199539412278682;0.000103299054899253;-0.000177630077814683;-0.000312542397296056;-0.000643609091639519;-0.000825906288810074;-0.00121542694978416;-0.00144696224015206;-0.00188780634198338;-0.00215205084532499;-0.00261905975639820;-0.00288073718547821;-0.00333431642502546;-0.00355887529440224;-0.00397714320570231;-0.00415017036721110;-0.00452080322429538;-0.00463182153180242;-0.00495035620406270;-0.00499616703018546;-0.00526099419221282;-0.00525003019720316;-0.00548633141443133;-0.00544595252722502;-0.00567658245563507;-0.00563712883740664;-0.00591288879513741;-0.00592574477195740;-0.00627183075994253;-0.00633022561669350;-0.00673343287780881;-0.00682275742292404;-0.00724709080532193;-0.00729048391804099;-0.00763439247384667;-0.00753090716898441;-0.00770306773483753;-0.00734320469200611;-0.00725770648568869;-0.00666713807731867;-0.00645426893606782;-0.00575022399425507;-0.00549881858751178;-0.00485049793496728;-0.00484966160729528;-0.00446466309949756;-0.00473765702918172;-0.00463294005021453;-0.00530872726812959;-0.00543307187035680;-0.00601134542375803;-0.00563309807330370;-0.00569157488644123;-0.00486574694514275;-0.00479240156710148;-0.00398568809032440;-0.00407510157674551;-0.00306554441340268;-0.00259233312681317;-0.000892437994480133;-0.000552843906916678;0.000961960002314299;0.00178868125658482;0.00466188322752714;0.00644186418503523;0.0106202047318220;0.0130523610860109;0.0175840351730585;0.0206472221761942;0.0302094258368015;0.0343338213860989;0.0430484861135483;0.0434842184185982;0.0665084943175316;0.0665084943175316;0.0434842184185982;0.0430484861135483;0.0343338213860989;0.0302094258368015;0.0206472221761942;0.0175840351730585;0.0130523610860109;0.0106202047318220;0.00644186418503523;0.00466188322752714;0.00178868125658482;0.000961960002314299;-0.000552843906916678;-0.000892437994480133;-0.00259233312681317;-0.00306554441340268;-0.00407510157674551;-0.00398568809032440;-0.00479240156710148;-0.00486574694514275;-0.00569157488644123;-0.00563309807330370;-0.00601134542375803;-0.00543307187035680;-0.00530872726812959;-0.00463294005021453;-0.00473765702918172;-0.00446466309949756;-0.00484966160729528;-0.00485049793496728;-0.00549881858751178;-0.00575022399425507;-0.00645426893606782;-0.00666713807731867;-0.00725770648568869;-0.00734320469200611;-0.00770306773483753;-0.00753090716898441;-0.00763439247384667;-0.00729048391804099;-0.00724709080532193;-0.00682275742292404;-0.00673343287780881;-0.00633022561669350;-0.00627183075994253;-0.00592574477195740;-0.00591288879513741;-0.00563712883740664;-0.00567658245563507;-0.00544595252722502;-0.00548633141443133;-0.00525003019720316;-0.00526099419221282;-0.00499616703018546;-0.00495035620406270;-0.00463182153180242;-0.00452080322429538;-0.00415017036721110;-0.00397714320570231;-0.00355887529440224;-0.00333431642502546;-0.00288073718547821;-0.00261905975639820;-0.00215205084532499;-0.00188780634198338;-0.00144696224015206;-0.00121542694978416;-0.000825906288810074;-0.000643609091639519;-0.000312542397296056;-0.000177630077814683;0.000103299054899253;0.000199539412278682;0.000438275397755206;0.000504163792356849;0.000711402331944555;0.000757536618039012;0.000941256934311241;0.000970508321188390;0.00112940196413547;0.00113707792479545;0.00126522965729237;0.00124500622041523;0.00133683520834893;0.00128354446496815;0.00133613229263574;0.00125107879284769;0.00127227860502899;0.00116881751455367;0.00117518915794790;0.00106938590761274;0.00107458082493395;0.000976325361989439;0.000986958388239145;0.000899786071386188;0.000917050696443766;0.000839823391288519;0.000861107720993459;0.000789527839515358;0.000810209428891540;0.000739923561923206;0.000757268513552845;0.000686733983457089;0.000699920696206391;0.000627973699010909;0.000635892385616899;0.000561159511562437;0.000563361681997776;0.000486156553961337;0.000483527459437028;0.000404926977353171;0.000397741649067029;0.000318128237267956;0.000307041103951633;0.000227590193389915;0.000214448809856549;0.000137767725391313;0.000125900493003428;5.60838416276965e-05;5.02412294736132e-05;-8.73047702043550e-06;-4.72333249490475e-06;-5.05654606968164e-05;-3.52301030943636e-05;-6.79179793223739e-05;-4.28442799602635e-05;-6.51490190648474e-05;-3.42213788826484e-05;-5.04889976582490e-05;-1.83406900760019e-05;-3.33562711603008e-05;-4.23987921749358e-06;-2.17281467485009e-05;1.91679600902717e-06;-1.92453007912263e-05;-1.05431081465213e-06;-2.50278062594589e-05;-1.07119876702200e-05;-3.57107892341446e-05;-2.32034999498865e-05;-4.73243890155572e-05;-3.50936279573943e-05;-5.71364471397828e-05;-4.43620847363491e-05;-6.41348160570487e-05;-5.08896991959773e-05;-6.89247826812789e-05;-5.58257561351638e-05;-7.27825463400222e-05;-6.03502958256286e-05;-7.66037774155848e-05;-6.51030204608105e-05;-8.06597963673994e-05;-6.96814458933659e-05;-8.41807996039279e-05;-7.29246894479729e-05;-8.55386242619716e-05;-7.36052024876699e-05;-8.39765125419945e-05;-7.12400506017730e-05;-7.98882028902881e-05;-6.66903506498784e-05;-7.42190532037057e-05;-6.14576274529100e-05;-6.83813268551603e-05;-5.65567206649575e-05;-6.36259210295975e-05;-5.28606069565285e-05;-6.01769424974918e-05;-5.08340817759745e-05;-5.82771608605981e-05;-4.99990164826158e-05;-5.77343744225800e-05;-5.03129376738798e-05;-5.77224382141139e-05;-5.09027122461703e-05;-5.77728533244226e-05;-5.09671081090346e-05;-5.70885822526179e-05;-5.02208677062299e-05;-5.53727368242107e-05;-4.84385309391655e-05;-5.27869597135577e-05;-4.56896414107177e-05;-4.92710241815075e-05;-4.21609947807156e-05;-4.49967046733946e-05;-3.78563345293514e-05;-4.01357101509348e-05;-3.30817747453693e-05;-3.49014189851005e-05;-2.81163574982202e-05;-2.97294609481469e-05;-2.33804348681588e-05;-2.49658369284589e-05;-1.91569379239809e-05;-2.07951434276765e-05;-1.55583547893912e-05;-1.72637646755902e-05;-1.25738288261346e-05;-1.43436000143993e-05;-1.00863753687008e-05;-1.18510579341091e-05;-8.04422688815976e-06;-9.74719387158984e-06;-6.19943784840871e-06;-7.86714645073516e-06;-4.68824191557360e-06;-6.23310143055278e-06;-3.27740440297930e-06;-4.83521307614865e-06;-2.21534014599456e-06;-3.66121844308509e-06;-1.28170449897880e-06;-2.76812897936907e-06;-6.86418673012668e-07;-2.07402172236471e-06;-2.23946727828661e-07;-1.63528102348209e-06;-4.08994687006725e-08;-1.37933045607497e-06;-2.66167110574145e-09;-1.34562446874043e-06;-1.89131796446418e-07;-1.48076890127413e-06;-4.96577740705106e-07;-1.74156264165504e-06;-8.98607140698005e-07;-2.06227764465439e-06;-1.30694149902411e-06;-2.35923926084070e-06;-1.66116649324977e-06;-2.59222156273609e-06;-1.91947469829756e-06;-2.69845577349770e-06;-2.03905187845521e-06;-2.69159659183060e-06;-2.05027185984363e-06;-2.56730595538102e-06;-1.95255347534840e-06;-2.37868403019093e-06;-1.79389826371335e-06;-2.14470992432325e-06;-1.59106241426343e-06;-1.87913667559769e-06;-1.37454082960176e-06;-1.60782531111181e-06;-1.14838780973514e-06;-1.34701951992611e-06;-9.42195356401498e-07;-1.10881467207946e-06;-7.53207928028132e-07;-8.96838230346475e-07;-5.91990499287931e-07;-4.97213795824791e-06;];
    Nsl=[conv(Nsl(:,1),speech_spectrum_filter), conv(Nsl(:,2),speech_spectrum_filter)];    
    Nsf=[conv(Nsf(:,1),speech_spectrum_filter), conv(Nsf(:,2),speech_spectrum_filter)];
    Nsr=[conv(Nsr(:,1),speech_spectrum_filter), conv(Nsr(:,2),speech_spectrum_filter)];
    Nll=[conv(Nll(:,1),speech_spectrum_filter), conv(Nll(:,2),speech_spectrum_filter)];    
    Nlf=[conv(Nlf(:,1),speech_spectrum_filter), conv(Nlf(:,2),speech_spectrum_filter)];
    Nlr=[conv(Nlr(:,1),speech_spectrum_filter), conv(Nlr(:,2),speech_spectrum_filter)];
    Nslp=[conv(Nslp(:,1),speech_spectrum_filter), conv(Nslp(:,2),speech_spectrum_filter)];    
    Nsfp=[conv(Nsfp(:,1),speech_spectrum_filter), conv(Nsfp(:,2),speech_spectrum_filter)];
    Nsrp=[conv(Nsrp(:,1),speech_spectrum_filter), conv(Nsrp(:,2),speech_spectrum_filter)];
    Nllp=[conv(Nllp(:,1),speech_spectrum_filter), conv(Nllp(:,2),speech_spectrum_filter)];    
    Nlfp=[conv(Nlfp(:,1),speech_spectrum_filter), conv(Nlfp(:,2),speech_spectrum_filter)];
    Nlrp=[conv(Nlrp(:,1),speech_spectrum_filter), conv(Nlrp(:,2),speech_spectrum_filter)];

    %energy level equalisation of the BRIR and SEIR, respectively 
    %(no need to equalize BRIR to SEIR because they are never directly compared)
    %In Lavandier et al 2012, the mean left-right rms of the signals was equalized across condition, target/masker
    %reference for equalisation: frontal source at short distance (Nsf & Nsfp)
    Nsl=Nsl*10^((local_mean_nrj(Nsf)-local_mean_nrj(Nsl))/20);
    Nsr=Nsr*10^((local_mean_nrj(Nsf)-local_mean_nrj(Nsr))/20);
    Nll=Nll*10^((local_mean_nrj(Nsf)-local_mean_nrj(Nll))/20);
    Nlf=Nlf*10^((local_mean_nrj(Nsf)-local_mean_nrj(Nlf))/20);
    Nlr=Nlr*10^((local_mean_nrj(Nsf)-local_mean_nrj(Nlr))/20);
    Nslp=Nslp*10^((local_mean_nrj(Nsfp)-local_mean_nrj(Nslp))/20);
    Nsrp=Nsrp*10^((local_mean_nrj(Nsfp)-local_mean_nrj(Nsrp))/20);
    Nllp=Nllp*10^((local_mean_nrj(Nsfp)-local_mean_nrj(Nllp))/20);
    Nlfp=Nlfp*10^((local_mean_nrj(Nsfp)-local_mean_nrj(Nlfp))/20);
    Nlrp=Nlrp*10^((local_mean_nrj(Nsfp)-local_mean_nrj(Nlrp))/20);

    %target BRIRs (short distance on the right)
    Tsr=Nsr;
    Tsrp=Nsrp;

    %compute predictions
    binauralSNR(1)=jelfs2011(Tsr,Nsl,Fs(1), 'single');
    binauralSNR(2)=jelfs2011(Tsr,Nsf,Fs(1), 'single');
    binauralSNR(3)=jelfs2011(Tsr,Nsr,Fs(1), 'single');
    binauralSNR(4)=jelfs2011(Tsr,Nll,Fs(1), 'single');
    binauralSNR(5)=jelfs2011(Tsr,Nlf,Fs(1), 'single');
    binauralSNR(6)=jelfs2011(Tsr,Nlr,Fs(1), 'single');
    binauralSNR(7)=jelfs2011(Tsrp,Nslp,Fs(1), 'single');
    binauralSNR(8)=jelfs2011(Tsrp,Nsfp,Fs(1), 'single');
    binauralSNR(9)=jelfs2011(Tsrp,Nsrp,Fs(1), 'single');
    binauralSNR(10)=jelfs2011(Tsrp,Nllp,Fs(1), 'single');
    binauralSNR(11)=jelfs2011(Tsrp,Nlfp,Fs(1), 'single');
    binauralSNR(12)=jelfs2011(Tsrp,Nlrp,Fs(1), 'single');


    %transformbinauralSNR into predicted SRT
    %to compare data and prediction, mean of predictions is set to mean of SRT
    predictSRT= mean(SRT)-(binauralSNR-mean(binauralSNR));  
    
    %compare data and prediction
    Corr_Pearson = corr(SRT, predictSRT','type','Pearson');
    Corr_Spearman = corr(SRT, predictSRT','type','Spearman');
    MeanError = mean(abs(SRT-predictSRT'));
    LargestError = max(abs(SRT-predictSRT'));
    RMSError = sqrt(mean((SRT-predictSRT').^2));

    %Plots
    fig1=figure; errorbar((1:1:6),SRT(1:6),StdError(1:6), 'bo', 'LineWidth', 2.5),
    xlabel('SSN position (target @ right/near)'), ylabel('SRT (dB)'), set(gca, 'XTick', 1:6), set(gca, 'XTickLabel', {'left/near'; 'front/near' ; 'right/near' ; 'left/far'; 'front/far' ; 'right/far'}),
    ylim([-11 -4]); xlim([.5 6.5]); grid, hold on,
    text(0.7, -4.2, 'model: jelfs2011')
    errorbar((1:1:6),SRT(7:12),StdError(7:12), 'ro', 'LineWidth', 2.5)
    plot((1:1:6),predictSRT(1:6), 'b-', 'LineWidth', 2.5),
    plot((1:1:6),predictSRT(7:12), 'r-', 'LineWidth', 2.5),
    legend('data (BRIR)','data (SEIR)', 'model (BRIR)', 'model (SEIR)','Location', 'NorthEast')
    temp=['r= ' num2str(round(Corr_Pearson*100)/100)]; text(4.5, -9.7, temp ), temp=['MeanErr= ' num2str(round(MeanError*10)/10) ' dB']; text(4.5, -10.2, temp ), temp=['MaxErr= ' num2str(round(LargestError*10)/10) ' dB']; text(4.5, -10.7, temp ), clear temp
        
end

if flags.do_fig3
    %Prediction of exp 4 from Collin & Lavandier 2013 with model vicente2020nh
    
    %TARGET WAS ALWAYS THE SAME IN ALL CONDITIONS (close, in front)
    % 8 interferers were tested (3 modulations in 2 configs, and 2 modulations in one additional config)
    % 1=C03	0	1v (in front, 1-voice modulated)
    % 2=C01	0	2v (in front, 2-voice modulated)
    % 3=C08	0	st (in front, stationary)
    % 4=C06	p25	1v (azimut +25, 1-voice modulated)
    % 5=C04	p25	2v (azimut +25, 2-voice modulated)
    % 6=C02	p25	st (azimut +25, stationary)
    % 7=C07	mp25	2v (azimuts + and + 25, 2 noise interferers 1-voice modulated)
    % 8=C05	mp25	st (azimuts + and + 25, 2 noise interferers stationary)
    
    %Prog parameters
    FS=48000;                %ATTENTION FS=48kHz in Collin13
    nb_interferer=36;        %number of masker excerpts used for the predictions
    cut_begin=0.15*FS;       %cut silence at beginning of target/masker: 150ms are ok
    duration_mask=3.5*FS;    %masker duration used for the predictions
    ref_level=-30.3090;      %reference level for the rms equalization = equalized stimuli level (used for the experiment) before manipulation of these stimuli
    x=amt_load('collin2013', 'data_COL13_EXP4.mat');
    SRT = x.SRT;
    StdError = x.StdError;
    
    y = amt_load('collin2013', 'target.mat');
    target = y.target;
    %RMS equalisation
    target=target*10^((ref_level-local_meanrms(target))/20);

    % 1=C03	0	1v
    x1=amt_load('collin2013', 'sig_mask_1v_0.mat');
    sig_mask_1v_0 = x1.sig_mask_1v_0;
    % 2=C01	0	2v
    x2=amt_load('collin2013', 'sig_mask_2v_0.mat');
    sig_mask_2v_0 = x2.sig_mask_2v_0;
    % 3=C08	0	st
    x3=amt_load('collin2013', 'sig_mask_st_0.mat');
    sig_mask_st_0 = x3.sig_mask_st_0;
    % 4=C06	p25	1v
    x4=amt_load('collin2013', 'sig_mask_1v_p25.mat');
    sig_mask_1v_p25 = x4.sig_mask_1v_p25;
    % 5=C04	p25	2v
    x5=amt_load('collin2013', 'sig_mask_2v_p25.mat');
    sig_mask_2v_p25 = x5.sig_mask_2v_p25;
    % 6=C02	p25	st
    x6=amt_load('collin2013', 'sig_mask_st_p25.mat');
    sig_mask_st_p25 = x6.sig_mask_st_p25;
    % 7=C07	mp25 2v
    x7=amt_load('collin2013', 'sig_mask_2v_mp25.mat');
    sig_mask_2v_mp25 = x7.sig_mask_2v_mp25;
    % 8=C05	mp25 st
    x8=amt_load('collin2013', 'sig_mask_st_mp25.mat');
    sig_mask_st_mp25 = x8.sig_mask_st_mp25;
    %compute predictions, need to average the prediction across masker excerpts
    binauralSNR_1v_0=zeros(1,nb_interferer); BE_1v_0=zeros(1,nb_interferer); BU_1v_0=zeros(1,nb_interferer);
    binauralSNR_2v_0=zeros(1,nb_interferer); BE_2v_0=zeros(1,nb_interferer); BU_2v_0=zeros(1,nb_interferer);
    binauralSNR_st_0=zeros(1,nb_interferer); BE_st_0=zeros(1,nb_interferer); BU_st_0=zeros(1,nb_interferer);
    binauralSNR_1v_p25=zeros(1,nb_interferer); BE_1v_p25=zeros(1,nb_interferer); BU_1v_p25=zeros(1,nb_interferer);
    binauralSNR_2v_p25=zeros(1,nb_interferer); BE_2v_p25=zeros(1,nb_interferer); BU_2v_p25=zeros(1,nb_interferer);
    binauralSNR_st_p25=zeros(1,nb_interferer); BE_st_p25=zeros(1,nb_interferer); BU_st_p25=zeros(1,nb_interferer);
    binauralSNR_2v_mp25=zeros(1,nb_interferer); BE_2v_mp25=zeros(1,nb_interferer); BU_2v_mp25=zeros(1,nb_interferer);
    binauralSNR_st_mp25=zeros(1,nb_interferer); BE_st_mp25=zeros(1,nb_interferer); BU_st_mp25=zeros(1,nb_interferer);
    for i=1:nb_interferer
        %100*i/nb_interferer         %indicator of position within the loop, in percent
    %store the mean predictions across time frames (for each masker excerpt)   
    [binauralSNR_1v_0(i), BE_1v_0(i), BU_1v_0(i)]= vicente2020nh(target,cell2mat(sig_mask_1v_0(i)),FS);
    [binauralSNR_2v_0(i), BE_2v_0(i), BU_2v_0(i)]= vicente2020nh(target,cell2mat(sig_mask_2v_0(i)),FS);
    [binauralSNR_st_0(i), BE_st_0(i), BU_st_0(i)]= vicente2020nh(target,cell2mat(sig_mask_st_0(i)),FS);
    [binauralSNR_1v_p25(i), BE_1v_p25(i), BU_1v_p25(i)]= vicente2020nh(target,cell2mat(sig_mask_1v_p25(i)),FS);
    [binauralSNR_2v_p25(i), BE_2v_p25(i), BU_2v_p25(i)]= vicente2020nh(target,cell2mat(sig_mask_2v_p25(i)),FS);
    [binauralSNR_st_p25(i), BE_st_p25(i), BU_st_p25(i)]= vicente2020nh(target,cell2mat(sig_mask_st_p25(i)),FS);
    [binauralSNR_2v_mp25(i), BE_2v_mp25(i), BU_2v_mp25(i)]= vicente2020nh(target,cell2mat(sig_mask_2v_mp25(i)),FS);
    [binauralSNR_st_mp25(i), BE_st_mp25(i), BU_st_mp25(i)]= vicente2020nh(target,cell2mat(sig_mask_st_mp25(i)),FS);
    end
    %averaging across masker excerpt 
    binauralSNR(1)=mean(binauralSNR_1v_0);
    binauralSNR(2)=mean(binauralSNR_2v_0);
    binauralSNR(3)=mean(binauralSNR_st_0);
    binauralSNR(4)=mean(binauralSNR_1v_p25);
    binauralSNR(5)=mean(binauralSNR_2v_p25);
    binauralSNR(6)=mean(binauralSNR_st_p25);
    binauralSNR(7)=mean(binauralSNR_2v_mp25);
    binauralSNR(8)=mean(binauralSNR_st_mp25);

    %transform binauralSNR into predicted SRT
    predictSRT= mean(SRT)-(binauralSNR-mean(binauralSNR));       %to compare data and prediction, mean of predictions is set to mean of SRT
    %compare data and prediction
    Corr_Pearson = corr(SRT, predictSRT','type','Pearson');
    Corr_Spearman = corr(SRT, predictSRT','type','Spearman');
    MeanError = mean(abs(SRT-predictSRT'));
    LargestError = max(abs(SRT-predictSRT'));
    RMSError = sqrt(mean((SRT-predictSRT').^2));

    %Plot
    fig1=figure;  
    errorbar(1:3,SRT(1:3),StdError(1:3),'bo','Markersize',9,'linewidth',2.5), grid, hold on
    errorbar(1:3,SRT(4:6),StdError(4:6),'bv','Markersize',9,'linewidth',2.5)
    errorbar(2:3,SRT(7:8),StdError(7:8),'rs','Markersize',9,'linewidth',2.5)
    xlim([0.7,3.2]), ylim([-8.5,-2])
    set(gca,'FontSize',12,'XTickLabel',{'1-voice mod.','2-voice mod.','stationary'},'XTick',1:3,'TickLength',[0,0])
    xlabel({'Type of interfering noise'}), ylabel('SRT (dB)')
    legend('Interferer @ 0�','Interferer @ +25�','Interferer @ +/-25�','location','northwest'),
    text(0.8, -3.7, 'model (lines): vicente2020nh')
    temp=['r= ' num2str(round(Corr_Pearson*100)/100)]; text(2.5, -7, temp ), temp=['MeanErr= ' num2str(round(MeanError*10)/10) ' dB']; text(2.5, -7.5, temp ), temp=['MaxErr= ' num2str(round(LargestError*10)/10) ' dB']; text(2.5, -8, temp ), clear temp
    plot(1:3,predictSRT(1:3), 'b', 'LineWidth', 2.5,'HandleVisibility','off'),
    plot(1:3,predictSRT(4:6), 'b', 'LineWidth', 2.5,'HandleVisibility','off'),
    plot(2:3,predictSRT(7:8), 'r', 'LineWidth', 2.5,'HandleVisibility','off'),
end

if flags.do_fig4
    % Predictions of experiment 3 of Lavandier & Culling 2008, with model leclere2015.m

    fs = 20000;

    %%%%----------------BRIR/SRT loading------------------------------------%%%%%
    x = amt_load('lavandier2008', 'SRTLavandier2008.mat');
    SRT = x.SRT;
    y = amt_load('lavandier2008', 'BRIR.mat');
    BRIR = y.BRIR;
    filtering = 'on'; %load filters used for rms equalization of the brirs
    filepath1 = amt_load('lavandier2008', 'shape.ir');
    fid = fopen(filepath1);

    [shape] = fread(fid,inf,'float');
    fclose(fid);
    filepath2 = amt_load('lavandier2008', 'feminized_shape.ir');
    fid = fopen(filepath2);

    [feminizedShape] = fread(fid,inf,'float');
    fclose(fid);

    NRJCalculationMode = 1; %%--- 0 : energy equalization in temporal domain ; 1 : energy equalization in spectral domain (between 20Hz & Fs/2)
    EqMode = 1; %%% 1: applied to filtered BRIR, 2 : applied to raw BRIR, 3 : applied to spectrum of filtered BRIR , 4: applied to spectrum of raw BRIR
    interfererFilter = 'feminizedShape';
    targetFilter = 'shape';
    coupeBas = 'off';
    f_lim = 20;

    %%%%%%%%-----------Conditions and measurements------------%%%%%%%%%%%%%%

    measuredSRT = cell2mat(SRT(2,:));
    errorBar = cell2mat(SRT(3,:));
    absorption = {'1', '07', '05', '02'};
    canal = {'left', 'right'};
    position = {'L', 'R'};


    %%%%%------------------Plot settings-----------------------%%%%%%%%%%%%%

    plotResults = 'on';
    saveResults = 'on';

    %%%%%-------------------Variable declarations----------------%%%%%%%%%%%

    NRJ = zeros(4,4);
    %%%%*********************************************************************

    %%

    %%%%%----Load IR and energy computations
    for i = 1:4%%% Absorptions
        eval(['interferer' cell2mat(absorption(i)) '= cell2mat(BRIR(2,i+1));'])
        eval(['target' cell2mat(absorption(i)) '= cell2mat(BRIR(3,i+1));'])
        
        if i~= 0
            if strcmp(coupeBas, 'on')
                [b,a] = butter(8,0.02,'high');
                eval(['target' cell2mat(absorption(i)) '= filter(b,a, target' cell2mat(absorption(i)) ');'])
            end
        end
        
        if strcmp(filtering, 'on')
            eval(['interferer' cell2mat(absorption(i)) 'Filtered = [conv(' interfererFilter ',interferer' cell2mat(absorption(i)) '(:,1)) conv(' interfererFilter ', interferer' cell2mat(absorption(i)) '(:,2))];'])
            eval(['target' cell2mat(absorption(i)) 'Filtered = [conv(' targetFilter ',target' cell2mat(absorption(i)) '(:,1)) conv(' targetFilter ', target' cell2mat(absorption(i)) '(:,2))];'])
        else
            eval(['interferer' cell2mat(absorption(i)) 'Filtered = interferer' cell2mat(absorption(i)) ';'])
            eval(['target' cell2mat(absorption(i)) 'Filtered = target' cell2mat(absorption(i)) ';'])
        end
        
        eval(['SpectrumInterferer' cell2mat(absorption(i)) '= fft(interferer' cell2mat(absorption(i)) 'Filtered, pow2(nextpow2(length(interferer' cell2mat(absorption(i)) '))));'])
        eval(['SpectrumTarget' cell2mat(absorption(i)) '= fft(target' cell2mat(absorption(i)) 'Filtered, pow2(nextpow2(length(target' cell2mat(absorption(i)) '))));'])
        
        Nint = length(eval(['SpectrumInterferer' cell2mat(absorption(i))]));
        Ntar = length(eval(['SpectrumTarget' cell2mat(absorption(i))]));
        
        k20int = ceil(f_lim*Nint/fs);
        k20tar = ceil(f_lim*Ntar/fs);
            
        switch NRJCalculationMode
            case 1

                eval(['energyInt' cell2mat(absorption(i)) '= local_energy(SpectrumInterferer' cell2mat(absorption(i)) '(k20int:Nint-k20int,:))/(length(SpectrumInterferer' cell2mat(absorption(i)) ')-2*k20int+1);'])
                eval(['energyTar' cell2mat(absorption(i)) '= local_energy(SpectrumTarget' cell2mat(absorption(i)) '(k20tar:Ntar-k20tar,:))/(length(SpectrumTarget' cell2mat(absorption(i)) ')-2*k20tar+1);'])
                
            case 0
                
                eval(['energyInt' cell2mat(absorption(i)) '= local_energy(interferer' cell2mat(absorption(i)) 'Filtered);'])%%%%----Calculs d'energie sur BRIRs filtr�es
                eval(['energyTar' cell2mat(absorption(i)) '= local_energy(target' cell2mat(absorption(i)) 'Filtered);'])
                
            otherwise
                amt_disp('Wrong EqMode chosen')
        end
        NRJ(i,1:4) = [eval(['energyInt' cell2mat(absorption(i))]) eval(['energyTar' cell2mat(absorption(i))])];
    end
    interferer1 = [target1(:,2) target1(:,1)];
    
    clear i, clear j, clear k
    %%
    averagedEnergyAcrossConditions = 1;
    delta = [];
   

for i = 1:4 %%%% Absorption
       
            eval(['alpha(i,:) = sqrt(averagedEnergyAcrossConditions ./ energyTar' cell2mat(absorption(i)) ');']);
            eval(['beta(i,:) = sqrt(averagedEnergyAcrossConditions ./ energyInt' cell2mat(absorption(i)) ');']);
        
    switch EqMode
        case 1 %%%%%%% BRIReq = alpha * BRIR_Filtered---------------------------------------------------------------------------------
            
            eval(['target' cell2mat(absorption(i)) 'Eq(:,1) = alpha(i,1) * target' cell2mat(absorption(i)) 'Filtered(:,1);'])
            eval(['target' cell2mat(absorption(i)) 'Eq(:,2) = alpha(i,2) * target' cell2mat(absorption(i)) 'Filtered(:,2);'])
            eval(['interferer' cell2mat(absorption(i)) 'Eq(:,1) = beta(i,1) * interferer' cell2mat(absorption(i)) 'Filtered(:,1);'])
            eval(['interferer' cell2mat(absorption(i)) 'Eq(:,2) = beta(i,2) * interferer' cell2mat(absorption(i)) 'Filtered(:,2);'])
            
        case 2 %%%%% BRIReq = alpha * BRIR-----------------------------------------------------------------------------------------
            
            eval(['target' cell2mat(absorption(i)) 'Eq(:,1) = alpha(i,1) * target' cell2mat(absorption(i)) '(:,1);'])
            eval(['target' cell2mat(absorption(i)) 'Eq(:,2) = alpha(i,2) * target' cell2mat(absorption(i)) '(:,2);'])
            eval(['interferer' cell2mat(absorption(i)) 'Eq(:,1) = beta(i,1) * interferer' cell2mat(absorption(i)) '(:,1);'])
            eval(['interferer' cell2mat(absorption(i)) 'Eq(:,2) = beta(i,2) * interferer' cell2mat(absorption(i)) '(:,2);'])
            
        case 3 %%%%% FFT(BRIReq) = alpha * FFT(BRIR_Filtered)--------------------------------------------------------------------------
            
            %%% TARGET
            for k =1:2
                eval(['klim = ceil(f_lim*length(SpectrumTarget' cell2mat(absorption(i)) ')/fs);'])
                for n =1:length(eval(['SpectrumTarget' cell2mat(absorption(i))]))
                    if klim <= n && n <= length(eval(['SpectrumTarget' cell2mat(absorption(i))]))-klim
                        eval(['SpectrumTarget' cell2mat(absorption(i)) 'eq(n,k) = alpha(i,k) * SpectrumTarget' cell2mat(absorption(i)) '(n,k);'])
                    else
%                         eval(['SpectrumTarget' cell2mat(absorption(i)) 'eq(n,k) = SpectrumTarget' cell2mat(absorption(i)) '(n,k);'])
                        eval(['SpectrumTarget' cell2mat(absorption(i)) 'eq(n,k) = 0;'])
                    end
                end
                
                
            %%% INTERFERER
                eval(['klim = ceil(f_lim*length(SpectrumInterferer' cell2mat(absorption(i)) ')/fs);'])
                for n =1:length(eval(['SpectrumInterferer' cell2mat(absorption(i))]))
                    if klim <= n && n <= length(eval(['SpectrumInterferer' cell2mat(absorption(i))]))-klim
                        eval(['SpectrumInterferer' cell2mat(absorption(i)) 'eq(n,k) = beta(i,k) * SpectrumInterferer' cell2mat(absorption(i)) '(n,k);'])
                    else
%                         eval(['SpectrumInterferer' cell2mat(absorption(i)) 'eq(n,k) = SpectrumInterferer' cell2mat(absorption(i)) '(n,k);'])
                    eval(['SpectrumInterferer' cell2mat(absorption(i)) 'eq(n,k) = 0;'])
                    end
                end
                
            end
            eval(['target' cell2mat(absorption(i)) 'Eq = real(ifft(SpectrumTarget' cell2mat(absorption(i)) 'eq, pow2(nextpow2(length(SpectrumTarget' cell2mat(absorption(i)) 'eq)))));'])
            eval(['interferer' cell2mat(absorption(i)) 'Eq = real(ifft(SpectrumInterferer' cell2mat(absorption(i)) 'eq, pow2(nextpow2(length(SpectrumInterferer' cell2mat(absorption(i)) 'eq)))));'])
            
        case 4 %%%%% FFT(BRIReq) = alpha * FFT(BRIR)--------------------------------------------------------------------------------------
            %%%% TARGET
            eval(['TARGET' cell2mat(absorption(i)) ' = fft(target' cell2mat(absorption(i)) ',pow2(nextpow2(length(target' cell2mat(absorption(i)) '))));'])
            eval(['klim = ceil(f_lim*length(TARGET' cell2mat(absorption(i)) ')/fs);'])
            for k = 1:2

                eval(['energyTARGET' cell2mat(absorption(i)) '(:,k) = energy(TARGET' cell2mat(absorption(i)) '(klim:length(TARGET' cell2mat(absorption(i)) ')-klim,k))/(length(TARGET' cell2mat(absorption(i)) ')-2*klim+1);'])

%                 eval(['alpha = sqrt(averagedEnergyAcrossConditions / energyTARGET' cell2mat(absorption(i)) '(:,k));'])
                eval(['target' cell2mat(absorption(i)) 'Eq(:,k) = alpha(i,k) * target' cell2mat(absorption(i)) 'Filtered(:,k);'])
                for n =1:length(eval(['TARGET' cell2mat(absorption(i))]))
                    if klim <= n && n <= length(eval(['TARGET' cell2mat(absorption(i))]))-klim
                        eval(['TARGET' cell2mat(absorption(i)) 'eq(n,k) = alpha(i,k) * TARGET' cell2mat(absorption(i)) '(n,k);'])
                    else
                        eval(['TARGET' cell2mat(absorption(i)) 'eq(n,k) = TARGET' cell2mat(absorption(i)) '(n,k);'])
                    end
                end
            end
            
            %%%% INTERFERER
            eval(['INTERFERER' cell2mat(absorption(i)) ' = fft(interferer' cell2mat(absorption(i)) ',pow2(nextpow2(length(interferer' cell2mat(absorption(i)) '))));'])
            eval(['klim = ceil(f_lim*length(INTERFERER' cell2mat(absorption(i)) ')/fs);'])
            for k = 1:2

                eval(['energyINTERFERER' cell2mat(absorption(i)) '(:,k) = local_energy(INTERFERER' cell2mat(absorption(i)) '(klim:length(INTERFERER' cell2mat(absorption(i)) ')-klim,k))/(length(INTERFERER' cell2mat(absorption(i)) ')-2*klim+1);'])

%                 eval(['alpha = sqrt(averagedEnergyAcrossConditions / energyINTERFERER' cell2mat(absorption(i)) '(:,k));'])
                eval(['interferer' cell2mat(absorption(i)) 'Eq(:,k) = alpha(i,k) * interferer' cell2mat(absorption(i)) 'Filtered(:,k);'])
                for n =1:length(eval(['INTERFERER' cell2mat(absorption(i))]))
                    if klim <= n && n <= length(eval(['INTERFERER' cell2mat(absorption(i))]))-klim
                        eval(['INTERFERER' cell2mat(absorption(i)) 'eq(n,k) = alpha(i,k) * INTERFERER' cell2mat(absorption(i)) '(n,k);'])
                    else
                        eval(['INTERFERER' cell2mat(absorption(i)) 'eq(n,k) = INTERFERER' cell2mat(absorption(i)) '(n,k);'])
                    end
                end
            end
            
            
            eval(['interferer' cell2mat(absorption(i)) 'Eq = real(ifft(INTERFERER' cell2mat(absorption(i)) 'eq, pow2(nextpow2(length(INTERFERER' cell2mat(absorption(i)) 'eq)))));'])
        otherwise
            
    end
    
end
%%    
    clear i, clear j, clear k

%APPLICATION OF THE MODEL in the different conditions    

    TIR = zeros(1,16);
    BMLD = [];
    betterEar = [];
    for i = 1:4
        for j = 1:4
            %%%----TIR computations
            eval(['[TIR(4*(i-1)+j), BMLD(:,4*(i-1)+j), betterEar(:,4*(i-1)+j)] = leclere2015(target' cell2mat(absorption(j)) 'Eq, interferer' cell2mat(absorption(i)) 'Eq, fs);'])
        end
    end
    
    TIR = -TIR;
    fittedTIR = TIR - mean(TIR) + mean(mean(measuredSRT));%%%----TIR alinement
        
%%%%%--------------STATS RESULTS---------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    coefficientCorrelation = corr(measuredSRT', fittedTIR','type','Pearson');
    erreur = sum(abs(measuredSRT-fittedTIR));
    MeanErr = erreur / length(measuredSRT);
    RMSE = sqrt(mean((measuredSRT-fittedTIR).^2));
    largestError = max(abs(measuredSRT - fittedTIR));
    spearCoeff = corr(measuredSRT', fittedTIR', 'type', 'Spearman');
    
%%%%%--------------PLOTS--------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(plotResults, 'on')
        fig1=figure;
        errorbar(measuredSRT, errorBar, 'bo', 'LineWidth', 2.5),
        ylim([-10.5 -3.5]); xlim([.5 16.5]); grid, hold on,
        plot((1:1:4),fittedTIR(1:4), 'b-', 'LineWidth', 2.5)
        h=plot((5:1:8),fittedTIR(5:8), 'b-', 'LineWidth', 2.5); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h=plot((9:1:12),fittedTIR(9:12), 'b-', 'LineWidth', 2.5); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h=plot((13:1:16),fittedTIR(13:16), 'b-', 'LineWidth', 2.5); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        legend('data','model', 'Location', 'SouthEast')
        text(1, -5.2, 'model: leclere2015')
        temp=['r= ' num2str(round(coefficientCorrelation*100)/100)]; text(8.7, -9.2, temp ), temp=['MeanErr= ' num2str(round(MeanErr*10)/10) ' dB']; text(8.7, -9.7, temp ), temp=['MaxErr= ' num2str(round(largestError*10)/10) ' dB']; text(8.7, -10.2, temp ), clear temp
        h=line([4.5 4.5],[-3.5 -10.5], 'color', 'k'); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h=line([8.5 8.5],[-3.5 -10.5], 'color', 'k'); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h=line([12.5 12.5],[-3.5 -10.5], 'color', 'k'); h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        xlabel('Absorption coefficient used for the target'), ylabel('SRT (dB)'),
        set(gca, 'XTick', 1:16)
        set(gca, 'XTickLabel', {'1'; '0.7' ; '0.5' ; '0.2'})
        text(1, -3.8, 'Interferer'), text(5, -3.8, 'Interferer'), text(9, -3.8, 'Interferer'), text(13, -3.8, 'Interferer')
        text(1, -4.1, 'absorp. 1'), text(5, -4.1, 'absorp. 0.7'), text(9, -4.1, 'absorp. 0.5'), text(13, -4.1, 'absorp. 0.2')
                     
     end
end

if flags.do_fig5
    %predictions of the broadband conditions measured by Rana & Buchholz(2018) with the model vicente2020
    
    % Conditions:
    % target always in front (same in all conditions), not diotic because ear-specific amplification is applied on stimuli at each ear
    % maskers (noise-vocoded speech) spatially co-located: both distractors presented to both ears > not diotic because ear-specific amplification is applied on stimuli at each ear
    % maskers (noise-vocoded speech) spatially separated: one distrator send to each ear only (no cross-talk, "infinite ILD")
    % 4 noise level conditions: 0, 10, 20 and 30 dB SL 
    FS = 44100;

    %% ------------------------ MODEL PARAMETERS --------------------------- %%
    % Parameter for the better-ear SNR computation
    Ceiling = 20; % Maximum SNR allowed
    % Time-resolution parameters
    BE_duration = 0.024; % set the time frame duration to compute the SNR at the better-ear
    BU_duration = 0.300; % set the time frame duration to compute the binaural unmasking advantage. Take into account binaural sluggishness.
    HannWindowBE = hann(floor(BE_duration*FS)); % create a Hann window, that is the window used for the computation of the SNR at the better ear
    HannWindowBU = hann(floor(BU_duration*FS)); % create a Hann window, that is the window used for the computation of binaural unmasking advantage
    WindowOverlap = 0.5; % fifty percent of overlap
    % SET FREQUENCY SAMPLING of MODEL
    nerbs_step = 0.5; %Number of ERB between 2 gammatone filters
    nerbs = 1:nerbs_step:round(f2erbrate(FS/2));
    fc = zeros(1,length(nerbs));
    for i=1:length(fc)
        fc(i) = round(erbrate2f(nerbs(i)));
    end
    % SII weightings used in the model
    weightings = f2siiweightings(fc);
    % Parameter to set the internal noise Formula
    B = -10; % Parameter setting the internal noise floor of the listener
    Nlim = 83; % Parameter setting the external level at which the internal noise level starts to increase  
    ratio_OHC = 0.7; % set the percentage of the hearing threshold elevation due to the OHC loss (HLtotal = OHC_loss + IHC_loss)
    max_OHC = 57.6; % maximum allowed for the Outer Hair Cell loss, the loss above this value contribute to the Inner Hair Cell loss



    %% ----------------- EXPERIMENT-SPECIFIC METHODS & DATA ---------------- %%
    filepath1 = amt_load('rana2018', 'NH_names.txt');
    NH_names = importdata(filepath1);
    filepath2 = amt_load('rana2018', 'HI_names.txt');
    HI_names = importdata(filepath2);
    NBlistener_NH = length(NH_names); 
    NBlistener_HI = length(HI_names); 
    NBcdtns = 8; 
    % Load SRTs
    %Some HI listeners was not tested at all sensation levels to prevent 
    %them from loudness discomfort. This is show by a "NaN" in the data matrix.
    filepath3 = amt_load('rana2018', 'SRT_HI.txt');
    data_HI = importdata(filepath3);
    filepath4 = amt_load('rana2018', 'SRT_NH.txt');
    data_NH = importdata(filepath4);

    % Load audiograms
    f0 = [250 500 1000 2000 3000 4000 6000 8000]; % central frequency at which the audiograms have been measured
    filepath5 = amt_load('rana2018', 'audiog_R_HI.txt');
    audTonRdBHL_HI=dlmread(filepath5);
    filepath6 = amt_load('rana2018', 'audiog_L_HI.txt');
    audTonLdBHL_HI=dlmread(filepath6);
    filepath7 = amt_load('rana2018', 'audiog_R_NH.txt');
    audTonRdBHL_NH=dlmread(filepath7);
    filepath8 = amt_load('rana2018', 'audiog_L_NH.txt');
    audTonLdBHL_NH=dlmread(filepath8);
    % Load SRT in quiet (used to define the 0 dB SL)
    filepath9 = amt_load('rana2018', 'SRT_quiet_SPL_left_NH.txt');
    SRT_quiet_left_NH=dlmread(filepath9);
    filepath10 = amt_load('rana2018', 'SRT_quiet_SPL_right_NH.txt');
    SRT_quiet_right_NH=dlmread(filepath10);
    filepath11 = amt_load('rana2018', 'SRT_quiet_SPL_left_HI.txt');
    SRT_quiet_left_HI=dlmread(filepath11);
    filepath12 = amt_load('rana2018', 'SRT_quiet_SPL_right_HI.txt');
    SRT_quiet_right_HI=dlmread(filepath12);

    %% ---------------------------- VARIABLES ------------------------------ %%
    % Model Outputs
    BinauralRatio_NH = zeros(NBlistener_NH,NBcdtns);
    BE_SNR_NH = zeros(NBlistener_NH,NBcdtns);
    BUAdv_NH = zeros(NBlistener_NH,NBcdtns);
    BinauralRatio_HI = zeros(NBlistener_HI,NBcdtns);
    BE_SNR_HI = zeros(NBlistener_HI,NBcdtns);
    BUAdv_HI = zeros(NBlistener_HI,NBcdtns);
    % Audiogram related 
    aud_L_OHC_SPL_NH = zeros(NBlistener_NH,length(fc)); % hearing loss related to OHC loss in dB SPL 
    aud_R_OHC_SPL_NH = zeros(NBlistener_NH,length(fc));
    aud_L_IHC_HL_NH = zeros(NBlistener_NH,length(fc));  % hearing loss related to IHC loss in dB HL
    aud_R_IHC_HL_NH = zeros(NBlistener_NH,length(fc)); 
    aud_L_OHC_SPL_HI = zeros(NBlistener_HI,length(fc)); % hearing loss related to OHC loss in dB SPL 
    aud_R_OHC_SPL_HI = zeros(NBlistener_HI,length(fc));
    aud_L_IHC_HL_HI = zeros(NBlistener_HI,length(fc));  % hearing loss related to IHC loss in dB HL
    aud_R_IHC_HL_HI = zeros(NBlistener_HI,length(fc)); 
    % Internal Noise 
    InternalNoise_R_NH = cell(NBlistener_NH,1);
    InternalNoise_L_NH = cell(NBlistener_NH,1);
    InternalNoise_R_HI = cell(NBlistener_HI,1);
    InternalNoise_L_HI = cell(NBlistener_HI,1);

    %% ---------------------- Predictions HI listeners --------------------- %%
    [BinauralRatio_HI,BE_SNR_HI , BUAdv_HI] = amt_cache('get','HI_predictions',flags.cachemode);
    if isempty(BinauralRatio_HI)
        for listenerNB = 1:NBlistener_HI
            amt_disp(num2str(listenerNB))    %To see the progress
            CurrentListener = HI_names{listenerNB};  %listener code/name as a string
        %    ListenerPath = [RootPath '\signals\HI\' CurrentListener]; %Path to get the signals of the listener
            Participate = ~isnan(data_HI(listenerNB,:));

            % load input signals 
            %[MaskerColocSigs, MaskerSeparSigs, TargetSigs] = load_signals_RANA18a(RootPath,ListenerPath,[SRT_quiet_left_HI(listenerNB) SRT_quiet_right_HI(listenerNB)]);
            x = amt_load('rana2018', ['HI_', CurrentListener,'.mat']);
            MaskerColocSigs = x.MaskerColocSigs;
            MaskerSeparSigs = x.MaskerSeparSigs;
            TargetSigs = x.TargetSigs;

            % Internal noise implementation (external level is approximated by the known external noise level)
            Aud_L = audTonLdBHL_HI(listenerNB,:); %listener audiogram
            Aud_R = audTonRdBHL_HI(listenerNB,:);
            N = mean([SRT_quiet_left_HI(listenerNB) SRT_quiet_right_HI(listenerNB)]) + [0 10 20 30]; % SRT in quiet = 0 dB SL + [0 10 20 30] dB SL to get the level for each level condition
            [InternalNoise_L_HI{listenerNB}, InternalNoise_R_HI{listenerNB}, ~, ~, aud_L_OHC_SPL_HI(listenerNB,:), aud_R_OHC_SPL_HI(listenerNB,:), aud_L_IHC_HL_HI(listenerNB,:), aud_R_IHC_HL_HI(listenerNB,:), Gamma] = ...
                vicente2020_internalnoise(fc, f0, Aud_L, Aud_R, N, B, Nlim, ratio_OHC, max_OHC);

            % Apply the Vicente2020 model (model outputs is listener-dependent because it is based on the individual audiograms, and the signals can be listener dependent such as here)
            [BinauralRatio_HI(listenerNB,:), BE_SNR_HI(listenerNB,:), BUAdv_HI(listenerNB,:)] = ...
                local_predict_RANA18a_Broad_listener_Vicente20(MaskerColocSigs, MaskerSeparSigs, TargetSigs, InternalNoise_L_HI{listenerNB}, InternalNoise_R_HI{listenerNB}, Ceiling, HannWindowBE, HannWindowBU, WindowOverlap, weightings, FS, fc, Participate);

        end
        amt_cache('set','HI_predictions',BinauralRatio_HI,BE_SNR_HI , BUAdv_HI);
    end
    %% ---------------------- Predictions NH listeners --------------------- %%
     [BinauralRatio_NH,BE_SNR_NH , BUAdv_NH] = amt_cache('get','NH_predictions',flags.cachemode);
    if isempty(BinauralRatio_NH)   
    
        for listenerNB = 1:NBlistener_NH
            disp(num2str(listenerNB))    %To see the progress
            CurrentListener = NH_names{listenerNB};  %listener code/name as a string
            %ListenerPath = [RootPath '\signals\NH\' CurrentListener]; %Path to get the signals of the listener
            Participate = ~isnan(data_NH(listenerNB,:));
            %dir(['*',R])

            % load input signals 
            %[MaskerColocSigs, MaskerSeparSigs, TargetSigs] = load_signals_RANA18a(RootPath,ListenerPath,[SRT_quiet_left_NH(listenerNB) SRT_quiet_right_NH(listenerNB)]);
            %x = amt_load('rana2018', 'NH_target.mat'); 
            x = amt_load('rana2018', ['NH_', CurrentListener,'.mat']);
            MaskerColocSigs = x.MaskerColocSigs;
            MaskerSeparSigs = x.MaskerSeparSigs;
            TargetSigs = x.TargetSigs;



            % Internal noise implementation (external level is approximated by the known external noise level)
            Aud_L = audTonLdBHL_NH(listenerNB,:); %listener audiogram
            Aud_R = audTonRdBHL_NH(listenerNB,:);
            N = mean([SRT_quiet_left_NH(listenerNB) SRT_quiet_right_NH(listenerNB)]) + [0 10 20 30]; % SRT in quiet = 0 dB SL + [0 10 20 30] dB SL to get the level for each level condition
            [InternalNoise_L_NH{listenerNB}, InternalNoise_R_NH{listenerNB}, ~, ~, aud_L_OHC_SPL_NH(listenerNB,:), aud_R_OHC_SPL_NH(listenerNB,:), aud_L_IHC_HL_NH(listenerNB,:), aud_R_IHC_HL_NH(listenerNB,:), Gamma] = ...
                vicente2020_internalnoise(fc, f0, Aud_L, Aud_R, N, B, Nlim, ratio_OHC, max_OHC);

            % Apply the Vicente2020 model (model outputs is listener-dependent because it is based on the individual audiograms, and the signals can be listener dependent such as here)
            [BinauralRatio_NH(listenerNB,:), BE_SNR_NH(listenerNB,:), BUAdv_NH(listenerNB,:)] = ...
                local_predict_RANA18a_Broad_listener_Vicente20(MaskerColocSigs, MaskerSeparSigs, TargetSigs, InternalNoise_L_NH{listenerNB}, InternalNoise_R_NH{listenerNB}, Ceiling, HannWindowBE, HannWindowBU, WindowOverlap, weightings, FS, fc, Participate);
        end
        amt_cache('set','NH_predictions',BinauralRatio_NH,BE_SNR_NH , BUAdv_NH);
    end
    %% -------------- CONVERT BinauralSNR into predicted SRTs -------------- %%
    MappingReference = mean([mean(data_NH) nanmean(data_HI)]) - (- mean([mean(BinauralRatio_NH) nanmean(BinauralRatio_HI)]));
    PredictedSRT_HI = - BinauralRatio_HI + MappingReference;
    PredictedSRT_NH = - BinauralRatio_NH + MappingReference;

    %% --------------------- Average and Standard Error -------------------- %%

    AvPredictedSRT_HI = nanmean(PredictedSRT_HI);
    AvPredictedSRT_NH = mean(PredictedSRT_NH);
    PredictedStdErr_NH = std(PredictedSRT_NH,0,1) ./ sqrt(NBlistener_NH);
    PredictedStdErr_HI = std(PredictedSRT_HI,0,1) ./ sqrt(sum(~isnan(PredictedSRT_HI)));

    AvMeasuredSRT_HI = nanmean(data_HI);
    AvMeasuredSRT_NH = mean(data_NH);
    MeasuredStdErr_NH = std(data_NH,0,1) ./ sqrt(NBlistener_NH);
    MeasuredStdErr_HI = nanstd(data_HI,0,1) ./ sqrt(sum(~isnan(PredictedSRT_HI)));

    %% ----------------------- Performance statistics ---------------------- %%
    Corr_Pearson = corr([AvMeasuredSRT_NH AvMeasuredSRT_HI]', [AvPredictedSRT_NH AvPredictedSRT_HI]','type','Pearson');
    MeanError = mean(abs([AvMeasuredSRT_NH AvMeasuredSRT_HI]'- [AvPredictedSRT_NH AvPredictedSRT_HI]')); 
    MaximumError = max(abs([AvMeasuredSRT_NH AvMeasuredSRT_HI]'- [AvPredictedSRT_NH AvPredictedSRT_HI]'));

    %% -------------------------------- PLOT ------------------------------- %%
    fig1 = figure;
    xStart = 0.08; xlength = 0.85;yStart = 0.1; ylength = 0.8;xShiftBetweenPanel = 0.03;yShiftBetweenPanel = 0.02;
    NbRow = 1;HeightPanel = ylength/NbRow;
    NbPanel = 2; % one for each listener group
    NbColumn = round(NbPanel/NbRow);
    LengthPanel = xlength / NbColumn;
    PanelPos = zeros(NbPanel,4);
    rank = 1;
    for n = 1:NbRow
        for i = 1:NbColumn
            PanelPos(rank,:) = [xStart+(i-1)*(LengthPanel+xShiftBetweenPanel) yStart+(n-1)*(HeightPanel+yShiftBetweenPanel) LengthPanel HeightPanel];
            rank = rank + 1;
        end
    end
    YMax = max([AvPredictedSRT_NH AvPredictedSRT_HI AvMeasuredSRT_NH AvMeasuredSRT_HI]) + max([MeasuredStdErr_NH MeasuredStdErr_HI PredictedStdErr_NH PredictedStdErr_HI]) + 0.5;
    YMin = min([AvPredictedSRT_NH AvPredictedSRT_HI AvMeasuredSRT_NH AvMeasuredSRT_HI]) + min([MeasuredStdErr_NH MeasuredStdErr_HI PredictedStdErr_NH PredictedStdErr_HI]) - 0.5;
    subplot('Position',PanelPos(1,:))
    errorbar(1:4, AvMeasuredSRT_NH(1:4), MeasuredStdErr_NH(1:4),'bs','LineWidth', 2.5), grid on, hold on,
    errorbar(1:4, AvMeasuredSRT_NH(5:8), MeasuredStdErr_NH(5:8),'ro','LineWidth', 2.5)
    set(gca, 'FontSize',12,'XTick', 1:4, 'XTickLabel', {'0'; '10'; '20'; '30'}, 'xlim', [0.8 4.2])
    ylabel('SRT (dB)','FontSize',14,'visible','on'),
    errorbar(1:4, AvPredictedSRT_NH(1:4), PredictedStdErr_NH(1:4),'b--', 'LineWidth', 2.5);
    errorbar(1:4, AvPredictedSRT_NH(5:8), PredictedStdErr_NH(5:8),'r-', 'LineWidth', 2.5);
    ylim([YMin YMax])
    legend('data (co-loc)','data (separ)', 'model (co-loc)', 'model (separ)','Location', 'NorthEast')
    text(1.8, -18, 'NH listeners','FontSize',12)
    subplot('Position',PanelPos(2,:))
    errorbar(1:4, AvMeasuredSRT_HI(1:4), MeasuredStdErr_HI(1:4),'bs','MarkerFaceColor','None','LineWidth', 2.5), grid on, hold on,
    errorbar(1:4, AvMeasuredSRT_HI(5:8), MeasuredStdErr_HI(5:8),'ro','MarkerFaceColor','Auto','LineWidth', 2.5)
    set(gca,'FontSize',12, 'XTick', 1:4, 'XTickLabel', {'0'; '10'; '20'; '30'}, 'TickLabelInterpreter', 'latex', 'xlim', [0.8 4.2], 'YtickLabel','')
    errorbar(1:4, AvPredictedSRT_HI(1:4), PredictedStdErr_HI(1:4),'b--', 'LineWidth', 2.5);
    errorbar(1:4, AvPredictedSRT_HI(5:8), PredictedStdErr_HI(5:8),'r-', 'LineWidth', 2.5);
    ylim([YMin YMax])
    text(2.1, 3.5, 'model: vicente2020')
    temp=['r = ' num2str(round(Corr_Pearson,2))]; text(2.1, 2, temp), 
    temp=['MeanErr = ' num2str(round(MeanError,1)) ' dB']; text(2.1, 1, temp), 
    temp=['MaxErr = ' num2str(round(MaximumError,1)) ' dB']; text(2.1, 0, temp), clear temp
    text(1.8, -18, 'HI listeners','FontSize',12)
    aOverall = axes('Position',[PanelPos(1,1) PanelPos(1,2)-0.025  PanelPos(end,1)+PanelPos(end,3)-PanelPos(1,1) 0.825],'visible','off');
    xlabel('Masker level (dB SL)','FontSize',14,'visible','on')
    
    
    
end

if flags.do_fig6
    % Prediction of exp 1 from Deroche et al (2014) with model prudhomme2020
    
    %Conditions : 
    %Maskers : harmonic masker, F0 50;100;200;400
    %          inharmonic masker, F0 50;100;200;400

    %Program parameters
    Fs = 44100;
    %cut_begin = 0.15*Fs;       %cut silence at beginning of target\masker
    %duration_mask = 2.5*Fs;    %masker duration used for the predictions
    %ref_level = 30;        

    nb_reps = 800;
    r = local_randomset(nb_reps); %create a distribution of 800 for the jitter
    jitter = 0.25.*r;

    %% Load data and target_stats
    %load data from the experiment
    x = amt_load('deroche2014', 'results_DER14a_exp1.mat');
    SRT = x.SRTmean;

    x=amt_load('deroche2014', 'target.mat');
    targetSpec = x.targetSpec;
    target_fc = x.target_fc;
    target = x.target;
    
    %% load maskers
    Mask = {'H50','H100','H200','H400','I50','I100','I200','I400'};
    %conditionName = {'**/H50*.wav','**/H100*.wav','**/H200*.wav','**/H400*.wav','**/I50*.wav','**/I100*.wav','**/I200*.wav','**/I400*.wav'};
    %f0 = [50 100 200 400 50 100 200 400];

    for j = 1:length(Mask)

        [SNR] = amt_cache('get', ['DER14tempPredictions_',Mask{j}], flags.cachemode);
    
        if ~isempty(SNR)
            predicted_SNR(:,j) = mean(SNR,2);
        else
            y=amt_load('deroche2014', 'masker.mat');
            maskers = y.maskers;
            for j = 1:length(Mask)
                SNR = zeros(1,nb_reps);
                masker = amt_load('deroche2014', ['maskers_',Mask{j},'.mat']);
                for k = 1:nb_reps
                    masker_temp = masker.maskers_temp{k};
                    parfor kk=k
                        SNR(:,kk) = prudhomme2020(target,targetSpec,masker_temp,maskers.(Mask{j}).f0,Fs, jitter(kk),target_fc);
                    end
                    amt_cache('set', ['DER14tempPredictions_',Mask{j}], SNR);
                    predicted_SNR(:,j) = mean(SNR,2);
                end
            end
        end
    end
    

    predictSRT = mean(SRT)-(predicted_SNR-mean(predicted_SNR));       %to compare data and prediction, mean of predictions is set to mean of SRT
    Corr_Pearson = corr(SRT', predictSRT','type','Pearson');
    MeanError = mean(abs(SRT-predictSRT));
    LargestError = max(abs(SRT-predictSRT));

    %%
    figure
    plot(1:4,SRT(1:4),'bo','markerfacecolor','b'), grid, hold on
    plot(1:4,SRT(5:8),'ro')
    plot(1:4,predictSRT(1:4),'b-')
    plot(1:4,predictSRT(5:8),'r--')
    xlim([0.7,4.2]), ylim([-22.5,-8.5])
    set(gca,'XTick',[1 2 3 4],'XTicklabel',[50 100 200 400]);
    set(findall(gca, 'Type', 'Line'),'LineWidth',2.5,'MarkerSize',9);
    text(0.8, -16.5, 'model: prudhomme2020')
    temp=['r= ' num2str(round(Corr_Pearson*100)/100)]; text(0.8, -18, temp ), temp=['MeanErr= ' num2str(round(MeanError*10)/10) ' dB']; text(0.8, -19, temp ), temp=['MaxErr= ' num2str(round(LargestError*10)/10) ' dB']; text(0.8, -20, temp ), clear temp
    xlabel('Fundamental frequency (Hz)')
    ylabel('SRT (dB)')
    legend('Harmonic (data)','Inharmonic (data)','Harmonic (model)','Inharmonic (model)')

end

end

function r = local_randomset(nb)
%returns an array of random numbers chosen from a two-parameter 
%probability distribution where the first parameter is a truncated
%normal distribution and the second is given by the input parameter 'nb'

rng('shuffle')
pd = makedist('Normal');
t = truncate(pd,-3,3);
r = random(t,nb,1);
end

function [energy] = local_energy(signal)

energy = sum(abs(signal).^2);

end

function [mean_level, rms_g, rms_d] = local_meanrms(sig)
%compute (in the time domain) the mean of rms level across the left and right channel of a
%stereo file (2-column matrix)

rms_g=20*log10(sqrt(mean(sig(:,1).*sig(:,1))));
rms_d=20*log10(sqrt(mean(sig(:,2).*sig(:,2))));

mean_level=(rms_g+rms_d)/2;
end

function [mean_nrj, nrj_g, nrj_d] = local_mean_nrj(brir)
%compute (in the time domain) the mean energy level across the
%left and right channels of a BRIR

nrj_g=20*log10(sqrt(sum(brir(:,1).*brir(:,1))));
nrj_d=20*log10(sqrt(sum(brir(:,2).*brir(:,2))));

mean_nrj=(nrj_g+nrj_d)/2;
end

function [BinauralRatio, BE_SNR, BU_Advantage, BE_SNR_TFFB, BU_Advantage_TFFB] = ...
    local_predict_RANA18a_Broad_listener_Vicente20(MaskerColocSigs, MaskerSeparSigs, TargetSigs, InternalNoise_L, InternalNoise_R, Ceiling, HannWindowBE, HannWindowBU, WindowOverlap, weightings, FS, fc, Participate)

% HERE PREDICTION FOR 1 LISTENER
% To be used after load_signals_data_BALJ2_HI_Broad_listener_Vicente19 to load data for each listner,

%%% INPUTS
%......

%%% OUTPUTS
% BinauralSNR = model outputs (BE_SNR + BU_Advantage)
% BE_SNR = time-average broadband better-ear SNR
% BU_Advantage = time-average broadband binaural unmasking advantage
% BE_SNR_TFFB = better-ear SNR per time frame and frequency band
% BU_Avantage_TFFB = binaural unmasking advantage per time frame and frequency band

% Conditions:
% target always in front (same in all conditions), not diotic because ear-specific amplification is applied on stimuli at each ear
% maskers (noise-vocoded speech) spatially co-located: both distractors presented to both ears > not diotic because ear-specific amplification is applied on stimuli at each ear
% maskers (noise-vocoded speech) spatially separated: one distrator send to each ear only (no cross-talk, "infinite ILD")
% 4 noise level conditions: 0, 10, 20 and 30 dB SL 

%PREDICTION IN THE DIFFERENT CONDITIONS
% 8 conditions, 1 per column:
% 1-coll 0dBSL = ref condition 
% 2-coll 10 dBSL
% 3-coll 20 dBSL
% 4-coll 30 dBSL
% 5-sep 0 dBSL
% 6-sep 10 dBSL
% 7-sep 20 dBSL
% 8-sep 30 dBSL


BinauralRatio = zeros(1,length(MaskerColocSigs) + length(MaskerSeparSigs));
BE_SNR = zeros(1,length(MaskerColocSigs) + length(MaskerSeparSigs));
BU_Advantage = zeros(1,length(MaskerColocSigs) + length(MaskerSeparSigs));
BE_SNR_TFFB = cell(1,length(MaskerColocSigs) + length(MaskerSeparSigs));
BU_Advantage_TFFB = cell(1,length(MaskerColocSigs) + length(MaskerSeparSigs));
% 
%% Colloc conditions
for i = 1:length(MaskerColocSigs)
    if Participate(i)
        [BinauralRatio(i), BE_SNR(i), BU_Advantage(i), ~, ~] = vicente2020(FS, fc, MaskerColocSigs{i}, TargetSigs{i}, InternalNoise_L(i,:),InternalNoise_R(i,:), Ceiling, HannWindowBE, HannWindowBU, WindowOverlap, weightings);
    else
        BinauralRatio(i) = NaN;
        BE_SNR(i) = NaN;
        BU_Advantage(i) = NaN;
    end
end
%% Separ conditions
for i = 1:length(MaskerSeparSigs)
    if Participate(i)
        [BinauralRatio(i+4), BE_SNR(i+4), BU_Advantage(i+4), ~, ~] = vicente2020(FS, fc, MaskerSeparSigs{i}, TargetSigs{i}, InternalNoise_L(i,:),InternalNoise_R(i,:), Ceiling, HannWindowBE, HannWindowBU, WindowOverlap, weightings);
    else
        BinauralRatio(i+4) = NaN;
        BE_SNR(i+4) = NaN;
        BU_Advantage(i+4) = NaN;
    end
end

%% Model output
BinauralRatio = BE_SNR + BU_Advantage;
end


