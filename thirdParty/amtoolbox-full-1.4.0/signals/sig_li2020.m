function out = sig_li2020(varargin)
% SIG_LI2020 modification of HRTFs/BRIRs for the Li2020 experiments
%
%   Usage: HRTFs = sig_li2020(exp)
%
%   Input parameters:
%     exp   : Experiment (exp A, exp B, exp C, exp D, exp E)
%
%   Output parameters:
%     Obj   : modified HRTFs/BRIRs
%
%   SIG_LI2020 generates modified HRTFs or BRIRs from
%   Li et al. (2020) as specified by the exp flag.
%
%   exp1 - exp4 (size of out matrix): 
%   5x256x5x2 <---> Nr. subjects x HRIR length x conditions x left/right
%   exp 5 (size of out matrix): 
%   5x16384x5x5x2 <---> Nr. subjects x BRIR length x smoothing condition x compression condition x left/right 
%
%   Examples:
%   ---------
%   To get measured modified HRTFs for experiment A:
%
%     sig_li2020('exp1');
%
%   To get measured modified HRTFs for experiment B:
%
%     sig_li2020('exp2');
%
%   To get measured modified HRTFs for experiment C:
%
%     sig_li2020('exp3');
%
%   To get measured modified HRTFs for experiment D:
%
%     sig_li2020('exp4');
%
%   To get measured modified HRTFs for experiment E:
%
%     sig_li2020('exp5');
% 
%
%   Requirements:
%   -------------
%
%   1) Data in hrtf/li2020 and auxdata/li2020 (Raw_BRIR_HRTF_dataset)
%
%
%   References:
%     S. Li, R. Baumgartner, and J. Peissig. Modeling perceived
%     externalization of a static, lateral sound image. Acta Acustica, 4(5),
%     2021.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_li2020.php


%   #Author: Song Li, Institute of Communications Technology, Leibniz
%   University of Hannover, Germany

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.


%% Check input
definput.keyvals.Obj = [];

definput.flags.experiment = {'exp1','exp2','exp3','exp4','exp5'};

[flags,kv]=ltfatarghelper({'Obj'},definput,varargin);

%% HRTF/BRIR manipulation
fs =44100;
hrir_dataset = load('AMT Implementation\data\Raw_BRIR_HRTF_dataset.mat');

if flags.do_exp1 || flags.do_exp2 || flags.do_exp3 || flags.do_exp4
    % Load original HRTFs
    % size: 5 x 256 x 2
    %      5: number of subject  
    %    256: HRIR length    
    %      2:   left & right HRIR 
    dim = size(hrir_dataset.hrir_90_deg);

    hrir_smooth =1; % Opt: smoohting to reduce the irreverent spectral details
    hrir_template = zeros(dim);

    if (hrir_smooth ==1)
        % Opt: smoohting to reduce the irreverent spectral details
        for i = 1: dim(1) 
            hrir_template(i,:,1:2) = sig_hassager2016(squeeze(hrir_dataset.hrir_90_deg(i,:,:)),.2,fs);    
        end    
    else    
        hrir_template = hrir_dataset.hrir_90_deg;
    end
elseif flags.do_exp5


    dim_pre = size(hrir_dataset.brir_90_deg);
    brir_padded=zeros(dim_pre(1),2^nextpow2(length(hrir_dataset.brir_90_deg)),dim_pre(3));
    for i=1:dim_pre(1)
        brir_padded(i,:,:) = padarray(squeeze(hrir_dataset.brir_90_deg(i,:,:)),2^nextpow2(length(hrir_dataset.brir_90_deg))-dim_pre(2),0,'post');
    end
    hrir_dataset.brir_90_deg = brir_padded;
    dim = size(hrir_dataset.brir_90_deg);
    
    brir_smooth =1; % Opt: smoohting to reduce the irreverent spectral details
 
    %window for extraction the direct sound part
    wlen = 44*2; %length of window for raw IR measurements
    wend1=110; %2.5 ms
    win = tukeywin(wlen);

    win_direct=[ones((wend1-wlen/4),1);win(length(win)*3/4+1:end);zeros((dim(2)-wend1),1)];
    win_reverb=1-win_direct;

    %3 direct sound extraction
    
    brir_template_direct_raw = zeros(dim); 
    brir_template_direct = zeros(dim);
    brir_template_reverb = zeros(dim);

    
    
    for i = 1: dim(1)

        brir_template_direct_raw(i,:,1:2)    =    squeeze(hrir_dataset.brir_90_deg(i,:,1:2)) .* win_direct;
        brir_template_reverb(i,:,1:2)        =    squeeze(hrir_dataset.brir_90_deg(i,:,1:2)) .* win_reverb;
    end

    if (brir_smooth ==1)
        for i = 1: dim(1)
            brir_template_direct(i,:,1:2) = sig_hassager2016(squeeze(brir_template_direct_raw(i,:,:)),.2,fs);
        end
    else
        brir_template_direct = brir_template_direct_raw;
    end
  
end


    %% Experiments  
    %% Experiment 1 
    if flags.do_exp1

        %% Exp1: Modification of ILDs
        ILD_increase=[0 5 10 15 20]; % ILD increase level
        len_hrir=256;
        hrir_low_modified_Exp1  = zeros(5,len_hrir,length(ILD_increase),2);    %ILD compression at low frequencies
        hrir_high_modified_Exp1 = zeros(5,len_hrir,length(ILD_increase),2);   %ILD compression at high frequencies
        hrir_bb_modified_Exp1  = zeros(5,len_hrir,length(ILD_increase),2);    %ILD compression at broadband 

            for i = 1: dim(1)

                individual_hrir_template=squeeze(hrir_template(i,:,1:2));

                for k=1:length(ILD_increase)
                    hrir_low_modified_Exp1(i,:,k,:)  = fun_ILD_increase(individual_hrir_template,fs,ILD_increase(k),200,3000);  
                    hrir_high_modified_Exp1(i,:,k,:) = fun_ILD_increase(individual_hrir_template,fs,ILD_increase(k),3000,16000);
                    hrir_bb_modified_Exp1(i,:,k,:)   = fun_ILD_increase(individual_hrir_template,fs,ILD_increase(k),200,16000);
                end

            end

            out.HRTF_low  = hrir_low_modified_Exp1;
            out.HRTF_high = hrir_high_modified_Exp1;    
            out.HRTF_bb   = hrir_bb_modified_Exp1;    

    end
    
    %% Experiment 2 
    if flags.do_exp2 

        %% Exp2: Reaming the ILD, while reducing the magnitude spectral of one ear
        bdwidth=[0 1 4 16 64]; % specra compression level 
        len_hrir=256;
        hrir_modified_Exp2  = zeros(5,len_hrir,length(bdwidth),2);   


        for i = 1: dim(1)

            individual_hrir_template=squeeze(hrir_template(i,:,1:2));

            for k=1:length(bdwidth)
              hrir_modified_Exp2(i,:,k,:) = fun_ipsi_compress(individual_hrir_template,fs,bdwidth(k));
            end
        end

            out.HRTF  = hrir_modified_Exp2;   

    end
    
    
    %% Experiment 3 
    if flags.do_exp3 

        %% Exp3: Reducing the ILD contrast, while keeping the magnitude spectral of one ear
        compression_factor=[1,0.75,0.5,0.25,0];  % ILD compression factor 
        len_hrir=256;
        hrir_ipsi_modified_Exp3    = zeros(5,len_hrir,length(compression_factor),2);    %modified ipsilateral HRTF 
        hrir_contra_modified_Exp3  = zeros(5,len_hrir,length(compression_factor),2);    %modified contralateral HRTF 

        for i = 1: dim(1)
            individual_hrir_template=squeeze(hrir_template(i,:,1:2));

            for k=1:length(compression_factor)

                [hrir_ipsi_modified_Exp3(i,:,k,:), hrir_contra_modified_Exp3(i,:,k,:)]=...
                fun_ILD_compress(individual_hrir_template,fs,compression_factor(k),3000,16000);

            end
        end

            out.HRTF_ipsi  = hrir_ipsi_modified_Exp3;  
            out.HRTF_contra  = hrir_contra_modified_Exp3;   

    end
    
     
    %% Experiment 4 
    if flags.do_exp4 

        %% Exp4: reducing the spectral magnitude of one ear
        bdwidth=[0 1 4 16 64]; 
        len_hrir=256;

        hrir_both_modified_Exp4       = zeros(5,len_hrir,length(bdwidth),2);    %HRTF compression at both ear
        hrir_ipsi_modified_Exp4       = zeros(5,len_hrir,length(bdwidth),2);    %HRTF compression at ipsilateral ear
        hrir_contra_modified_Exp4     = zeros(5,len_hrir,length(bdwidth),2);    %HRTF compression at contralateral ear

        for i = 1: dim(1)
            individual_hrir_template=squeeze(hrir_template(i,:,1:2));   

            for k=1:length(bdwidth)
                     if (k==1)

                             y_modified=sig_hassager2016(individual_hrir_template,0.001,fs); 
                     else         
                             y_modified=sig_hassager2016(individual_hrir_template, bdwidth(k),fs);
                     end   

                     hrir_both_modified_Exp4(i,:,k,:)   = [y_modified(:,1), y_modified(:,2)];
                     hrir_ipsi_modified_Exp4(i,:,k,:)   = [y_modified(:,1), individual_hrir_template(:,2)];
                     hrir_contra_modified_Exp4(i,:,k,:) = [individual_hrir_template(:,1), y_modified(:,2)];

            end
        end

            out.HRTF_ipsi  = hrir_ipsi_modified_Exp4;  
            out.HRTF_contra  = hrir_contra_modified_Exp4;   
            out.HRTF_both  = hrir_both_modified_Exp4; 
    end    
    
    
    %% Experiment 5 
    if flags.do_exp5 

        %% Exp5: reducing the spectral magnitude of one ear
        bdwidth=[0,1,4,16,64];
        compression_factor=[1,0.75,0.5,0.25,0]; 
        len_hrir=256;

        brir_both_modified_Exp5       =    zeros(5,dim(2),length(bdwidth),length(compression_factor), 2);    %BRIR compression at both ear

        for i = 1:dim(1)

                individual_brir_direct_template=squeeze(brir_template_direct(i,:,1:2));
                individual_brir_reverb_template=squeeze(brir_template_reverb(i,:,1:2));    
            for k=1:length(bdwidth)

                for j=1:length(compression_factor)

                    if (k==1)
                      brir_both_modified_Exp5(i,:,k,j,:) = ...
                          [sig_hassager2016(individual_brir_direct_template(1:len_hrir,:),0.01,fs);zeros(dim(2)-len_hrir,2)]+compression_factor(j)*individual_brir_reverb_template;
                    else                     
                      brir_both_modified_Exp5(i,:,k,j,:) = ...
                          [sig_hassager2016(individual_brir_direct_template(1:len_hrir,:),bdwidth(k),fs);zeros(dim(2)-len_hrir,2)]+compression_factor(j)*individual_brir_reverb_template;


                    end

                end

            end

        end

            out.BRIR = brir_both_modified_Exp5;
    end    
    

end


%% functions

function modified_hrir=fun_ILD_increase(hrir_part,fs,ILD_increase,f_low,f_high)

    Nfft=length(hrir_part);
    HRTF_part=fftreal(hrir_part); 
    HRTF_mag_part=db(abs(HRTF_part));

    [~,excessPhaseTF] = RB_minphase(hrir_part,1,'freq');
    
    ILD_ref=HRTF_mag_part(:,1)-HRTF_mag_part(:,2);
    freq = 0:fs/Nfft:fs/2;
    idf = freq >= f_low-1 & freq <= f_high+1; % indices for dedicated frequency range


    modILD_mag=ILD_ref;   % Initiate the modified ILD magnitude
    modILD_mag(idf,:)=ILD_ref(idf,:)+ILD_increase; %compression_factor*meanmag + varmag;

    HRTF_modmag_part=HRTF_mag_part; % Initiate the modified HRTF magnitude

    
    HRTF_modmag_part(idf,1)=HRTF_modmag_part(idf,1);  % left HRTF remained unchanged
    HRTF_modmag_part(idf,2)=HRTF_modmag_part(idf,1)-modILD_mag(idf,:);  % right HRTF spectra change with increased ILD

    
    yRandPhase = ifftreal(10.^(HRTF_modmag_part./20),Nfft);
    YminPhase = RB_minphase(yRandPhase,1,'freq');
    y = ifft(YminPhase.*excessPhaseTF,Nfft);

    modified_hrir=y;

      
end

function modified_hrir=fun_ipsi_compress(hrir_part,fs,compression_factor)


    Nfft=length(hrir_part);
    HRTF_part=fftreal(hrir_part); 
    HRTF_mag_part=db(abs(HRTF_part));

    [~,excessPhaseTF] = RB_minphase(hrir_part,1,'freq');
    
    ILD_ref=HRTF_mag_part(:,1)-HRTF_mag_part(:,2); % reference ILD

    if (compression_factor==0)
        modified_hrir(:,1) = sig_hassager2016(hrir_part(:,1),0.001,fs);  % smoothing process according to Hassager et al. 2016 
    else
        modified_hrir(:,1) = sig_hassager2016(hrir_part(:,1),compression_factor,fs);   
    end
    mag_right = db(abs(fftreal(modified_hrir(:,1))))-ILD_ref;
    yRandPhase = ifftreal(10.^(mag_right./20),Nfft);
    YminPhase = RB_minphase(yRandPhase,1,'freq');
    modified_hrir(:,2) = ifft(YminPhase.*excessPhaseTF(:,2),Nfft);


end

function [modified_ipsi_hrir, modified_contra_hrir]=fun_ILD_compress(hrir_part,fs,compression_factor,f_low,f_high)


    Nfft=length(hrir_part);
    HRTF_part=fftreal(hrir_part,Nfft); 
    HRTF_mag_part=db(abs(HRTF_part));

    [~,excessPhaseTF] = RB_minphase(hrir_part,1,'freq');
    
    ILD_ref=HRTF_mag_part(:,1)-HRTF_mag_part(:,2); % reference ILD
    freq = 0:fs/Nfft:fs/2;
    idf = freq >= f_low-1 & freq <= f_high+1; % indices for dedicated frequency range
    idwf = idf(:) | circshift(idf(:),[1,0]); % include one neighbouring position for evaluation of frequency weighting
    wf = diff(freqtoerb(freq(idwf))); % frequency weighting according to differentiated ERB scale
    wf = wf(:)/sum(wf);
    
    
    mag = ILD_ref(idf,:); % HRTF magnitudes in dB       
    meanmag = sum(wf.*mag,1);
    varmag = mag - meanmag;

    modmag = ILD_ref;
    modmag(idf,:)= meanmag + compression_factor*varmag;

    HRTF_modmag_part_contra=HRTF_mag_part;
    HRTF_modmag_part_ipsi=HRTF_mag_part;
    
 
    HRTF_modmag_part_contra(idf,1)=HRTF_modmag_part_contra(idf,1);
    HRTF_modmag_part_contra(idf,2)=HRTF_modmag_part_contra(idf,1)-modmag(idf,:);

    HRTF_modmag_part_ipsi(idf,2)=HRTF_modmag_part_ipsi(idf,2);
    HRTF_modmag_part_ipsi(idf,1)=HRTF_modmag_part_ipsi(idf,2)+modmag(idf,:);


    yRandPhase = ifftreal(10.^(HRTF_modmag_part_contra./20),Nfft);
    YminPhase = RB_minphase(yRandPhase,1,'freq');
    modified_contra_hrir = ifft(YminPhase.*excessPhaseTF,Nfft);

    yRandPhase = ifftreal(10.^(HRTF_modmag_part_ipsi./20),Nfft);
    YminPhase = RB_minphase(yRandPhase,1,'freq');
    modified_ipsi_hrir = ifft(YminPhase.*excessPhaseTF,Nfft);
            

end

function [minPhase,excessPhase] = RB_minphase(IR,dim,TFdomainFlag)
% RB_minphase - create minimum-phase filter via causal cepstrum
%
% Usage: [minPhase,excessPhase] = RB_minphase(IR,dim,TFdomainFlag)
% amt toolbox
% RB, 2016/6/3

Nfft = 2.^nextpow2(size(IR,dim));

TF = fft(IR,Nfft,dim);
logTF = log(abs(TF)+eps);

cep = ifft(logTF,Nfft,dim);
Nshift = mod(dim-1,ndims(cep));
cep1 = shiftdim(cep,Nshift);
cep1(Nfft/2+2:Nfft,:,:,:,:) = 0;    % set non-causal part to zero and 
cep1(2:Nfft/2,:,:,:,:) = 2*cep1(2:Nfft/2,:,:,:,:);    % multiply causal part by 2 (due to symmetry)
cepMinPhase = shiftdim(cep1,ndims(cep)-Nshift);

logTFminPhase = fft(cepMinPhase,Nfft,dim);
TFminPhase = exp(logTFminPhase);

switch TFdomainFlag
  case 'freq'
    minPhase = TFminPhase;
  case 'time'
    minPhase = ifft(TFminPhase,Nfft,dim);
end


if nargout == 2
  switch TFdomainFlag
    case 'freq'
      excessPhase = TF./TFminPhase;
    case 'time'
      excessPhase = ifft(TF./TFminPhase,Nfft,dim);
  end
end

end



