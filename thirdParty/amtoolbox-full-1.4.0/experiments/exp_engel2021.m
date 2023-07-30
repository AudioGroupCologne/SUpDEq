function exp_engel2021(varargin)
%EXP_ENGEL2021 Simulations of Engel et al. (2021)
%
%   Usage: data = exp_engel2021(flag)
%
%
%   Input parameters:
%     fig1   : spherical harmonics spectra
%
%     fig2   : mag/phase errors per spatial order
%
%     fig3   : ITD/ILD per spatial order
%
%     fig4   : mag/phase errors per method
%
%     fig5   : ITD/ILD per method
%
%     fig6   : loudness per direction
%
%     fig7   : models' outputs per spatial order
%
%   EXP_ENGEL2021(flag) reproduces figures of the study from 
%   Engel et al. (2021).
%
%   This script will look for cached data of hnm and results obtained from
%   running the 3 models baumgartner2021, jelfs2011, and reijniers2014. 
%
%   If these are not found, the data will be generated from a set of HRTF.
%
%   WARNINGS: The HRTF data (hnm) is fairly large (~9GB) and may take several 
%   minutes to generate. Most of it is not needed for the plots.
%
%   The results data is smaller (~300MB) but may take LONG to generate 
%   (it took ~12 hours on a consumer-grade laptop. It is recommended to
%   download it if available online.
%
%
%   Examples:
%   ---------
%
%   To display Fig.1 of Engel(2021) use :
%
%     exp_engel2021('fig1');
%
%   To display Fig.2 of Engel(2021) use :
%
%     exp_engel2021('fig2');
%
%   To display Fig.3 of Engel(2021) use :
%
%     exp_engel2021('fig3');
%
%   To display Fig.4 of Engel(2021) use :
%
%     exp_engel2021('fig4');
%
%   To display Fig.5 of Engel(2021) use :
%
%     exp_engel2021('fig5');
%
%   To display Fig.6 of Engel(2021) use :
%
%     exp_engel2021('fig6');
%
%   To display Fig.7 of Engel(2021) use :
%
%     exp_engel2021('fig7');
%
%   See also: baumgartner2021 jelfs2011 reijniers2014
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_engel2021.php


%   #Author: Isaac Engel (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.flags.experiment = {'missingflag',...
  'fig1','fig2','fig3','fig4','fig5','fig6','fig7'};
definput.import={'amt_cache'};
[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.experiment{2:end-2}),...
             sprintf('%s or %s',definput.flags.experiment{end-1},definput.flags.experiment{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

%% Parameters
N_vec = 1:44; % spatial orders to be tested
Nref = 44; % reference (very high) spatial order
r = 0.0875; % nominal head radius = 8.75 cm
itdJND = 20; % ITD JND for plots, according to Klockgether 2016
ildJND = 0.6; % ILD JND for plots, according to Klockgether 2016
pagewidth = 17.9112; % in cm
colwidth = 8.6472; % in cm

% HRTF file name (without the .sofa)
hrirname = 'FABIAN_HRIR_measured_HATO_0';

%% Load HRTF

SOFA_obj = amt_load('engel2021',[hrirname,'.sofa']);
[h,fs,az,el] = sofa2hrtf(SOFA_obj); % get HRTF data in a convenient format

% Zero-pad the HRIRs to increase frequency resolution
taps = 2048;
h = [h;zeros(taps-size(h,1),size(h,2),size(h,3))];

H = ffth(h); % to frequency domain
nfreqs = size(H,1); % number of frequency bins
f = linspace(0,fs/2,nfreqs).'; % vector with frequencies

%% Defining test conditions

ncond = 9; 
test_conditions = cell(ncond,1);
% 1. Trunc (no preprocessing, no tapering, no EQ)
test_conditions{1}.name = 'Trunc';
test_conditions{1}.preproc = 'Trunc';
test_conditions{1}.tap = 0;
test_conditions{1}.dualBand = false;
test_conditions{1}.EQ = 0;
% 2. EQ (no preprocessing, no tapering, HRF EQ)
test_conditions{2}.name = 'EQ';
test_conditions{2}.preproc = 'Trunc';
test_conditions{2}.tap = 0;
test_conditions{2}.dualBand = false;
test_conditions{2}.EQ = 2;
% 3. Tap (no preprocessing, dual band tapering with Hann weights, HRF EQ)
test_conditions{3}.name = 'Tap';
test_conditions{3}.preproc = 'Trunc';
test_conditions{3}.tap = 2;
test_conditions{3}.dualBand = true;
test_conditions{3}.EQ = 2;
% 4. TA (time-aligned HRTF, no tapering, no EQ)
test_conditions{4}.name = 'TA';
test_conditions{4}.preproc = 'TA';
test_conditions{4}.tap = 0;
test_conditions{4}.dualBand = false;
test_conditions{4}.EQ = 0;
% 5. MagLS (MagLS, no tapering, no EQ)
test_conditions{5}.name = 'MagLS';
test_conditions{5}.preproc = 'MagLS';
test_conditions{5}.tap = 0;
test_conditions{5}.dualBand = false;
test_conditions{5}.EQ = 0;
% 6. MagLS+CC (MagLS, no tapering, CovConst)
test_conditions{6}.name = 'MagLSCC';
test_conditions{6}.preproc = 'MagLS';
test_conditions{6}.tap = 0;
test_conditions{6}.dualBand = false;
test_conditions{6}.EQ = 1;
% 7. SpSub (SpSub, no tapering, no EQ)
test_conditions{7}.name = 'SpSub';
test_conditions{7}.preproc = 'SpSub';
test_conditions{7}.tap = 0;
test_conditions{7}.dualBand = false;
test_conditions{7}.EQ = 0;
% 8. SpSubMod (SpSub, dual band tapering with Hann weights, HRF EQ)
test_conditions{8}.name = 'SpSubMod';
test_conditions{8}.preproc = 'SpSubMod';
test_conditions{8}.tap = 2;
test_conditions{8}.dualBand = true;
test_conditions{8}.EQ = 2;
% 9. BiMagLS
test_conditions{9}.name = 'BiMagLS';
test_conditions{9}.preproc = 'BiMagLS';
test_conditions{9}.tap = 0;
test_conditions{9}.dualBand = false;
test_conditions{9}.EQ = 0;

%% Define a few direction subsets

% 1. Nearest neighbours to 110-point Lebedev grid
gridTestReq = sofia_lebedev(110,0);
[~,indLeb] = getGridSubset([az,el],gridTestReq,0);
indLeb = find(indLeb); % get indices
azLeb = az(indLeb);
elLeb = el(indLeb);
wLeb = gridTestReq(:,3);
% 2. Median plane (lat = 0 +/- 1 deg)
[lat,~] = sph2hor(az*180/pi,90-el*180/pi); % lat/pol coordinates
indMP = find(abs(lat)<1); 
azMP = az(indMP);
elMP = el(indMP);
% 3. Horizontal plane (el = 90 +/- 1 deg)
indHP = find(abs(el-pi/2)<(pi/180));
azHP = az(indHP);
% elHP = el(indHP);
indFront = find(abs(el-pi/2)<(pi/180) & abs(az)<(pi/180)); % front

%% Get high-order HRTF

amt_disp(' ');
amt_disp('Getting reference HRTF...');

hnm_filename = 'hnm_ref';
[Hnm_mag_ref, hnm_ref, Hnm_TA_ref] = amt_cache('get',hnm_filename,flags.cachemode);

if isempty(hnm_ref)

    % Standard
    hnm_ref = toSH(h,Nref,'az',az,'el',el,'fs',fs,'mode','Trunc');
    % Time-aligned
    hnm_TA_ref = toSH(h,Nref,'az',az,'el',el,'fs',fs,'mode','TA');
    Hnm_TA_ref = ffth(hnm_TA_ref);
    % Mag-only
    Y = AKsh(Nref,[],az*180/pi,el*180/pi,'real').';
    Hnm_mag_ref = mult3(abs(H),pinv(Y));
    amt_cache('set',hnm_filename,Hnm_mag_ref,hnm_ref,Hnm_TA_ref,fs);
    
end
Hnm_ref = ffth(hnm_ref); % to frequency domain

amt_disp('Got reference HRTF!');

%% Get results for high-order HRTF

amt_disp(' ');
amt_disp('Getting results for reference HRTF...');

results_filename = 'res_ref';
results = amt_cache('get', results_filename, flags.cachemode);

if isempty(results)
    
    % --- Numerical analysis --- %
    % Magnitude and phase delay 110-point Lebedev grid
    results.mag = 20*log10(abs(H(:,indLeb,:)));
    results.pd = div2( -unwrap(angle(H(:,indLeb,1)))+unwrap(angle(H(:,indLeb,2))) , 2*pi*f )*1e6;
    % Loudness across directions
    [L,wERB,calibrationGain] = perceptualSpectrum(H,fs);
    Lavgd = sum( mult2(L,wERB) , 1 );
    results.L = Lavgd;
    results.L_perDirection = L; % save to calculate PSD
    results.Lavg = sum( mult2(Lavgd(:,indLeb,:),wLeb.') , 2 ); % weighted avg on the Lebedev grid
    results.calibrationGain = calibrationGain; % save to calibrate loudness
    % Interaural differences for horizontal plane
    results.itd = itdestimator(permute(h(:,indHP,:),[2,3,1]),'fs',fs,...
        'MaxIACCe','lp','upper_cutfreq', 3000,...
        'butterpoly', 10)*1e6;
    results.ild = getILD(h(:,indHP,:),fs);

    % --- Reijniers 2014 --- %
    amt_disp('Running reijniers2014...');

    % Make DTF and SOFA object for Lebedev grid directions only
    dtf = getDTF(h(:,indLeb,:),fs);
    SOFA_obj = hrtf2sofa(dtf,fs,azLeb,elLeb);
    % Preprocessing source information (demo_reijniers2014)
    [template_loc, target] = reijniers2014_featureextraction(SOFA_obj);
    % Run virtual experiments (demo_reijniers2014)
    num_exp = 100;
    [doa, ~] = reijniers2014(template_loc, target, 'num_exp', num_exp);       
    % Calculate performance measures (demo_reijniers2014)
    results.lat_acc = reijniers2014_metrics(doa, 'accL'); % mean lateral error
    results.lat_prec = reijniers2014_metrics(doa, 'precL'); % lateral std
    results.pol_acc = reijniers2014_metrics(doa, 'accP'); % mean polar error
    results.pol_prec = reijniers2014_metrics(doa, 'precP'); % polar std
    results.template_loc = template_loc; % template DTF
    
    % --- Baumgartner 2021 --- %
    amt_disp('Running baumgartner2021...');

    % Make DTF for median plane directions only
    dtf = getDTF(h(:,indMP,:),fs);
    ndirs = numel(elMP);
    results.ext = nan(ndirs,1);
    template_ext = cell(ndirs,1);
    for j=1:ndirs
        template_ext{j} = hrtf2sofa(dtf(:,j,:),fs,azMP(j),elMP(j));
        % Get externalisation values
        results.ext(j) = baumgartner2021(template_ext{j},template_ext{j});
    end
    results.template_ext = template_ext; % template DTFs
    
    % --- Jelfs 2011 --- %
    amt_disp('Running jelfs2011...');

    ndirs = numel(azHP);
    srm = nan(ndirs,1);
    target = squeeze(h(:,indFront,:)); % target fixed at front
    for j = 1:ndirs
        interferer = squeeze(h(:,indHP(j),:)); % interferer moves around the HP
        srm(j) = jelfs2011(target,interferer,fs);
    end 
    results.srm = srm;

    % --- Save results --- %
    amt_cache('set', results_filename, results);
    
end

% Some parameters needed for the calculations below:
mag_ref = results.mag; % to calculate magnitude error
pd_ref = results.pd; % to calculate phase delay error
Lref = results.L_perDirection; % to calculate PSD
Lavg_ref = results.Lavg; % to level-match loudness
calibrationGain = results.calibrationGain; % to calibrate perceptual spectrum
template_loc = results.template_loc; % template localisation model
template_ext = results.template_ext; % template externalisation model
clear results

amt_disp('Got results for reference HRTF!');

%% Get results for all low-order HRTFs

amt_disp(' ');
amt_disp('Getting results for all low-order HRTFs...');

nhrtfs = numel(N_vec)*ncond;
hrtfcount = 0;

for N=N_vec % iterate through spatial orders

    for i=1:ncond % iterate through conditions
        
        hrtfcount = hrtfcount+1;
                
        name = test_conditions{i}.name;
        results_filename = sprintf('res_ord%0.2d_%s',N,name);
        results = amt_cache('get',results_filename,flags.cachemode);
        
        if isempty(results) % generate results if not available
            
            amt_disp(' ');
            amt_disp(sprintf('Processing HRTF %d/%d...',hrtfcount,nhrtfs));
            
            hnm_filename = sprintf('hnm_ord%0.2d_%s',N,name);
            hnm = amt_cache('get',hnm_filename,flags.cachemode);
            
            if isempty(hnm) % generate hnm if not available
                
                amt_disp(['Generating ',num2str(hnm_filename),'...']);
                [hnm,fs,varOut] = toSH(h,N,'az',az,'el',el,'fs',fs,...
                    'mode' ,test_conditions{i}.preproc,...
                    'tapering',test_conditions{i}.tap,...
                    'dualBand',test_conditions{i}.dualBand,...
                    'EQ',test_conditions{i}.EQ,...
                    'Hnm_ref',Hnm_ref); % reference Hnm only used for EQ
                amt_cache('set', hnm_filename, hnm, fs, varOut);
                
            end
            
            % --- Interpolate HRTF --- %
            isaligned = strcmp(name,'TA') || strcmp(name,'BiMagLS');
            hInterp = fromSH(hnm,fs,az,el,isaligned,r);
            HInterp = ffth(hInterp);

            % --- Numerical analysis --- %
            % Magnitude and phase delay error 110-point Lebedev grid
            mag = 20*log10(abs(HInterp(:,indLeb,:)));
            pd = div2( -unwrap(angle(HInterp(:,indLeb,1)))+unwrap(angle(HInterp(:,indLeb,2))) , 2*pi*f )*1e6;
            results.err_mag = mean(abs(mag-mag_ref),2); % avg abs difference across directions
            results.err_pd = mean(abs(pd-pd_ref),2); % avg abs difference across directions
            % Loudness across directions
            [L,wERB] = perceptualSpectrum(HInterp,fs,calibrationGain);
            Lavgd = sum( mult2(L,wERB) , 1 ); % ERB-weighted avg loudness per direction
            Lavg = sum( mult2(Lavgd(:,indLeb,:),wLeb.') , 2 ); % weighted avg on the Lebedev grid
            % L = L - Lavg + Lavg_ref; % level-match with reference
            L = L - repmat(Lavg - Lavg_ref,size(L,1),size(L,2),1); % level-match with reference
            Lavgd = sum( mult2(L,wERB) , 1 ); % recalculate avg loudness per direction
            PSD = sum( mult2(abs(L-Lref),wERB) , 1 ); % avg abs difference per direction (PSD)
            PSDavg = sum( mult2(PSD(:,indLeb,:),wLeb.') , 2 ); % avg PSD over the Lebedev grid
            results.L = Lavgd;
            results.PSD = PSD;
            results.PSDavg = PSDavg;
            % Interaural differences for horizontal plane
            results.itd = itdestimator(permute(hInterp(:,indHP,:),[2,3,1]),'fs',fs,...
                'MaxIACCe','lp','upper_cutfreq', 3000,...
                'butterpoly', 10)*1e6;
            results.ild = getILD(hInterp(:,indHP,:),fs);

            % ---  Reijniers 2014 --- %
            amt_disp('Running reijniers2014...');

            % Make DTF and SOFA object for Lebedev grid directions only
            dtf = getDTF(hInterp(:,indLeb,:),fs);
            SOFA_obj = hrtf2sofa(dtf,fs,azLeb,elLeb);
            % Preprocessing source information (demo_reijniers2014)
            [~, target] = reijniers2014_featureextraction(SOFA_obj);
            % Run virtual experiments (demo_reijniers2014)
            num_exp = 100;
            [doa, ~] = reijniers2014(template_loc, target, 'num_exp', num_exp);       
            % Calculate performance measures (demo_reijniers2014)
            results.lat_acc = reijniers2014_metrics(doa, 'accL'); % mean lateral error
            results.lat_prec = reijniers2014_metrics(doa, 'precL'); % lateral std
            results.pol_acc = reijniers2014_metrics(doa, 'accP'); % mean polar error
            results.pol_prec = reijniers2014_metrics(doa, 'precP'); % polar std
            
            % --- Baumgartner 2021 --- %
            amt_disp('Running baumgartner2021...'); 

            % Make DTF for median plane directions only
            dtf = getDTF(hInterp(:,indMP,:),fs);
            ndirs = numel(elMP);
            results.ext = nan(ndirs,1);
            for j=1:ndirs
                target = hrtf2sofa(dtf(:,j,:),fs,azMP(j),elMP(j));
                % Get externalisation values
                results.ext(j) = baumgartner2021(target,template_ext{j});
            end
            
            % --- Jelfs 2011 --- %
            amt_disp('Running jelfs2011...'); 

            ndirs = numel(azHP);
            srm = nan(ndirs,1);
            target = squeeze(hInterp(:,indFront,:)); % target fixed at front
            for j = 1:ndirs
                interferer = squeeze(hInterp(:,indHP(j),:)); % interferer moves around the HP
                srm(j) = jelfs2011(target,interferer,fs);
            end 
            results.srm = srm;

            % --- Save results --- %
            amt_cache('set',results_filename,results);
            
        end
    end
end

amt_disp('Got results for all low-order HRTFs!');

%% Plot Fig. 1 (SH spectra)

if flags.do_fig1    
    
    fig1size = [pagewidth 4.5];
    fig1 = figure('units','centimeters','PaperUnits','centimeters',...
        'PaperSize',fig1size,'Renderer','painters',...
        'pos',[2 2 fig1size(1) fig1size(2)],...
        'paperposition',[0 0 fig1size(1) fig1size(2)]);
    % Tight subplot
    gap = [.02 .008]; % gap between subplots in norm units (height width)
    marg_h = [.18 .03]; % figure height margins in norm units (lower upper)
    marg_w = [.05 .08]; % figure width margins in norm units (left right)
    [ha, ~] = tight_subplot(1,3,gap,marg_h,marg_w);
    % Fig. 1a: original HRTF
    axes(ha(1)), plotSHenergy(Hnm_ref(:,:,1),fs); % left only
    xlabel('f (Hz)'), ylabel('Order (n)')
    axpos = get(ha(1),'pos');
    annotation(fig1,'textbox',...
        [axpos(1)+axpos(3)-0.2 axpos(2)+axpos(4)-0.09 0.2 0.09],...
        'String',{'a) Original'},...
        'HorizontalAlignment','right',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    set(gca,'fontsize',7)
    % Fig. 1b: time-aligned HRTF
    axes(ha(2)), plotSHenergy(Hnm_TA_ref(:,:,1),fs);
    xlabel('f (Hz)'), ylabel('');
    set(ha(2),'YTickLabel',{})
    axpos = get(ha(2),'pos');
    annotation(fig1,'textbox',...
        [axpos(1)+axpos(3)-0.2 axpos(2)+axpos(4)-0.09 0.2 0.09],...
        'String',{'b) Time-aligned'},...
        'HorizontalAlignment','right',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    set(gca,'fontsize',7)
    % Fig. 1c: magnitude-only HRTF
    axes(ha(3)), plotSHenergy(Hnm_mag_ref(:,:,1),fs);
    xlabel('f (Hz)'), ylabel('');
    set(ha(3),'YTickLabel',{})
    axpos = get(ha(3),'pos');
    annotation(fig1,'textbox',...
        [axpos(1)+axpos(3)-0.2 axpos(2)+axpos(4)-0.09 0.2 0.09],...
        'String',{'c) Magnitude only'},...
        'HorizontalAlignment','right',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    set(gca,'fontsize',7)
    c = colorbar; c.Label.String = 'Energy (dB)';
    c.Position = [0.9282    0.1784    0.0199    0.7934];

end

%% Plot Fig. 2 (mag/phase errors per spatial order)

if flags.do_fig2     
    % Load data from file
    names = {
        'ord01_Trunc'
        'ord05_Trunc'
        'ord10_Trunc'
        'ord20_Trunc'
        'ord30_Trunc'
        'ord44_Trunc'
    };
    labels = {
        'Trunc (N=1)'
        'Trunc (N=5)'
        'Trunc (N=10)'
        'Trunc (N=20)'
        'Trunc (N=30)'
        'Trunc (N=44)'  
    };
    n = numel(names);
    err_mag = zeros(nfreqs,n);
    err_pd = zeros(nfreqs,n);
    for i=1:n
        name = ['res_',names{i}];
        results = amt_cache('get', name, flags.cachemode); 
        err_mag(:,i) = results.err_mag(:,1,1); % left ear only
        err_pd(:,i) = results.err_pd;
    end

    % Tight subplot
    fig2size = [colwidth, 8.47];
    fig2 = figure('units','centimeters','PaperUnits','centimeters',...
        'PaperSize',fig2size,'Renderer','painters',...
        'pos',[2 2 fig2size(1) fig2size(2)],...
        'paperposition',[0 0 fig2size(1) fig2size(2)]);
    gap = [.02 .008]; % gap between subplots in norm units (height width)
    marg_h = [.09 .02]; % [.09 .1] figure height margins in norm units (lower upper)
    marg_w = [.1 .03]; % figure width margins in norm units (left right)
    [ha, ~] = tight_subplot(2,1,gap,marg_h,marg_w);

    % Top plot: magnitude error
    colors = parula(n+1);
    lsvec = {'-'};%,'-.'};
    lwvec = [0.5,0.5];
    mvec = {'^','v','x','s','o','d','p','h','>'};
    ms = 3; % marker size
    mi = int32(logspace(log10(1),log10(1025),10));

    axes(ha(1))
    for i=1:n
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        % For Matlab version >= R2016b:
        %semilogx(f,err_mag(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
        % For older versions:
        semilogx(f,err_mag(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls); hold on
        semilogx(f(mi),err_mag(mi,i),'Color',colors(i,:),'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end

    semilogx([f(2) 20000], [1 1],'k:')
    legend(labels,'location','west')

    xlabel(''), xlim([f(2) 20000])
    set(ha(1),'XTickLabel',{})
    ylabel('Error (dB)'), grid on
    axpos = get(ha(1),'pos');
    annotation(fig2,'textbox',...
        [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
        'String',{'a. Magnitude error'},...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FitBoxToText','off',...
        'EdgeColor','none');

    % Bottom plot: phase error
    axes(ha(2))
    for i=1:n
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        % For Matlab version >= R2016b:
        %semilogx(f,err_pd(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
        % For older versions:
        semilogx(f,err_pd(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls); hold on
        semilogx(f(mi),err_pd(mi,i),'Color',colors(i,:),'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end

    semilogx([f(2) 20000], [20 20],'k:')
    xlim([f(2) 20000]), grid on, ylabel('Error (\mus)')
    set(ha(2),'XTick',[100,1000,10000,20000],'XTickLabel',{'100','1k','10k','20k'})
    xlabel('f (Hz)')
    axpos = get(ha(2),'pos');
    annotation(fig2,'textbox',...
        [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
        'String',{'b. Phase delay error'},...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FitBoxToText','off',...
        'EdgeColor','none');

    set(ha,'fontsize',7)

end

%% Plot Fig. 3 (ITD/ILD per spatial order)

if flags.do_fig3    
    % Load data from file
    names = {
        'ord01_Trunc'
        'ord05_Trunc'
        'ord10_Trunc'
        'ord20_Trunc'
        'ord30_Trunc'
        'ord44_Trunc'
        'ref'
    };
    labels = {
        'N=1'
        'N=5'
        'N=10'
        'N=20'
        'N=30'
        'N=44'
        'Reference'
    };
    n = numel(names);
    ndirs = numel(azHP);
    itd = zeros(ndirs,n);
    ild = zeros(ndirs,n);
    for i=1:n
        name = ['res_',names{i}];
        results = amt_cache('get', name, flags.cachemode); 
        itd(:,i) = results.itd;
        ild(:,i) = results.ild;
    end

    % Tight subplot  
    fig3size = [pagewidth 7.3025];
    fig3 = figure('units','centimeters','pos',[2 2 fig3size(1) fig3size(2)],...
        'Renderer','painters','PaperSize',[fig3size(1) fig3size(2)],...
        'paperposition',[0 0 fig3size(1) fig3size(2)]);
    gap = [.04 .06]; % gap between subplots in norm units (height width)
    marg_h = [0 .05]; % [.09 .1] figure height margins in norm units (lower upper)
    marg_w = [.03 .23]; % figure width margins in norm units (left right)

    [ha, ~] = tight_subplot(1,2,gap,marg_h,marg_w); 

    colors = parula(n);
    lsvec = {'-'};%,'-.'};
    lwvec = [0.5,0.5];
    mvec = {'^','v','x','s','o','d','p','h','>'};
    ms = 3; % marker size
    step_big = round(numel(azHP)/4);
    step_small = round(step_big/n);

    % ITD
    axes(ha(1))
    for i=1:n % multiply by 0.001 to have it in ms (cleaner plot)
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
        color = colors(i,:);
        if i==n % reference
            ls = ':';
            lw = 0.5;
            m = 'none';
            color = [0 0 0];
        end
        % For Matlab version >= R2016b:
        %polarplot(azHP,abs(itd(:,i))*0.001,'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on
        % For older versions:
        polarplot(azHP,abs(itd(:,i))*0.001,'Color',color,'LineWidth',lw,'LineStyle',ls); hold on
        polarplot(azHP(mi),abs(itd(mi,i))*0.001,'Color',color,'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end
    set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
    set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',{'0','45','90','135','','-135','-90','-45'},'RAxisLocation',-90)
    title('ITD (ms)')

    % ILD
    axes(ha(2))
    for i=1:n
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
        color = colors(i,:);
        if i==n % reference
            ls = ':';
            lw = 0.5;
            m = 'none';
            color = [0 0 0];
        end
        % For Matlab version >= R2016b:
        % polarplot(azHP,abs(ild(:,i)),'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on
        % For older versions:
        polarplot(azHP,abs(ild(:,i)),'Color',color,'LineWidth',lw,'LineStyle',ls); hold on
        polarplot(azHP(mi),abs(ild(mi,i)),'Color',color,'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end
    set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
    set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',{'0','45','90','135','','-135','-90','-45'},'RAxisLocation',-90)
    title('ILD (dB)')
    legend(labels,'position',[0.8217    0.3176    0.1331    0.3094]);

    % Violins
    fig3bsize = [fig3size(1) 3.3];
    fig3b = figure('units','centimeters','pos',[2 2 fig3bsize(1) fig3bsize(2)],...
        'Renderer','painters','PaperSize',[fig3bsize(1) fig3bsize(2)],...
        'paperposition',[0 0 fig3bsize(1) fig3bsize(2)]);
    gap = [.04 .06]; % gap between subplots in norm units (height width)
    marg_h = [0.1 0.03];
    marg_w = [0.06 0.01];
    [ha, ~] = tight_subplot(1,2,gap,marg_h,marg_w); 

    itderr = abs(itd-repmat(itd(:,end), 1, size(itd, 2)));
    ilderr = abs(ild-repmat(ild(:,end), 1, size(ild, 2)));

    axes(ha(1)), violinplot(itderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
    set(gca,'YTickLabelMode','auto')
    hold on, plot([0 n],[itdJND itdJND],'k:')
    set(gca,'XTickLabel',labels(1:end-1))
    ylabel('Abs. ITD error (ms)')
    set(gca,'fontsize',7);

    axes(ha(2)), violinplot(ilderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
    set(gca,'YTickLabelMode','auto')
    hold on, plot([0 n],[ildJND ildJND],'k:')
    set(gca,'XTickLabel',labels(1:end-1))
    ylabel('Abs. ILD error (dB)')
    set(gca,'fontsize',7);

end

%% Plot Fig. 4 (mag/phase errors per method)

if flags.do_fig4    
    % Load data from file
    names = {
        'ord03_Trunc'
        'ord03_EQ'
        'ord03_Tap'
        'ord03_TA'
        'ord03_MagLS'
        'ord03_MagLSCC'
        'ord03_SpSub'
        'ord03_SpSubMod'
        'ord03_BiMagLS'
    };
    labels = {
        'Trunc'
        'EQ'
        'Tap'
        'TA'
        'MagLS'
        'MagLS+CC'
        'SpSub'
        'SpSubMod'
        'BiMagLS'
    };
    n = numel(names);
    err_mag = zeros(nfreqs,n);
    err_pd = zeros(nfreqs,n);
    for i=1:n
        name = ['res_',names{i}];
        results = amt_cache('get', name, flags.cachemode); 
        err_mag(:,i) = results.err_mag(:,1,1); % left ear only
        err_pd(:,i) = results.err_pd;
    end
    N = 3;
    c = 343;
    fa = N*c/(2*pi*r); % aliasing frequency

    % Tight subplot
    fig4size = [colwidth, 8.47];
    fig4 = figure('units','centimeters','PaperUnits','centimeters',...
        'PaperSize',fig4size,'Renderer','painters',...
        'pos',[2 2 fig4size(1) fig4size(2)],...
        'paperposition',[0 0 fig4size(1) fig4size(2)]);
    gap = [.02 .008]; % gap between subplots in norm units (height width)
    marg_h = [.09 .02]; % [.09 .1] figure height margins in norm units (lower upper)
    marg_w = [.1 .03]; % figure width margins in norm units (left right)
    [ha, ~] = tight_subplot(2,1,gap,marg_h,marg_w);

    % Top plot: magnitude error
    colors = parula(n+1);
    lsvec = {'-'};%,'-.'};
    lwvec = [0.5,0.5];
    mvec = {'^','v','x','s','o','d','p','h','>'};
    ms = 3; % marker size
    mi = int32(logspace(log10(1),log10(1025),10));

    axes(ha(1))
    for i=1:n
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        % For Matlab version >= R2016b:
        %semilogx(f,err_mag(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
        % For older versions:
        semilogx(f,err_mag(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls); hold on
        semilogx(f(mi),err_mag(mi,i),'Color',colors(i,:),'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end
    xlabel(''), xlim([f(2) 20000])
    set(ha(1),'XTickLabel',{})
    semilogx([fa fa],get(gca,'ylim'),'k--')
    semilogx(get(gca,'xlim'), [1 1],'k:')
    legend(labels,'position',[0.1167    0.5817    0.2827    0.3392])

    ylabel('Error (dB)'), grid on
    axpos = get(ha(1),'pos');
    annotation(fig4,'textbox',...
        [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
        'String',{'a. Magnitude error'},...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FitBoxToText','off',...
        'EdgeColor','none');

    % Bottom plot: phase error
    axes(ha(2))
    for i=1:n
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        % For Matlab version >= R2016b:
        %semilogx(f,err_pd(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
        % For older versions:
        semilogx(f,err_pd(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls); hold on
        semilogx(f(mi),err_pd(mi,i),'Color',colors(i,:),'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end
    xlim([f(2) 20000]), grid on, ylabel('Error (\mus)')
    semilogx([fa fa],get(gca,'ylim'),'k--')
    semilogx(get(gca,'xlim'), [20 20],'k:')
    set(ha(2),'XTick',[100,1000,10000,20000],'XTickLabel',{'100','1k','10k','20k'})
    xlabel('f (Hz)')
    axpos = get(ha(2),'pos');
    annotation(fig4,'textbox',...
        [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
        'String',{'b. Phase delay error'},...
        'HorizontalAlignment','left',...
        'FontWeight','bold',...
        'FontSize',7,...
        'FitBoxToText','off',...
        'EdgeColor','none');

    set(ha,'fontsize',7)

end

%% Plot Fig. 5 (ITD/ILD per method)

if flags.do_fig5    
    % Load data from file
    names = {
        'ord03_Trunc'
        'ord03_EQ'
        'ord03_Tap'
        'ord03_TA'
        'ord03_MagLS'
        'ord03_MagLSCC'
        'ord03_SpSub'
        'ord03_SpSubMod'
        'ord03_BiMagLS'
        'ref'
    };
    labels = {
        'Trunc'
        'EQ'
        'Tap'
        'TA'
        'MagLS'
        'MagLS+CC'
        'SpSub'
        'SpSubMod'
        'BiMagLS'
        'Reference'
    };
    n = numel(names);
    ndirs = numel(azHP);
    itd = zeros(ndirs,n);
    ild = zeros(ndirs,n);
    for i=1:n
        name = ['res_',names{i}];
        results = amt_cache('get', name, flags.cachemode); 
        itd(:,i) = results.itd;
        ild(:,i) = results.ild;
    end

    % Tight subplot  
    fig5size = [pagewidth 7.3025];
    fig5 = figure('units','centimeters','pos',[2 2 fig5size(1) fig5size(2)],...
        'Renderer','painters','PaperSize',[fig5size(1) fig5size(2)],...
        'paperposition',[0 0 fig5size(1) fig5size(2)]);
    gap = [.04 .06]; % gap between subplots in norm units (height width)
    marg_h = [0 .05]; % [.09 .1] figure height margins in norm units (lower upper)
    marg_w = [.03 .23]; % figure width margins in norm units (left right)
    [ha, ~] = tight_subplot(1,2,gap,marg_h,marg_w); 

    colors = parula(n+1);
    lsvec = {'-'};%,'-.'};
    lwvec = [0.5,0.5];
    mvec = {'^','v','x','s','o','d','p','h','>','+'};
    ms = 3; % marker size
    step_big = round(numel(azHP)/4);
    step_small = round(step_big/n);

    % ITD
    axes(ha(1))
    for i=1:n % multiply by 0.001 to have it in ms (cleaner plot)
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
        color = colors(i,:);
        if i==n % reference
            ls = ':';
            lw = 0.5;
            m = 'none';
            color = [0 0 0];
        end
        % For Matlab version >= R2016b:
        %polarplot(azHP,abs(itd(:,i))*0.001,'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on
        % For older versions:
        polarplot(azHP,abs(itd(:,i))*0.001,'Color',color,'LineWidth',lw,'LineStyle',ls); hold on
        polarplot(azHP(mi),abs(itd(mi,i))*0.001,'Color',color,'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end
    set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
    set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',{'0','45','90','135','','-135','-90','-45'},'RAxisLocation',-90)
    title('ITD (ms)')

    % ILD
    axes(ha(2))
    for i=1:n
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
        mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
        color = colors(i,:);
        if i==n % reference
            ls = ':';
            lw = 0.5;
            m = 'none';
            color = [0 0 0];
        end
        % For Matlab version >= R2016b:
        %polarplot(azHP,abs(ild(:,i)),'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on
        % For older versions:
        polarplot(azHP,abs(ild(:,i)),'Color',color,'LineWidth',lw,'LineStyle',ls); hold on
        polarplot(azHP(mi),abs(ild(mi,i)),'Color',color,'LineStyle','none','Marker',m,'MarkerSize',ms,'HandleVisibility','off')
    end
    set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
    set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',{'0','45','90','135','','-135','-90','-45'},'RAxisLocation',-90)
    title('ILD (dB)')
    legend(labels,'position',[0.8211    0.2543    0.1402    0.4362]);

    % Violins
    fig5bsize = [fig5size(1) fig5size(2)*0.55];
    fig5b = figure('units','centimeters','pos',[2 2 fig5bsize(1) fig5bsize(2)],...
        'Renderer','painters','PaperSize',[fig5bsize(1) fig5bsize(2)],...
        'paperposition',[0 0 fig5bsize(1) fig5bsize(2)]);
    gap = [.04 .06]; % gap between subplots in norm units (height width)
    marg_h = [0.3 0.03];
    marg_w = [0.06 0.01];
    [ha, ~] = tight_subplot(1,2,gap,marg_h,marg_w); 

    itderr = abs(itd-repmat(itd(:,end), 1, size(itd, 2)));
    ilderr = abs(ild-repmat(ild(:,end), 1, size(ild, 2)));

    axes(ha(1)), violinplot(itderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
    set(gca,'YTickLabelMode','auto')
    hold on, plot([0 n],[itdJND itdJND],'k:')
    set(ha(1),'XTickLabel',labels(1:end-1))
    ylabel('Abs. ITD error (ms)')
    set(gca,'fontsize',7,'xticklabelrotation',45)

    axes(ha(2)), violinplot(ilderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
    set(gca,'YTickLabelMode','auto')
    hold on, plot([0 n],[ildJND ildJND],'k:')
    set(ha(2),'XTickLabel',labels(1:end-1))
    ylabel('Abs. ILD error (dB)')
    set(gca,'fontsize',7,'xticklabelrotation',45)

end

%% Plot Fig. 6 (loudness per direction)

if flags.do_fig6    
    % Load data from file
    names = {
        'ref'
        'ord03_Trunc'
        'ord03_EQ'
        'ord03_Tap'
        'ord03_TA'
        'ord03_MagLS'
        'ord03_MagLSCC'
        'ord03_SpSub'
        'ord03_SpSubMod'
        'ord03_BiMagLS'
    };
    labels = {
        'Reference'
        'Trunc'
        'EQ'
        'Tap'
        'TA'
        'MagLS'
        'MagLS+CC'
        'SpSub'
        'SpSubMod'
        'BiMagLS'
    };
    n = numel(names);
    ndirs = numel(az);
    L = cell(n,1);
    PSD = zeros(ndirs,n);
    PSDavg = zeros(n,1);
    for i=1:n
        name = ['res_',names{i}];
        results = amt_cache('get', name, flags.cachemode); 
        L{i} = results.L;
        if isfield(results,'PSD')
            PSD(:,i) = results.PSD(:,:,1).'; % left ear only
            PSDavg(i) = results.PSDavg(:,:,1); % left ear only
        end
    end

    % Tight subplot
    fig6size = [pagewidth 6.9003];
    fig6 = figure('units','centimeters','pos',[2 2 fig6size(1) fig6size(2)],...
        'Renderer','painters','PaperSize',[fig6size(1) fig6size(2)],...
        'paperposition',[0 0 fig6size(1) fig6size(2)]);
    gap = [.07 .008]; % gap between subplots in norm units (height width)
    marg_h = [.13 .06]; % [.09 .1] figure height margins in norm units (lower upper)
    marg_w = [.06 .08]; % figure width margins in norm units (left right)
    [ha, ~] = tight_subplot(2,5,gap,marg_h,marg_w);
    clims = [1 5];

    colormap parula
    for i=1:n
        axes(ha(i)) 
        plotSph(az,el,L{i}(:,:,1)) % left ear only
        if i==1 || i==6
            set(ha(i),'YTick',-60:30:60)
            ylabel('Elevation (deg)')
        else
           set(ha(i),'YTickLabel',{})
        end
        if i>=6
            xlabel('Azimuth (deg)')
        else
          set(ha(i),'XTickLabel',{})
        end
        title(labels{i})

        set(gca,'fontsize',7,'XDir','reverse')
        grid(gca,'off')
        caxis(clims)
        if i==n
            c = colorbar; c.Label.String = 'Loudness (sones)';
            c.Position = [0.9373    0.3497    0.0111    0.3804];
        end
    end 
    for i=(n+1):10
        axis(ha(i),'off')
    end

    % Violins
    fig6bsize = [fig6size(1) fig6size(2)/2];
    fig6b = figure('units','centimeters','pos',[2 2 fig6bsize(1) fig6bsize(2)],...
        'Renderer','painters','PaperSize',[fig6bsize(1) fig6bsize(2)],...
        'paperposition',[0 0 fig6bsize(1) fig6bsize(2)]);

    violinplot(PSD(:,2:end),[],'ShowData',false,'BoxWidth',0.03); grid on
    set(gca,'XTickLabel',labels(2:end))
    ylabel('PSD (sones)')
    set(gca,'fontsize',7)

end

%% Plot Fig. 7 (models' outputs per spatial order)

if flags.do_fig7    
    % Load data from file
    names = {
        'Trunc'
        'EQ'
        'Tap'
        'TA'
        'MagLS'
        'MagLSCC'
        'SpSub'
        'SpSubMod'
        'BiMagLS'
        'ref'
    };
    labels = {
        'Trunc'
        'EQ'
        'Tap'
        'TA'
        'MagLS'
        'MagLS+CC'
        'SpSub'
        'SpSubMod'
        'BiMagLS'
        'Reference'
    };
    n = numel(names);
    m = numel(N_vec);
    PSD = nan(n,m);
    lat_prec = nan(n,m);
    pol_prec = nan(n,m);
    ext = nan(n,m);
    srm = nan(n,m);
    for i=1:n
        name = names{i};
        if strcmp(name,'ref')
            results = amt_cache('get', ['res_',name], flags.cachemode);
            PSD(i,:) = 0;
            lat_prec(i,:) = results.lat_prec;
            pol_prec(i,:) = results.pol_prec;
            ext(i,:) = mean(results.ext); 
            srm(i,:) = mean(results.srm);
        else
            for j=1:m
                N=N_vec(j);
                results = amt_cache('get', sprintf('res_ord%0.2d_%s',N,name), flags.cachemode); 
                PSD(i,j) = results.PSDavg(:,:,1); % left ear only
                lat_prec(i,j) = results.lat_prec;
                pol_prec(i,j) = results.pol_prec;
                ext(i,j) = mean(results.ext); 
                srm(i,j) = mean(results.srm);
            end
        end
    end
    % Tight subplot
    fig7size = [pagewidth 9.7790];
    fig7 = figure('units','centimeters','pos',[2 2 fig7size(1) fig7size(2)],...
        'Renderer','painters','PaperSize',[fig7size(1) fig7size(2)],...
        'paperposition',[0 0 fig7size(1) fig7size(2)]);
    gap = [.06 .05]; % gap between subplots in norm units (height width)
    marg_h = [.1 .03]; % [.09 .1] figure height margins in norm units (lower upper)
    marg_w = [.07 .02]; % figure width margins in norm units (left right)
    [ha, ~] = tight_subplot(2,3,gap,marg_h,marg_w); 

    colors = parula(n); 
    lsvec = {'-'};%,'-.'};
    mvec = {'^','v','x','s','o','d','p','h','>'};
    lwvec = [0.5,0.5];
    ms = 3; % marker size
    mi = [1:5:44]; % marker indices

    colors(end,:) = [0 0 0]; % last one is reference

    for i=1:n
        if i<n
            ls = lsvec{mod(i-1,numel(lsvec))+1};
            lw = lwvec(mod(i-1,numel(lwvec))+1);
            m = mvec{i};
        else
            ls = ':'; % reference
            lw = 0.5;
            m = 'None';
        end
        % For Matlab version >= R2016b:
        %plot(ha(1),PSD(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(1),'on')
        %plot(ha(2),lat_prec(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(2),'on')
        %plot(ha(3),pol_prec(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(3),'on')
        %plot(ha(4),ext(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(4),'on')
        %plot(ha(5),srm(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(5),'on')
        % For older versions:
        plot(ha(1),N_vec,PSD(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls),hold(ha(1),'on')
        plot(ha(2),N_vec,lat_prec(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls),hold(ha(2),'on')
        plot(ha(3),N_vec,pol_prec(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls),hold(ha(3),'on')
        plot(ha(4),N_vec,ext(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls),hold(ha(4),'on')
        plot(ha(5),N_vec,srm(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls),hold(ha(5),'on')
        plot(ha(1),N_vec(mi),PSD(i,mi),'Color',colors(i,:),'Marker',m,'MarkerSize',ms,'LineStyle','none','HandleVisibility','off')
        plot(ha(2),N_vec(mi),lat_prec(i,mi),'Color',colors(i,:),'Marker',m,'MarkerSize',ms,'LineStyle','none','HandleVisibility','off')
        plot(ha(3),N_vec(mi),pol_prec(i,mi),'Color',colors(i,:),'Marker',m,'MarkerSize',ms,'LineStyle','none','HandleVisibility','off')
        plot(ha(4),N_vec(mi),ext(i,mi),'Color',colors(i,:),'Marker',m,'MarkerSize',ms,'LineStyle','none','HandleVisibility','off')
        plot(ha(5),N_vec(mi),srm(i,mi),'Color',colors(i,:),'Marker',m,'MarkerSize',ms,'LineStyle','none','HandleVisibility','off')       
    end

    legend(ha(1),labels,'position',[0.7678 0.1414 0.1367 0.2944]);

    ylabel(ha(1),'PSD (sones)')
    ylabel(ha(2),'Lateral precision (deg)') 
    ylabel(ha(3),'Polar precision (deg)')
    ylabel(ha(4),'Externalisation')
    ylabel(ha(5),'SRM (dB)') 
    for i=1:5
        grid(ha(i),'on')
        xlim(ha(i),[1 44])
        set(ha(i),'XTick',[1,5:5:44])
        if i>3
            xlabel(ha(i),'Spatial order (N)')
        end
        set(ha(i),'fontsize',7)
    end
    set(ha(6),'visible','off') % hide last axis
    set(ha(1),'yscale','log') % set some plots' Y axis to log scale
    set(ha(2),'yscale','log','YTick',[1:5,7,10:10:50,70,100])
    set(ha(3),'yscale','log','YTick',[1:5,7,10:10:50,70,100])

end
 
end


