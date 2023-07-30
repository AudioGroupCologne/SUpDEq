%DEMO_BAUMGARTNER2014_BLOCKPROCESSING Demo for sagittal-plane localization model from Baumgartner et al. (2014)
%
%   DEMO_BAUMGARTNER2014_BLOCKPROCESSING demonstrates how to compute and visualize 
%   the baseline prediction (localizing broadband sounds with own ears) 
%   for a listener of the listener pool and the median plane using the 
%   sagittal-plane localization model from Baumgartner et al. (2014) within a
%   blockprocessing framework
%
%   Figure 1: Baseline prediction
% 
%      This demo computes the baseline prediction (localizing broadband 
%      sounds with own ears) for an exemplary listener (NH58).
%
%      Predicted polar response angle probability of subject NH58 as a  
%      function of the polar target angle with probabilities encoded by
%      brigthness.
%
%   See also: baumgartner2014 exp_baumgartner2014 baumgartner2014_virtualexp
%   localizationerror demo_baumgartner2014
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_baumgartner2014_blockprocessing.php


%   #Author: Fabian Brinkmann (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% Load HRIRs --------------------------------------------------------------
% HRIRs from the FABIAN head-related transfer function database
hrirname = 'FABIAN_HRIR_measured_HATO_0';
HRIRs = amt_load('engel2021',[hrirname,'.sofa']);


% use only median plane HRIRs for faster computation

id = HRIRs.SourcePosition(:,1) == 0 | HRIRs.SourcePosition(:,1) == 180;

HRIRs.Data.IR = HRIRs.Data.IR(id,:,:);

HRIRs.SourcePosition = HRIRs.SourcePosition(id,:); 

HRIRs = SOFAupdateDimensions(HRIRs);



% Parameters for blockwise processing -------------------------------------

% audio content. We use noise  and low-pass filtered noise for this demo.

% In other cases this might be raw and processed audio files.

fsstim = 44100;

content_raw = randn(2^13, 1);

% 4th order Butterworth low-pass with 3 kHz cut-off frequency

sos = [1.26277170e-03  2.52554339e-03  1.26277170e-03 1.00000000e+00 -1.31605254e+00  4.46155785e-01
       1.00000000e+00  2.00000000e+00  1.00000000e+00 1.00000000e+00 -1.57087561e+00  7.26170328e-01];

if ~isoctave
    content_proc = sosfilt(sos, content_raw);
else
    [b, a] = butter(4, 3000 * 2/fsstim);
    content_proc = filter(b, a, content_raw);
end

N_block = 2048; % block length in samples

N_hop = 2048/2; % hop size in samples

window = true;  % apply Hann window to each block



% Parameters for baumgartner 2014 -----------------------------------------

errflag = 'QE_PE_EB';      % (currently only working with this errflag)

regular_flag = 'regular';

S = 0.76;

polsamp = -30:5:210;

fs = HRIRs.Data.SamplingRate;



%% Run baumgartner 2014 ---------------------------------------------------

% check equal length of audio content
if size(content_raw, 1) ~= size(content_proc, 1)
    error('audio contents must be of the same length')
end

% check equal number of channels of audio content
if size(content_raw, 2) ~= size(content_proc, 2)
    error('audio contents must have the same number of channels')
end

% check number of channels of audio content
if size(content_raw, 2) > 2
    error('audio contents must have 1 or 2 channels')
end

if ~isstruct(HRIRs)
    error('HRIRs must be a SOFA file')
end

if ~isfield(HRIRs, 'GLOBAL_SOFAConventions')
    error('HRIRs must be a SOFA file')
end

if ~strcmp(HRIRs.GLOBAL_SOFAConventions, 'SimpleFreeFieldHRIR')
    error('HRIRs must be a SOFA file of the ''SimpleFreeFieldHRIR'' convention')
end

if N_block > size(content_raw, 1)
    error('Block lengths exceeds length of audio content')
end

if N_block < size(HRIRs.Data.IR, 3)
    error('N_block is smaller than the HRIR length')
end

if N_block < N_hop
    warning('Hop size exceeds block length')
end


%% loop across chunks of audio content ------------------------------------

% get the window function
if window
    w = hann(N_block);
end

% audio content length and number of channels
N_samples = size(content_raw, 1);
N_channels = size(content_raw, 2);

% number of blocks
if N_block == N_samples
    N_blocks = 1;
else
    N_blocks = ceil((N_samples-N_block) / N_hop) + 1;
end

% zero pad audio content
N_samples = N_block + (N_blocks - 1) * N_hop;
content_raw(end+1:N_samples, :) = 0;
content_proc(end+1:N_samples, :) = 0;

% force two channel audio content for convenience
if N_channels == 1
    content_raw = [content_raw content_raw];
    content_proc = [content_proc content_proc];
end

% copy SOFA file
H_raw = HRIRs;
H_proc = HRIRs;

% get HRIRs
h = shiftdim(HRIRs.Data.IR, 2);

% zero pad for easy overlap and add fft convolution
N_HRIR = size(h, 1);
h(end+1:N_block+N_HRIR, :, :) = 0;


% time axis for blocks
t = nan(N_blocks, 1);

% arrays for overlap and add
ola_target = zeros(N_HRIR, size(h,2), size(h,3));
ola_template = ola_target;



N_wait = ceil(N_blocks/100);
amt_disp('Processing audio blocks' );
for nn = 1:N_blocks


    amt_disp(['block nr: ', num2str(nn)] ,'volatile');
    % get current block of audio
    nn_start = (nn-1) * N_hop + 1;
    nn_end = nn_start + N_block - 1;
    raw = content_raw(nn_start:nn_end, :, :);
    proc = content_proc(nn_start:nn_end, :, :);
   
    % time of current block
    t(nn) = mean([nn_start nn_end]) / HRIRs.Data.SamplingRate;
    
    % zero pad for easy overlap and add fft convolution
    raw(end+1:N_block+N_HRIR, :) = 0;
    proc(end+1:N_block+N_HRIR, :) = 0;
    
    % convolve HRIRs with audio content
    raw_l = zeros(size(h(:,:,1)));
    raw_r = zeros(size(h(:,:,1)));
    proc_l = zeros(size(h(:,:,1)));
    proc_r = zeros(size(h(:,:,1)));
    for ii = 1:size(h, 2)
        raw_l(:,ii) = fft(raw(:,1)) .* fft(h(:,ii,1));
        raw_r(:,ii) = fft(raw(:,2)) .* fft(h(:,ii,2));
        proc_l(:,ii) = fft(proc(:,1)) .* fft(h(:,ii,1));
        proc_r(:,ii) = fft(proc(:,2)) .* fft(h(:,ii,2));
    end
    
    

    % join left and right channels
    if ~isoctave
        raw = cat(3, ifft(raw_l, 'symmetric'), ifft(raw_r, 'symmetric'));
        proc = cat(3, ifft(proc_l, 'symmetric'), ifft(proc_r, 'symmetric'));
    else
        raw = cat(3, ifft(raw_l), ifft(raw_r));
        proc = cat(3, ifft(proc_l), ifft(proc_r));
    end
    
    % overlap and add (checked and working :)
    raw(1:N_HRIR, :, :) = raw(1:N_HRIR, :, :) + ola_target;
    proc(1:N_HRIR, :, :) = proc(1:N_HRIR, :, :) + ola_template;
    
    % save overlap for next step
    ola_target = raw(N_block+1:end, :, :);
    ola_template = proc(N_block+1:end, :, :);
    
    % remove overlap from current block
    raw = raw(1:N_block, :, :);
    proc = proc(1:N_block, :, :);   

    % window
    if window
        for jj=1:size(raw, 2)
            raw(:,jj, 1) = raw(:,jj,2) .* w;
            proc(:,jj, 1) = proc(:,jj,2) .* w;
            raw(:,jj, 2) = raw(:,jj,2) .* w;
            proc(:,jj, 2) = proc(:,jj,2) .* w;
        end
    end    

    % put in SOFA objects
    H_raw.Data.IR = shiftdim(raw, 1);
    H_proc.Data.IR = shiftdim(proc, 1);
    
    % update dimensions of SOFA files
    if nn == 1
        H_raw = SOFAupdateDimensions(H_raw);
        H_proc = SOFAupdateDimensions(H_proc);
    end
    
    % run the loaclization model with raw and processed audio content
    [err_current, pred_current] = baumgartner2014(H_proc, HRIRs, errflag, 'S', S, 'polsamp', polsamp, 'fs', fs, 'fsstim', fsstim);
    [err_base, pred_base] = baumgartner2014(H_raw, HRIRs, errflag, 'S', S, 'polsamp', polsamp, 'fs', fs, 'fsstim', fsstim );    

    % collect the output

        % allocate space
        if nn == 1
            clear err
            err.qe = nan(N_blocks, 1);
            err.pe = err.qe;
            err.pb = err.qe;            

            pred.p = nan(size(pred_current.p,1), size(pred_current.p,2), N_blocks);
            pred.rang = pred_current.rang;
            pred.tang = pred_current.tang;
            
            err.qe_baseline = err.qe;
            err.pe_baseline = err.qe;
            err.pb_baseline = err.qe;
            pred.p_baseline = pred.p;
        end        

        % save current results
        err.qe(nn) = err_current.qe;
        err.pe(nn) = err_current.pe;
        err.pb(nn) = err_current.pb;
        
        pred.p(:,:,nn) = pred_current.p;
        
        err.qe_baseline(nn) = err_base.qe;
        err.pe_baseline(nn) = err_base.pe;
        err.pb_baseline(nn) = err_base.pb;
        pred.p_baseline(:,:,nn) = pred_base.p;
 
end
amt_disp();
    % add time vector to the output
    err.t = t;
    pred.t = t;    

    varargout{1} = err;
    varargout{2} = pred;





%% Plot results -----------------------------------------------------------



% time axis for audio content

t = (0:size(content_raw,1) - 1) / fs;



figure()



subplot(3,1,1)

plot(t, content_raw)

title 'Template audio content'

xlabel 'Time in seconds'

ylabel 'Amplitude'

xlim([0 t(end)])

if size(content_raw, 2) > 1

    legend('left', 'right', 'Location', 'SouthEast')

end



subplot(3,1,2)

plot(err.t, err.pe, 'r.-')

hold on

plot(err.t, err.pe_baseline, 'k.-')

title 'Polar error'

xlabel 'Time in seconds'

ylabel({'Polar error' 'in degrees'})

xlim([0 t(end)])

ylim([0 1.1 * max([err.pe; err.pe_baseline])])



subplot(3,1,3)

plot(err.t, err.qe, 'r.-')

hold on

plot(err.t, err.qe_baseline, 'k.-')

legend('Processing', 'Baseline', 'location', 'best')

title 'Quadrant error'

xlabel 'Time in seconds'

ylabel({'Quadrant error' 'in percent'})

xlim([0 t(end)])

ylim([0 1.1 * max([err.qe; err.qe_baseline])])




