%% SARITA - Spherical Array Interpolation by Time Alignment
%
% function [drirs_upsampled, DRTFs_ups] = Sarita_upsampling(drirs, source_grid, target_grid, radius, varargin)
%
% Perform spatial upsampling of sparse spherical microphone array signals
% to dense spherical microphone array singals
% A detailed description of the algorithm is presented in [1]
%
% Parameters:
% ------------
%   drirs:         array irs 
%   source_grid:   array sampling grid of the array to be upsampled in RAD
%                  [M X (az, col)] M: number of sampling points, 
%                                az: horizontal angle ranging from 0 to 2pi,
%                                col: colatitude angle ranging from 0 to pi
%   target_grid:   target sampling grid. Same convention as source_grid
%   radius:        array radius     
%   f_transit:     {default: 0}  
%   frame_length:  {default: 32}
%   frame_overlap: {default: frame_length/2}
%   c:             {default: 340}
%   fs:            {default: 48 kHz}
%
% Returns:
% -----------
%
%   drirs_upsampled: upsampled array signals in time domain
%   DRTFs_upsampled: upsampled array signals in frequency domain (single-sided)
%
% Dependencies:
% -------------
%
% SFS Toolbox >=2.5.0 https://sfs-matlab.readthedocs.io/en/2.5.0/#              
%
% References:
% -------------
%
% [1] 
%
%  (C)   2020 CP,  Christoph Pörschmann
%        TH Köln - University of Applied Sciences
%        Institute of Communications Engineering
%        Department of Acoustics and Audio Signal Processing
% 
%
% versioning:
%
% 03.11.2020 (TL):
% 28.09.2021 (CP): performance improvements: change loop order(frame-based/channel based. 
%                                            precalculate neighbor ids and weights 
% 22.12.2021 (CP): avoid negative time shifts
% 30.12.2021 (CP): add frame_overlap parameter required for streaming-based processing
% 03.01.2021 (TL): clean up, remove struct, implement varagin parameter
% 07.01.2021 (TL): pass SMA signals in time domain and coordinates in radiants
% 31.01.2021 (TL): ---- version for listening experiment -----
% 07.04.2021 (TL): refactoring

function [drirs_upsampled, DRTFs_ups] = Sarita_upsampling(drirs, source_grid, target_grid, radius, varargin)
% check inputs
if nargin < 4
    error('Not enough input arguments');
end
input_params = {'f_transit', 'frame_length', 'frame_overlap', 'c', 'fs'}; 

default.f_transit     = 0;
default.frame_length  = 32; 
default.frame_overlap = 0;
default.c             = 340;
default.fs            = 48e3;

if mod(size(varargin, 2), 2)
    error('Check passed parameter value pairs.')
end

% check inputs
for p_IDX = 1 : 2 : size(varargin, 2)-1
    is_valid = 0;
    for d_IDX = 1:size(input_params, 2)
        if strcmp(input_params{d_IDX}, varargin{p_IDX})
            [~] = evalc([input_params{d_IDX}, '=', 'varargin{p_IDX+1}']);
            is_valid = 1;
        end
    end
    if ~is_valid
        warning('%s is not valid parameter. It will be ignored ', varargin{p_IDX});
    end
end 

% fill all params which were not passed with default values 
for d_IDX = 1:size(input_params, 2)
    value = ['default.' input_params{d_IDX}];
    if ~exist(input_params{d_IDX}, 'var')
        [~] = evalc([input_params{d_IDX}, '=', value]); 
    end
end

frame_length = floor(frame_length/2)*2; % force frameLength to be even
if frame_overlap == 0 || frame_overlap > frame_length/2
    frame_overlap = frame_length/2;
end
clear default varargin

%% end checking arguments 
%disp('SARITA -- perform upsampling of array impulse responses') %JMA-Edit: Suppress prints
hannsize = frame_overlap*2;
hannwin = hann(hannsize).';
win = [hannwin(1:length(hannwin)/2), ...
       ones(1,frame_length-length(hannwin)), ...
       hannwin(length(hannwin)/2+1:length(hannwin))];

[sparse_grid_cart(:, 1), sparse_grid_cart(:, 2), sparse_grid_cart(:, 3)] = ... 
    sph2cart(source_grid(:, 1), pi/2 - source_grid(:, 2), ones(length(source_grid(:, 2)), 1));    


% allocate some memory

num_neighbors_dense_grid       = zeros(length(target_grid));    % Number of nearest neighbors for each sampling point
idx_neighbors_dense_grid       = zeros(length(target_grid), 4); % Indices of neighbors of each sampling point (max 4)
weights_neighbors_dense_grid   = zeros(length(target_grid), 4); % Weights of neighbors of each sampling point 
neighbors_combinations         = zeros(2, 0);                   % Array containing all combinations of nearest neighbors
combination_ptr                = zeros(2, 0);                   % Array describing which neighborsCombination is required for each the cross correlations   
 
for dense_idx = 1:size(target_grid, 1)         
    % find nearest neighbors and calculate interpolation weights
    [x, y, z] = sph2cart(target_grid(dense_idx, 1), ...
                         pi/2 - target_grid(dense_idx, 2), ...
                         1);
    [neighbor_idx, weights] = findvoronoi(sparse_grid_cart, [x, y, z]);
  
    % store them
    num_neighbors_dense_grid(dense_idx) = length(neighbor_idx);
    idx_neighbors_dense_grid(dense_idx, 1:num_neighbors_dense_grid(dense_idx)) = neighbor_idx;
    weights_neighbors_dense_grid(dense_idx, 1:num_neighbors_dense_grid(dense_idx)) = weights;
       
    % calculate the maximum timeshifts between all next neighbors
    % candidates and all combination next neighbors
    for sparse_node_idx = 2:length(neighbor_idx)   
        angle = acos(dot(sparse_grid_cart(neighbor_idx(1), :), sparse_grid_cart(neighbor_idx(sparse_node_idx), :))); % get angle between sampling points
        max_shifts(dense_idx, sparse_node_idx-1) = abs(ceil(angle * radius*fs/c)); % calculate the maximum shift                        

        % Determine combinations of nearest neighbors
        compareEntries = ismember(neighbors_combinations, [neighbor_idx(1) neighbor_idx(sparse_node_idx)]');
        if ~isempty(compareEntries)            
            if (max(compareEntries(1, :) .* compareEntries(2, :))  == 0)                 
                % entry does not exist
                neighbors_combinations = [neighbors_combinations, [neighbor_idx(1) neighbor_idx(sparse_node_idx)]'];
                combination_ptr = [combination_ptr, [length(neighbors_combinations(1,:)) 1]'];
            else
                % entry already exist, just store the correct reference in
                % combination_ptr
               [~, position]=max(compareEntries(1, :).*compareEntries(2, :));
                if neighbors_combinations(:, position) == [neighbor_idx(1), neighbor_idx(sparse_node_idx)]' 
                    % entry exists in the same order 
                    combination_ptr = [combination_ptr, [position 1]'];
                else
                    % entry exists in the inverted order
                    combination_ptr = [combination_ptr, [position -1]'];
                end
            end
        else
            % no element in list yet. First element is put in list
            neighbors_combinations = [neighbors_combinations [neighbor_idx(1) neighbor_idx(sparse_node_idx)]'];
            combination_ptr = [combination_ptr [length(neighbors_combinations(1,:)) 1]'];
        end
        
    end
end
      
abs_max_shift = max(max(abs(max_shifts))); % get the maximal possible time shift, required to increase buffer at the end

% zeropad with 2 x maxShiftOverAllDir samples
drirs = [drirs, zeros(size(drirs, 1), abs_max_shift)]; 

% allocate memory for upsampled drirs
drirs_upsampled = zeros(length(target_grid), length(drirs(1, :)) + abs_max_shift*2);

% calculate frequency delay line
frequencies = linspace(0, fs, frame_length+abs_max_shift*2); 
freq_delay_line = -1i * 2*pi * frequencies/fs; 

num_frames = floor(length(drirs(neighbor_idx, :)) / (frame_length-frame_overlap))-1;

for frame_idx = 1:num_frames % loop over all frames                       
    startTab = (frame_idx-1)*(frame_length-frame_overlap)+1;
    endTab = startTab + frame_length - 1;
    irsFrame = repmat(win, [size(drirs, 1), 1]) .* drirs(:, startTab:endTab); % make compatible with older Matlabversions with repmat
    
    % calculate cross-correlation between each nn frame
    for cc_idx = 1:length(neighbors_combinations)  
         corr_frame(:, cc_idx) = xcorr(irsFrame(neighbors_combinations(1, cc_idx), :), ...
                                       irsFrame(neighbors_combinations(2, cc_idx), :));
    end   
    
    neighbor_idx_cnt = 0; % index of entry in combination_ptr 
    for dense_node_index = 1:size(target_grid, 1) % loop over all desired sampling points        
        time_shift_mean = 0;
        time_shift_current = 0;
        
        % get nearst neighbors, weights, and maxShift of actual sampling point, can
        % later be directly addressed in following lines if desired
        neighbor_idx = idx_neighbors_dense_grid(dense_node_index, 1:num_neighbors_dense_grid(dense_node_index));
        weights = weights_neighbors_dense_grid(dense_node_index, 1:num_neighbors_dense_grid(dense_node_index));
        maxShift = max_shifts(dense_node_index, 1:num_neighbors_dense_grid(dense_node_index)-1);            
        neighborsIRs = irsFrame(neighbor_idx, :); % get all next neighbour irs of one frame and perform windowing              
          
        for sparse_node_idx = 2:length(neighbor_idx) % loop over all next neighbours            
            neighbor_idx_cnt = neighbor_idx_cnt+1;
            correlation = corr_frame(:, combination_ptr(1, neighbor_idx_cnt));                    
            if combination_ptr(2, neighbor_idx_cnt) == -1
                correlation = correlation(end:-1:1);
            end   
                
            % look for maximal value in the crosscorrelated IRs only in the relevant area                            
            correlation = correlation(frame_length-maxShift(sparse_node_idx-1):frame_length+maxShift(sparse_node_idx-1));
            [~, maxpos] = max(correlation);
            if (and(maxpos > 1, maxpos < length(correlation))) % Make sure that the determined position is not the first or last value of the correlation
                % refine the value of maxpos by calculating
                % the peak of a parabel from the maxpos and the adjacant values
                maxpos = calcParabelPeak([maxpos-1, maxpos, maxpos+1], ...
                                         [correlation(maxpos-1), correlation(maxpos), correlation(maxpos+1)]);
            end
            % determine required shift between the nearest neighbors     
            time_shift_current(sparse_node_idx) = (maxpos - (length(correlation) + 1)/2);            
            time_shift_mean = time_shift_mean+time_shift_current(sparse_node_idx) * weights(sparse_node_idx);
        end
        clear correlation
        
        % align every block according to the calculated time shift, weight
        % in amplitude and sum up
        for sparse_node_idx = 1:length(neighbor_idx)
            current_block = neighborsIRs(sparse_node_idx, :) * weights(sparse_node_idx);
            
            % add zeros at the end of the block to prevent wraparounds
            current_block = [current_block, zeros(1, 2 * abs_max_shift)];           
            time_shift_final = -time_shift_mean + time_shift_current(sparse_node_idx) + abs_max_shift;   % As abs_max_shift is added, time_shift_final will always be positive                      
            if time_shift_final < 0 % avoid negative time shifts
                time_shift_final = 0; 
            end
            
            % apply time shift (approx. 40% of the performace).
            current_block = fft(current_block);
            current_block = current_block .* exp(freq_delay_line .* time_shift_final);
            current_block = real(ifft(current_block, 'symmetric')); % inverse dft    
            
            % add to output
            drirs_upsampled(dense_node_index, startTab:endTab + 2*abs_max_shift) = ...
                drirs_upsampled(dense_node_index,startTab:endTab + 2*abs_max_shift) + current_block;
        end
        clear current_block 
    end
end

% kill maxShiftAllDir leading zeros to compensate the delay. 
% kill 2 * maxShiftOverAllDir at the tail to ensure in- and output have
% same length
drirs_upsampled = drirs_upsampled(:, abs_max_shift+1 : length(drirs_upsampled(1, :)) - 2*abs_max_shift);

if f_transit > 0 % perform low-frequency extension. Combine spectra of upsampled DRIRs and
                 % DRIRs interpolated with voronoi weights on complex data
                 
    % determine frequency bin of transition
    f_vector = linspace(5*eps, fs/2, round(size(drirs_upsampled, 2)/+1));
    [~, f_transit_index] = min(abs(f_vector - f_transit));
    
    % add leading zeros according to the number of bins of transition 
    % to avoid wrap around due to non ideal brick-wall filtering
    num_tabs = size(drirs_upsampled, 2); % store original length
    drirs_upsampled = [zeros(size(drirs_upsampled, 1), f_transit_index+1), drirs_upsampled];
    drirs = [zeros(size(drirs, 1), f_transit_index+1), drirs];
    
    NFFT = length(drirs_upsampled)*2;
    
    % DFT of upsampled data
    DRTFs_ups = fft(drirs_upsampled, NFFT, 2);
    DRTFs_ups = DRTFs_ups(:, 1:round(end/2)+1);

    % DFT of measured array signals 
    DRTFs = fft(drirs, NFFT, 2);
    DRTFs = DRTFs(:, 1:round(end/2)+1);
    
    DRTFsInterpLin = zeros(size(DRTFs_ups));
    for dense_idx = 1:size(target_grid, 1)
        neighbor_idx = idx_neighbors_dense_grid(dense_idx, 1:num_neighbors_dense_grid(dense_idx));
        weights = weights_neighbors_dense_grid(dense_idx, 1:num_neighbors_dense_grid(dense_idx));    
        
        % linear complex interpolation below transition frequency
        for nn_idx = 1:length(neighbor_idx)
            DRTFsInterpLin(dense_idx, :) = DRTFsInterpLin(dense_idx, :) + ...
                                           weights(nn_idx) * DRTFs(neighbor_idx(nn_idx), :);   
        end 
    end
    
    % replace low_frequency part with linear interpolated data
    DRTFs_ups(:, 1:f_transit_index-1) = DRTFsInterpLin(:, 1:f_transit_index-1);
    
    drirs_upsampled = ifft([DRTFs_ups, conj(DRTFs_ups(:, end-1:-1:2))], [], 2, 'symmetric');
    drirs_upsampled = circshift(drirs_upsampled, -(f_transit_index+1), 2);
    drirs_upsampled = drirs_upsampled(:, 1:num_tabs);
    % window
    win = hann(33);
    drirs_upsampled(:, 1:ceil(length(win)/2)) = drirs_upsampled(:, 1:ceil(length(win)/2)) ...
                                                .* win(1:ceil(end/2)).';
    drirs_upsampled(:, end-floor(length(win)/2):end) = drirs_upsampled(:, end-floor(length(win)/2):end) ...
                                                       .* win(end-floor(length(win)/2):end).';
end
% compute dft of final drirs
DRTFs_ups = fft(drirs_upsampled, [], 2);
DRTFs_ups = DRTFs_ups(:, 1:round(end/2)+1);

%disp('... done') %JMA-Edit: Suppress prints
end 

function xp = calcParabelPeak(x, y)
    % this function determines the peak of the parabel described by 3 points
    % reference: https://www.arndt-bruenner.de/mathe/10/parabeldurchdreipunkte.html
    xp = (x(2)^2*(y(3) - y(1)) - x(1)^2*(y(3) - y(2)) - x(3)^2*(y(2) - y(1)))/(2*(x(2)*(y(3) - y(1)) - x(1)*(y(3) - y(2)) - x(3)*(y(2) - y(1))));
end 
