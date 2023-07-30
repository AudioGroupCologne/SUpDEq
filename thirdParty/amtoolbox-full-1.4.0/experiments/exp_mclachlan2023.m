function exp_mclachlan2023
%   
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_mclachlan2023.php

 
%   #Author: Glen McLachlan
%   #Requirements: circstat

subject = {'NH257','NH714','NH919','NH963','NH983','NH1016','NH1017'};

error_s=[]; error_y=[]; error_p=[];
for sub=1:length(subject)

    % load tracker data
    T_d=amt_load('mclachlan2023',[subject{sub} '_d_tracker.mat']); T_d=T_d.T;
    T_s=amt_load('mclachlan2023',[subject{sub} '_s_tracker.mat']); T_s=T_s.T;
    
    % load localisation data
    m_s=amt_load('mclachlan2023',[subject{sub} '_s_loca.mat']); m_s=m_s.m_s; %(:,:,1)=flat, (:,:,2)=full, (:,:,3)=VBAP 
    m_d=amt_load('mclachlan2023',[subject{sub} '_d_loca.mat']); m_d=m_d.m_d; %(:,:,1)=flat, (:,:,2)=full, (:,:,3)=fixed, (:,:,4)=VBAP

    %% omit data with incorrect movements
    for i=1:size(T_d,2)
        for j=1:size(T_d,1)
            tmpT = T_d{j,i};
            tmpM = m_d(j,:,i);
            if isempty(tmpT) % empty tracking data
               v(:) = 1;
               T_d{j,i} = nan;
               m_d(j,:,i) = nan;
            else % maximum rotation values: [1=start angle, 2=end angle, 3=max deviation from axis]
                if (tmpM(16)==0||tmpM(16)==180) %pitch rotation
                    v(1) = min(abs(tmpT(:,2)));     
                    v(2) = max(abs(tmpT(:,2)));
                    v(3) = max(abs(tmpT(:,1)));
                elseif (tmpM(16)==90||tmpM(16)==270) %yaw rotation
                    v(1) = min(abs(tmpT(:,1)));     
                    v(2) = max(abs(tmpT(:,1)));
                    v(3) = max(abs(tmpT(:,2)));
                end   
            end
    
            % check validity of rotations
            if (v(2)-v(1)<7 || v(3)>2 ) % minimum allowed rotation size and maximum allowed deviation
                T_d{j,i} = nan;
                m_d(j,:,i) = nan;
            end

            if ~any(isnan(T_d{j,i}),'all')
                % rotate in the opposite direction too much
                if ((tmpM(16)==0 && any(tmpT(:,2)<-2)) || (tmpM(16)==180 && any(tmpT(:,2)>2)) ||...
                        (tmpM(16)==90 && any(tmpT(:,1)>2)) || (tmpM(16)==270 && any(tmpT(:,1)<-2)))
                    T_d{j,i} = nan;
                    m_d(j,:,i) = nan;
                end
            end
        end
    %% divide into yaw and pitch data
    idx_y = m_d(:,16,i)==90|m_d(:,16,i)==270; % yaw rotations
    idx_p = m_d(:,16,i)==0|m_d(:,16,i)==180; % pitch rotations
    m_y = m_d(idx_y,:,i);
    m_p = m_d(idx_p,:,i);

    %% compute error metrics
    if i<4
    error_s = [error_s; computeError(m_s(:,:,i))];
    end
    error_y = [error_y; computeError(m_y)];
    error_p = [error_p; computeError(m_p)];

    T_y{sub,i} = T_d(idx_y,i);
    T_p{sub,i} = T_d(idx_p,i);
    end
end

    %% prepare input for box plots
    Subject = num2cell([repelem(1:7,3)';repelem(1:7,4)';repelem(1:7,4)']);
    Rot = [repmat({'Static'},size(error_s,1),1); repmat({'Yaw'},size(error_y,1),1);repmat({'Pitch'}, size(error_p,1),1)];
    Stim = [repmat({'Flat';'Full';'Free Field'},7,1);repmat({'Flat';'Full';'Frozen';'Free Field'},14,1)];
    error = array2table([error_s; error_y; error_p]);
    error.Properties.VariableNames ={'FBC', 'PrecP', 'PrecL'};
    boxinput = cell2table([Subject,Rot,Stim]);
    boxinput.Properties.VariableNames = {'Subject','Rotation', 'Stimulus'};
    boxinput = [boxinput, error];
    boxinput.Rotation = categorical(boxinput.Rotation,{'Static','Yaw','Pitch'});
    boxinput.Stimulus = categorical(boxinput.Stimulus,{'Free Field';'Flat';'Full';'Frozen'});
    boxinput.Subject = categorical(boxinput.Subject);


    %% Figure 3
    f = figure;
    f.Position = [100 100 640 500];
    
    subplot(2,2,1);
    b1 = boxchart(boxinput.Stimulus,boxinput.PrecL,'GroupByColor',boxinput.Rotation);
    ylabel('Lateral error (deg)')
    set(gca,'XTick',[])
    %ylim([0.1 25])
    xline( [1.5,2.5,3.5])
    title('(a)')

    % Create a tile on the right column to get its position
    subplot(2,2,2,'Visible','off');
    % Construct a Legend with the data from the sub-plots
    hL = legend(b1);
    % Move the legend to the position of the extra axes
    hL.Position = [0.66, 0.68, 0.20, 0.15];
    
    subplot(2,2,3);
    boxchart(boxinput.Stimulus,boxinput.PrecP,'GroupByColor',boxinput.Rotation);
    ylabel('Polar error (deg)')
    xline( [1.5,2.5,3.5])
    %ylim([10.1 55])
    title('(b)')

    subplot(2,2,4);
    boxchart(boxinput.Stimulus,boxinput.FBC,'GroupByColor',boxinput.Rotation);
    ylabel('FBC rate (%)')
    xline( [1.5,2.5,3.5])
    %ylim([0.1 50])
    title('(c)')

    set(subplot(2,2,1), 'Position', [0.1, 0.56, 0.40, 0.40])
    set(subplot(2,2,3), 'Position', [0.1, 0.1, 0.40, 0.40])
    set(subplot(2,2,4), 'Position', [0.58, 0.1, 0.40, 0.40])

    %% Figure 4
    f= figure;
    f.Position = [100 100 640 500];
    
    subplot(2,2,1);
    b1 = boxchart(boxinput.Rotation,boxinput.PrecL,'GroupByColor',boxinput.Stimulus);
    ylabel(['Lateral error (deg)'])
    set(gca,'XTick',[])
    %ylim([0.1 55])
    xline( [1.5,2.5])
    title('(a)')

    % Create a tile on the right column to get its position
    subplot(2,2,2,'Visible','off');
    % Construct a Legend with the data from the sub-plots
    hL = legend(b1);
    % Move the legend to the position of the extra axes
    hL.Position = [0.66, 0.68, 0.20, 0.15];
    
    b2 = subplot(2,2,3);
    boxchart(boxinput.Rotation,boxinput.PrecP,'GroupByColor',boxinput.Stimulus);
    ylabel(['Polar error (deg)'])
    xline( [1.5,2.5])
    %set(gca,'XTick',[])
    %ylim([0.1 55])
    title('(b)')

    b3 = subplot(2,2,4);
    boxchart(boxinput.Rotation,boxinput.FBC,'GroupByColor',boxinput.Stimulus);
    ylabel('FBC rate (%)')
    xline( [1.5,2.5])
    %ylim([0.1 50])
    title('(c)')

    set(subplot(2,2,1), 'Position', [0.1, 0.56, 0.40, 0.40])
    set(subplot(2,2,3), 'Position', [0.1, 0.1, 0.40, 0.40])
    set(subplot(2,2,4), 'Position', [0.58, 0.1, 0.40, 0.40])

end

function errors= computeError(m)
maxlat=60; 

idx_fb = find(acos(sum(m(:,9:11).*m(:,12:14),2))>... % fb reversals
acos(sum([-m(:,9) m(:,10) m(:,11)].*m(:,12:14),2))...
& (abs(m(:,8))<80 | (m(:,8)>100 & m(:,8)<260))); % drop dirs near midline

idx_ud = find(acos(sum(m(:,9:11).*m(:,12:14),2))>... % ud reversals
acos(sum([m(:,9) m(:,10) -m(:,11)].*m(:,12:14),2))...
& (abs(m(:,8))>10 & (m(:,8)<170 | m(:,8)>190))); % drop dirs near midline

idx_ms = find(abs(m(:,7))<=maxlat & abs(m(:,5))<=maxlat) ; % only take subset closer to midsagittal
idx_fbm = intersect(idx_ms,idx_fb);
idx_udm = intersect(idx_ms,idx_ud); 

% lateral accuracy and precision
accL=mean(m(:,7)-m(:,5));  
precL=sqrt(sum((m(:,7)-m(:,5)-accL).^2)/size(m,1));

FBC = size(idx_fbm,1)*100/size(idx_ms,1); %percent fb reversals

idx = ismember(1:size(m,1),[idx_fbm;idx_udm]); %remove fb and ud errors
tmp = m(~idx,:);
idx = find(abs(tmp(:,7))<=maxlat & abs(tmp(:,5))<=maxlat); % only take subset closer to midsagittal
precP=rad2deg(circ_std(deg2rad(tmp(idx,8)-tmp(idx,6))));

errors=[FBC,precP,precL];

end
