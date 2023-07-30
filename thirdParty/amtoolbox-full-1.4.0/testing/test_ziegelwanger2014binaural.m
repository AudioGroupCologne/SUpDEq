% Tests if ziegelwanger2014 yields different model parameters for different subjects
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/testing/test_ziegelwanger2014binaural.php


if ~exist('amt_stop','file'), amt_start; end

S=data_baumgartner2014;

Xm=nan(length(S),6); 
R=nan(length(S),2);
Tau=nan(length(S),2);
Phi=nan(length(S),2);
Theta=nan(length(S),2);
ID=nan(length(S),1);

for ii=1:length(S)
  
  [~,res]=ziegelwanger2014(S(ii).Obj,4,1,1); % use thr method on lp-filtered IRs, AMT 1.1.1
  id=S(ii).id;
  ID(ii)=str2num(id(3:end));
    % left ear
  p=res.p_offaxis(:,1);
  R(ii,1)=p(1)*100; % in cm
  Xm(ii,1:3)=p(2:4)*100; % in cm 
  Tau(ii,1)=p(5)*1000; % in ms
  Phi(ii,1)=rad2deg(p(6)); % in degrees
  Theta(ii,1)=rad2deg(p(7)); % in degrees  
  leftEar = sprintf('  Left: r=%0.3g cm; M=[%5.3g, %5.3g, %5.3g] cm; tau=%5.3g ms; phi=%5.3g�; theta=%5.3g�\n', ...
    p(1)*100, p(2)*100, p(3)*100, p(4)*100, p(5)*1000, rad2deg(p(6)), rad2deg(p(7)));
  
    % right ear
  p=res.p_offaxis(:,2);
  R(ii,2)=p(1)*100; % in cm
  Xm(ii,4:6)=p(2:4)*100; % in cm 
  Tau(ii,2)=p(5)*1000; % in ms
  Phi(ii,2)=rad2deg(p(6)); % in degrees
  Theta(ii,2)=rad2deg(p(7)); % in degrees
  rightEar = sprintf(' Right: r=%0.3g cm; M=[%5.3g, %5.3g, %5.3g] cm; tau=%5.3g ms; phi=%5.3g�; theta=%5.3g�\n', ...
    p(1)*100, p(2)*100, p(3)*100, p(4)*100, p(5)*1000, rad2deg(p(6)), rad2deg(p(7)));
    % display
  TOAcheck = ['***' 10 id ': ' leftEar id ': ' rightEar];
  disp(TOAcheck);

  S(ii).Obj = SOFAaddVariable(S(ii).Obj,'TOAModel','SR',res.p_offaxis);
  S(ii).Obj = SOFAupdateDimensions(S(ii).Obj);
end

%% Plot results
figure; 
plot(ID,R,'o'); 
title('Head Radius as a function of subject ID'); 
legend('Estimated from left-ear HRTFs', 'Estimated from right-ear HRTFs');
ylabel('Head radius (cm)');
xlabel('Subject ID');

figure; 
plot(ID,Xm,'o'); 
title('Coordinates of the head center as a function of subject ID'); 
legend('Left X', 'Left Y', 'Left Z', 'Right X', 'Right Y', 'Right Z');
ylabel('Coordinate (cm)');
xlabel('Subject ID');

figure; 
plot(ID,Phi,'o'); 
title('Ear position on the head'); 
legend('Left ear', 'Right ear');
ylabel('Azimuth angle (degrees)');
xlabel('Subject ID');

figure; 
plot(ID,Theta,'o'); 
title('Ear position on the head'); 
legend('Left ear', 'Right ear');
ylabel('Elevation angle (degrees)');
xlabel('Subject ID');

