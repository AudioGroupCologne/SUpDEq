%DEMO_JELFS2011  Binaural speech intelligibility advantage
%
%   DEMO_JELFS2011 will plot the output of the
%   jelfs2011 binaural speech intelligibility advantage model for a
%   target azimuth angle of 0 deg. The masker position will move over a
%   full circle in the horizontal plane, and the output is visualized on
%   a polar plot. The KEMAR HRTF dataset is used.
%
%
%   See also: jelfs2011, culling2004
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_jelfs2011.php


%   #Author: Piotr Majdak (2015)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

target_azim=0;
database='kemar';

step = 5;
n_op = 360/step+1;
benefit = zeros(n_op,3);
weighted_SNR=zeros(n_op,3);
weighted_bmld=zeros(n_op,3);
angles = (0:step:360)'*pi/180;
for ii = 1:n_op
  [benefit(ii,:),weighted_SNR(ii,:), weighted_bmld(ii,:)] = ...
      jelfs2011({target_azim, database}, {(ii-1)*step, database});
end
figure; polar([angles angles angles], benefit); title('Benefit');
figure; polar([angles angles angles], weighted_SNR); title('Weighted SNR');
figure; polar([angles angles angles], weighted_bmld); title('Weighted BMLD');


