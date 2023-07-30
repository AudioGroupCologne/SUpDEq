function [outsig_fine,fc_fine,outsig_env,fc_env,outsig_ild] = dietz2011_filterbank(insig,fs,fc,varargin)
%DIETZ2011_FILTERBANK  filterbank of Dietz 2011 binaural model
%   Usage: [...] = dietz2011_filterbank(insig,fs,fc);
%
%   Input parameters:
%       insig       : binaural signal for which values should be calculated
%       fs          : sampling rate (Hz)
%       fc          : center frequencies of gammatone filterbank
%
%   Output parameters:
%       outsig_fine : output signal of fine structure filter
%       fc_fine     : center frequencies processed with fine structure filter
%       outsig_env  : output signal of envelope filter
%       fc_env      : center frequencies processed with modulation filter
%       outsig_ild  : output signal of ILD filter
%
%   DIETZ2011_FILTERBANK(insig,fs,fc) filters all frequency channels
%   of insig with a modulation filterbank consisting of three filters.
%   One centered at the center frequencies for the fine structure of the
%   signals. One centered at a fixed frequency of 135 Hz for the envelope of
%   the signals. And one a lowpass filter with a cutoff frequency of 30 Hz
%   for the calculation of the interaural level difference.
%
%   DIETZ2011_FILTERBANK accepts the following optional parameters:
%
%     'filter_order',fo
%                    Filter order for the two gammatone filter used for the fine
%                    structure and envelope of the modulation filter bank. The
%                    default value is 2.
%
%     'filter_attenuation_db',fadb
%                    Filter attenuation for the two gammatone filter used for the fine
%                    structure and envelope of the modulation filter bank. The
%                    default value is 10.
%                     
%     'fine_filter_finesse',fff
%                    Filter finesse (determines the bandwidth with fc/finesse)
%                    for the fine structure gammatone filter. The defulat value
%                    is 3.
%
%     'mod_center_frequency_hz',mcf_hz
%                    Center frequency of the gammatone envelope filter. The
%                    default value is 135.
%
%     'mod_filter_finesse',mff
%                    Filter finesse (determines the bandwidth with fc/finesse)
%                    for the envelope gammatone filter. The defulat value is 8.
% 
%     'level_filter_cutoff_hz',lfc_hz
%                    Cutoff frequency off the low pass filter used for ILD
%                    calculation. The default value is 30.
%
%     'level_filter_order',lforder
%                    Order of low pass filter for the ILD calculation. The
%                    default value is 2.
%
%   See also: dietz2011, dietz2011_interauralfunctions
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592--605, 2011. [1]http ]
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S016763931000097X
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/dietz2011_filterbank.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Tobias Peters (2002)
%   #Author: Mathias Dietz (2006-2009)
%   #Author: Martin Klein-Hennig (2011)
%   #Author: Martin Klein-Hennig 
%   #Author: Hagen Wierstorf (2013): for AMT 

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
 
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig)
  error('%s: insig has to be a numeric signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

if ~isnumeric(fc)
  error('%s: fc has to be a numeric signal!',upper(mfilename));
end

definput.import = {'dietz2011_filterbank'};
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Model processing starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Split the input signal:
% - below 1.4 kHz the fine structure filter is applied
% - above 1.4 kHz the modulation filter is applied
fc_fine = fc(fc<=1400);
fc_env = fc(fc>1400);
insig_fine = insig(:,fc<=1400,:);
insig_env = insig(:,fc>1400,:);

% --- fine structur filter ---
outsig_fine = zeros(size(insig_fine));
% gammatone filter centered at the center frequency for every frequency channel
for ii=1:length(fc_fine)
  gammatone_filter = hohmann2002_filter (fs, fc_fine(ii), ...
    fc_fine(ii)/kv.fine_filter_finesse, ...
    kv.filter_attenuation_db, kv.filter_order);
  outsig_fine(:,ii,1) = hohmann2002_process(gammatone_filter, squeeze(insig_fine(:,ii,1)));
  outsig_fine(:,ii,2) = hohmann2002_process(gammatone_filter, squeeze(insig_fine(:,ii,2)));
end

% --- modulation/envelope filter ---
outsig_env = zeros(size(insig_env));
% gammatone filter centered at a fixed frequency for every frequency channel
gammatone_filter = hohmann2002_filter (fs, kv.mod_center_frequency_hz, ...
  kv.mod_center_frequency_hz/kv.mod_filter_finesse, ...
  kv.filter_attenuation_db, kv.filter_order);
for ii=1:length(fc_env)
  outsig_env(:,ii,1) = hohmann2002_process(gammatone_filter, squeeze(insig_env(:,ii,1)));
  outsig_env(:,ii,2) = hohmann2002_process(gammatone_filter, squeeze(insig_env(:,ii,2)));
end

% --- ILD filter ---
% low pass filter with a fixed cutoff frequency for every frequency channel
[b,a] = butter(kv.level_filter_order,kv.level_filter_cutoff_hz/(fs/2),'low');
outsig_ild = filter(b,a,insig);

% vim: set sw=2 ts=2 expandtab textwidth=80: 


