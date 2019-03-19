function [benefit, weighted_SNR, weighted_bmld] = jelfs2011(target,interferer,varargin)
%JELFS2011  Predicted binaural advantage for speech in reverberant conditions
%   Usage:  [benefit weighted_SNR weighted_bmld] = jelfs2011(target,interferer,fs)
%  
%   Input parameters:
%     target        : Binaural target impulse respone (or stimulus)
%     interfererer  : Binaural interferer impulse response (or stimulus)
%                     Multiple interfering impulse responses MUST be
%                     concatenated, not added.
%
%   Output parameters:
%     benefit       : spatial release from masking (SRM)in dB
%     weighted_SNR  : component of SRM due to better-ear listening (dB)
%     weighted_bmld : component of SRM due to binaural unmasking (dB)
%    
%   JELFS2011(target,interferer,fs) computes the increase in speech
%   intelligibility of the target when the target and interferer are
%   spatially separated. They are preferably represented by their impulse
%   responses, but can be represented by noise recordings of equivalent
%   spectral shape emitted from the same source locations (using the same
%   noise duration for target and interferer). The impulse responses are
%   assumed to be sampled at a sampling frequency of fs Hz.  If the
%   modelled sources differ in spectral shape, this can be simulated by
%   pre-filtering the impulse responses.
%
%   [benefit, weighted_SNR, weighted_bmld]=JELFS2011(...) additionaly
%   returns the benefit from the SII weighted SNR and the SII weighted BMLD.
%
%   If target or interferer are cell-arrays, the HRTF data will be loaded. The first
%   argument in the cell-array is the azimuth angle, and the second
%   parameter is the database type. The elevation is set to zero.
%   function. 
%
%   Example:
%   --------
%
%   The following code will load HRIRs from the 'kemar' database and
%   compute the binaural speech intelligibility advantage for a target
%   at 0 degrees and interferers at 300 and 90 degrees:
%
%     jelfs2011({0,'kemar'},{[330 90],'kemar'})
%
%   See also: culling2005_bmld, exp_jelfs2011
% 
%   References:
%     J. Culling, S. Jelfs, and M. Lavandier. Mapping Speech Intelligibility
%     in Noisy Rooms. In Proceedings of the 128th convention of the Audio
%     Engineering Society, Convention paper 8050, 2010.
%     
%     S. Jelfs, J. Culling, and M. Lavandier. Revision and validation of a
%     binaural model for speech intelligibility in noise. Hearing Research,
%     2011.
%     
%     M. Lavandier, S. Jelfs, J. Culling, A. Watkins, A. Raimond, and
%     S. Makin. Binaural prediction of speech intelligibility in reverberant
%     rooms with multiple noise sources. J. Acoust. Soc. Am., 131(1):218-231,
%     2012.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/jelfs2011.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  definput.flags.ears={'both','left','right'};
  definput.keyvals.fs=[];
  definput.keyvals.pad=1024;
  [flags,kv,fs]=ltfatarghelper({'fs'},definput,varargin);
  
  % If target or interferer are cell arrays, load HRTFs.
  if iscell(target)
    X=SOFAload(fullfile(amt_basepath,'hrtf',mfilename,[target{2} '.sofa']));
    idx=find(X.SourcePosition(:,1)==target{1} & X.SourcePosition(:,2)==0);
    target=squeeze(X.Data.IR(idx,:,:))';
    target=postpad(target,size(target,1)+kv.pad);
    fs=X.Data.SamplingRate;
  end;
  
  if iscell(interferer)
    azims=numel(interferer{1});
    X=SOFAload(fullfile(amt_basepath,'hrtf',mfilename,[interferer{2} '.sofa']));
    for ii=1:azims
      idx(ii)=find(X.SourcePosition(:,1)==mod(interferer{1}(ii),360) & X.SourcePosition(:,2)==0);
    end
    interferer=shiftdim(X.Data.IR(idx,:,:),2);
    interferer=postpad(interferer,size(interferer,1)+kv.pad);
    fs2=X.Data.SamplingRate;
    
    if fs2~=fs
      error('%s: Mis-match between target and interferer sampling rate.',upper(mfilename));
    end;
    % Old code compatibility
    if ndims(interferer)==3
      s=size(interferer);
      interferer=reshape(interferer,s(1)*s(2),s(3));
      interferer=interferer/sqrt(azims);
    end;
  end;  
  
  if isempty(fs)
    error('%s: You must specify the sampling rate, fs.',upper(mfilename));
  end;
  
  % Make sure that there is at least 1 erb per channel, and get
  % the gammatone filters.
  nchannels=ceil(freqtoerb(fs/2));
  fc=erbspace(0,fs/2,nchannels);
  [b,a] = gammatone(fc,fs,'complex');
  
  effective_SNR = zeros(nchannels,1);
  bmld_prediction = zeros(nchannels,1);
  
  targ_f = 2*real(ufilterbankz(b,a,target));
  int_f  = 2*real(ufilterbankz(b,a,interferer));
  
  for n = 1:nchannels
    % Calculate the effect of BMLD
    if flags.do_both
      % cross-correlate left and right signal in channel n for both the
      % target and the inteferer
      [phase_t, coher_t] = do_xcorr(targ_f(:,n,1),targ_f(:,n,2),fs,fc(n)); 
      [phase_i, coher_i] = do_xcorr( int_f(:,n,1), int_f(:,n,2),fs,fc(n)); 
      
      bmld_prediction(n) = culling2005_bmld(coher_i,phase_t,phase_i,fc(n));
    end
    
    % Calculate the effect of better-ear SNR
    left_SNR  = sum(targ_f(:,n,1).^2) / sum(int_f(:,n,1).^2);
    right_SNR = sum(targ_f(:,n,2).^2) / sum(int_f(:,n,2).^2);
    
    if flags.do_both
      SNR = max(left_SNR,right_SNR);
    end;
    
    if flags.do_left
      SNR = left_SNR;
    end;
    
    if flags.do_right
      SNR = right_SNR;
    end
    
    % combination
    effective_SNR(n) = 10*log10(SNR);
  end
  
  % Calculate the SII weighting
  weightings = siiweightings(fc);
  
  if flags.do_both
    weighted_bmld = sum(bmld_prediction.*weightings);
  else
    weighted_bmld = 0;
  end
  
  weighted_SNR = sum(effective_SNR.*weightings);
  
  benefit = weighted_SNR + weighted_bmld;
end

% Helper function to do the cross-correlation, and extract the delay of
% the peak (output parameter 'phase' and the coherence at the peak).
function [phase, coherence] = do_xcorr(left, right, fs, fc)
  
  % Use the LTFAT correlation function to avoid depending on xcorr, which
  % is not in core of Matlab or Octave.
  iacc = pxcorr(squeeze(left),squeeze(right),'normalize');

  % Find the position of the largest correlation coefficient.
  [coherence, delay_samp] = max(iacc);
  if delay_samp > length(iacc)/2
     delay_samp=delay_samp-length(iacc);
  end
  phase = fc*2*pi*delay_samp/fs;
end
  
