function [Cohc,Cihc,OHC_Loss]=carney2015_fitaudiogram(FREQUENCIES,dBLoss,species,Dsd_OHC_Loss)
%CARNEY2015_FITAUDIOGRAM Cohc and Cihc values that produce a desired threshold shift 
%
%   Usage:
%     [Cohc,Cihc,OHC_Loss]=carney2015_fitaudiogram(FREQUENCIES,dBLoss,species,Dsd_OHC_Loss)
%
%   Input parameters:
%     FREQUENCIES : vector containing audiogram frequencies
%     dBLoss      : loss [dB] per frequency in FREQUENCIES
%     species     : model species "1" for cat, "2" for human BM tuning from
%                   Shera et al. (PNAS 2002), or "3" for human BM tuning from Glasberg &
%                   Moore (Hear. Res. 1990)
%     Dsd_OHC_Loss: optional array giving the desired threshold shift in
%                   dB that is caused by the OHC impairment alone (for each frequency in
%                   FREQUENCIES). If this array is not given, then the default desired
%                   threshold shift due to OHC impairment is 2/3 of the entire threshold
%                   shift at each frequency.  This default is consistent with the
%                   effects of acoustic trauma in cats (see Bruce et al., JASA 2003, and
%                   Zilany and Bruce, JASA 2007) and estimated OHC impairment in humans
%                   (see Plack et al., JASA 2004).
%
%   Output parameters:
%     The output variables are arrays with values corresponding to each
%     frequency in the input array FREQUENCIES.
%
%     Cohc     : is the outer hair cell (OHC) impairment factor; a value of 1
%                corresponds to normal OHC function and a value of 0 corresponds to
%                total impairment of the OHCs.
%     Cihc     : is the inner hair cell (IHC) impairment factor; a value of 1
%                corresponds to normal IHC function and a value of 0 corresponds to
%                total impairment of the IHC.
%     OHC_Loss : is the threshold shift in dB that is attributed to OHC
%                impairment (the remainder is produced by IHC impairment). 
%
%
%   Cohc and Cihc values that produce a desired threshold shift 
%   for the cat & human auditory-periphery model of Zilany et
%   al. (J. Acoust. Soc. Am. 2009, 2014) and Bruce, Erfani & Zilany (Hear.Res.).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/carney2015_fitaudiogram.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Authors: University of Rochester (UR EAR) team
%   #Author: M. S. A. Zilany 
%   #Author: I. C. Bruce (ibruce@ieee.org), 2013
%   #Authors: Clara Hollomey (2020): integration in the AMT
%   #Authors: Piotr Majdak (2021): integration for the AMT 1.0
%   #Authors: Alejandro Osses (2021): extensions for the AMT 1.1

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% CHANGES to original version: load('THRESHOLD_XXX) => amt_load(XXX) (C.H., Apr 2021)
% finished function with an 'end' (C.H., Apr 2021)

% Keep the contents of the mat-file in memory unless species changes.
persistent last_species file
if ~isequal(species,last_species)
	switch species
		case 1
% 			disp('Analyzing audiogram for cat AN model')
			file = amt_load('bruce2018', 'THRESHOLD_ALL_CAT.mat');
		case 2
% 			disp('Analyzing audiogram for human AN model - BM tuning from Shera et al. (2002)')
            file = amt_load('bruce2018', 'THRESHOLD_ALL_HM_Shera.mat');
			%file = load('THRESHOLD_ALL_HM_Shera','CF','CIHC','COHC','THR');
		case 3
% 			disp('Analyzing audiogram for human AN model - BM tuning from Glasberg & Moore (1990)')
			%file = load('THRESHOLD_ALL_HM_GM','CF','CIHC','COHC','THR');
            file = amt_load('bruce2018', 'THRESHOLD_ALL_HM_GM.mat');
		otherwise
			error('Species #%d not known.',species)
	end
	last_species = species;
end


% Variables are
% CF: 125 Hz to 10 kHz                   [1*37]
% CIHC: varies from 1.0 to 0.0001        [1*55]
% COHC: varies from 1.0 to 0             [1*56]
% THR : absolute thresholds              [37*55*56]

% for k = 1:length(file.THR(:,1,1))
% 	dBShift(k,:,:) = file.THR(k,:,:) - file.THR(k,1,1);
% end
try
	dBShift = file.THR - file.THR(:,1,1);
catch
	dBShift = bsxfun(@minus,file.THR,file.THR(:,1,1));
end


if nargin < 4
	Dsd_OHC_Loss = 2/3*dBLoss;
end

num_freq = length(FREQUENCIES);
Cohc = zeros(1,num_freq);
OHC_Loss = zeros(1,num_freq);
Loss_IHC = zeros(1,num_freq);
Cihc = zeros(1,num_freq);
for m = 1:length(FREQUENCIES)
	[~,N] = min(abs(file.CF - FREQUENCIES(m)));
	n = N(1);
	
	if Dsd_OHC_Loss(m) > dBShift(n,1,end)
		Cohc(m) = 0;
	else
		[~,idx] = sort(abs(squeeze(dBShift(n,1,:)) - Dsd_OHC_Loss(m)));
		Cohc(m) = file.COHC(idx(1));
	end
	OHC_Loss(m) = interp1(file.COHC,squeeze(dBShift(n,1,:)),Cohc(m),'nearest');
	[~,ind] = sort(abs(file.COHC - Cohc(m)));
	
	Loss_IHC(m) = dBLoss(m) - OHC_Loss(m);
	
	if dBLoss(m) > dBShift(n,end,ind(1))
		Cihc(m) = 0;
	else
		[~,indx] = sort(abs(squeeze(dBShift(n,:,ind(1))) - dBLoss(m)));
		Cihc(m) = file.CIHC(indx(1));
	end
	
end
end


