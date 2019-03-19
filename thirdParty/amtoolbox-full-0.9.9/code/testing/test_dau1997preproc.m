insig=greasy;
fs=16000;

[outsig_ref, fc_ref, mfc_ref] = ref_dau1997_preproc(insig, fs);

% The reference must use a 'basef', so also do this here.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_dau1997preproc.php

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
[outsig, fc, mfc] = dau1997_preproc(insig, fs, 'basef', 1000);

norm(fc-fc_ref)
norm(mfc-mfc_ref)

for nfc=1:length(fc)
  for nmfc=1:size(outsig{nfc},2)
    res=norm(outsig{nfc}(:,nmfc)-outsig_ref(:,nfc,nmfc))/norm(outsig{nfc}(:,nmfc));
    
    amt_disp(sprintf('Freq: %f MFreq: %f res: %f',fc(nfc),mfc(nmfc),res));
  end;  
end;



