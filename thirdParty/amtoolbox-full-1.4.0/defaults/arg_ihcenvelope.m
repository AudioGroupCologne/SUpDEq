function definput=arg_ihcenvelope(definput)
%ARG_IHCENVELOPE
%   #License: GPL
%   #Author: Peter Soendergaard (2011): Initial version
%   #Author: Alejandro Osses (2020): Extensions
%   #Author: Piotr Majdak (2021): Adapted to AMT 1.0
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_ihcenvelope.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
 
  definput.flags.ihc = {'ihc','no_ihc'}; 
  definput.flags.ihc_type={'ihc_undefined','ihc_bernstein1999','ihc_breebaart2001','ihc_dau1996','hilbert', ...
                    'ihc_lindemann1986','ihc_meddis1990','ihc_king2019','ihc_relanoiborra2019'};

  definput.keyvals.ihc_minlvl=[];
  definput.keyvals.ihc_filter_order=1;
  definput.keyvals.ihc_scal_constant=[];
  
  definput.groups.ihc_breebaart2001={'ihc_filter_order',5};
  definput.groups.ihc_dau1996={'ihc_filter_order',1};



