function decision = breebaart2001_centralproc(EI_map,monol,monor,bimonostring,monofactor)
%BREEBAART2001_CENTRALPROC   Binaural and monoaural signal detector from Breebaart et. al. 2001
%   Usage: decision = breebaart2001_centralproc(EI_map,monol,monor,bimonostring);
%          decision = breebaart2001_centralproc(EI_map,monol,monor,bimonostring,monofactor);
%
%   Input parameters:
%        EI_map         : binaural representation of the signal (from breebaart2001_eicell)
%        monol          : internal representation of the left ear signal
%        monor          : internal representation of the right ear signal
%        bimonostring   : defines which representations are used for the decision
%        monofactor     : sets the parameter for monaural detection
%                         (optional)
%  
%   Output parameters:
%        decision       : index of the estimated target
%   
%   BREEBAART2001_CENTRALPROC(EI_map,monol,monor,bimonostring,monofactor) serves
%   as an artificial observer for signal detection purposes. The
%   central processor develops a template which consists of the average
%   internal representation of all masker-alone intervals and its variance.
%   Furthermore the average signal interval and the average distance 
%   between masker template and signal template is computed. Then the
%   weighted distance U between template and actual stimulus is calculated. 
%   The higher U, the greater the likelihood that a signal is present. 
%   Thus, in a 3-IFC procedure the model will choose the interval with the 
%   highest value of U. After each trial the model recieves feedback and
%   stores the avaluated intervals in the templates. 
%   
%   The parameter bimonostring must contain characters indicating which 
%   channels are used for the decision:
%
%     'l'     use left mono channel (from monol*)
%
%     'r'     use right mono channel (from monor*)
%
%     'b'     use binaural channel (from EI_map ) only if the binaural representation
%             yield a non-zero decision distance
%
%     'B'     use binaural channel in any case
%
%
%   See also: exp_breebaart2001 demo_breebaart2001 breebaart2001 breebaart2001_eicell
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074--1088, August 2001.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/breebaart2001_centralproc.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Peter L. Soendergaard (2011)
%   #Author: Martina Kreuzbichler (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if nargin == 4
    monofactor = 0.0003;
end

persistent maskernumber;
persistent signalnumber;

intnum = length(EI_map);
tempsize = size(EI_map{1});
intnoisevar = 1;
binauralset = 0;
monolset = 0;
monorset = 0;

if isempty(maskernumber)
    maskernumber = intnum-1;
    signalnumber = 1;
else
    maskernumber = maskernumber+intnum-1;
    signalnumber = signalnumber + 1;
end
 
% binaural wanted
if strfind(lower(bimonostring),'b')
    binauralset = 1;
    persistent template; 
    persistent templatesq;
    persistent signaltemplate;
    U_b = zeros(1,intnum);
    
    if maskernumber == 2
        template = zeros(tempsize);
        templatesq = zeros(tempsize);
        signaltemplate = zeros(tempsize);
    end
    
    noisevar = templatesq - template.^2;
    meandiff = signaltemplate - template;
    weight = meandiff./(noisevar + intnoisevar);
    Nuvar = intnoisevar*sum(sum(weight.^2));        
end

% left mono channel wanted
if strfind(lower(bimonostring),'l')
    monolset = 1;
    persistent template_ml;
    persistent templatesq_ml;
    persistent signaltemplate_ml;
    monol = cellfun(@(x) x*monofactor,monol,'un',0);
    U_ml = zeros(1,intnum);
    
    if maskernumber == 2
        template_ml = zeros(tempsize);
        templatesq_ml = zeros(tempsize);
        signaltemplate_ml = zeros(tempsize);
    end
    
    noisevar_ml = templatesq_ml - template_ml.^2;
    meandiff_ml = signaltemplate_ml - template_ml;
    weight_ml = meandiff_ml./(noisevar_ml + intnoisevar);
    Nuvar_ml = intnoisevar*sum(sum(weight_ml.^2));
end

% right mono channel wanted
if strfind(lower(bimonostring),'r')
    monorset = 1;
    persistent template_mr;
    persistent templatesq_mr;
    persistent signaltemplate_mr;
    monor = cellfun(@(x) x*monofactor,monor,'un',0);
    U_mr = zeros(1,intnum);
    
    if maskernumber == 2
        template_mr= zeros(tempsize);
        templatesq_mr = zeros(tempsize);
        signaltemplate_mr = zeros(tempsize);
    end
    
    noisevar_mr = templatesq_mr - template_mr.^2;
    meandiff_mr = signaltemplate_mr - template_mr;
    weight_mr = meandiff_mr./(noisevar_mr + intnoisevar);
    Nuvar_mr = intnoisevar*sum(sum(weight_mr.^2));
end


for intcount = 1:intnum
    noise = randn;
    if binauralset == 1
        U_b(intcount) = sum(sum(weight.*(EI_map{intcount}-template)))...
            + noise*sqrt(Nuvar);
    end
    if monolset == 1
        U_ml(intcount) = sum(sum(weight_ml.*(monol{intcount}-template_ml)))...
            + noise*sqrt(Nuvar_ml);
    end
    if monorset ==1
        U_mr(intcount) = sum(sum(weight_mr.*(monor{intcount}-template_mr)))...
            + noise*sqrt(Nuvar_mr);
    end
end

% take binaural, mono left and mono right
if binauralset && monolset && monorset
    % binaural decision is empty and no big B
    if signalnumber > 1 && any(U_b) == 0 && ~isempty(strfind(bimonostring,'b')) 
        U = mean([U_ml;U_mr]);
    else
        U = mean([U_b;U_ml;U_mr]);
    end
% take binaural and mono left
elseif binauralset && monolset
    U = mean([U_b;U_ml]);
% take binaural and mono right
elseif binauralset && monorset
    U = mean([U_b;U_mr]);
% take mono left and mono right
elseif monorset && monolset
    U = mean([U_ml;U_mr]);
% take binaural
elseif binauralset
    U = U_b;
% take mono left
elseif monolset
    U = U_ml;
% take mono right
elseif monorset
    U = U_mr;
end
    
% If centralproc is called the first time, all the templates are empty
% In this case the results of U will be zeros. Therefore the response
% of max(U) will always be 1 = first occurence of 0.
[~,response] = max(U);


if binauralset
    % update of templates
    signaltemplate = ((signaltemplate*(signalnumber-1))+...
        (EI_map{1}))./signalnumber;
    
    % binaural test
    [~,response_b] = max(U_b);
    if response_b ~= response && any(U_b) == 1
        amt_disp(sprintf(['Signal #%i: Monaural-based decision dominated the final decision: '...
            ' final decision = %i, binaural-based decision = %i \n'],...
            signalnumber,response,response_b));
    end
    adtemplate = zeros(tempsize);
    adtemplatesq = zeros(tempsize);
    
    for updatecounter = 2:intnum
        adtemplate = adtemplate + EI_map{updatecounter};
        adtemplatesq = adtemplatesq + EI_map{updatecounter}.^2;
    end
    
    template = ((template*(maskernumber-(intnum-1)))...
        +(adtemplate))./maskernumber;
    templatesq = ((templatesq*(maskernumber-(intnum-1)))+...
        (adtemplatesq).^2)./maskernumber; 
end

if monolset
   % update of templates 
    signaltemplate_ml = ((signaltemplate_ml*(signalnumber-1))+...
        (monol{1}))./signalnumber;
    adtemplate_ml = zeros(tempsize);
    adtemplatesq_ml = zeros(tempsize);
    
    for updatecounter = 2:intnum
        adtemplate_ml = adtemplate_ml + monol{updatecounter};
        adtemplatesq_ml = adtemplatesq_ml + monol{updatecounter}.^2;
    end
    
    template_ml = ((template_ml*(maskernumber-(intnum-1)))...
        +(adtemplate_ml))./maskernumber;
    templatesq_ml = ((templatesq_ml*(maskernumber-(intnum-1)))...
        +(adtemplatesq_ml).^2)./maskernumber; 
end

if monorset
   % update of templates
   signaltemplate_mr = ((signaltemplate_mr*(signalnumber-1))+...
       (monor{1}))./signalnumber;
   adtemplate_mr = zeros(tempsize);
   adtemplatesq_mr = zeros(tempsize);
   
   for updatecounter = 2:intnum
       adtemplate_mr = adtemplate_mr + monor{updatecounter};
       adtemplatesq_mr = adtemplatesq_mr + monor{updatecounter}.^2;
   end
   
   template_mr = ((template_mr*(maskernumber-(intnum-1)))+...
       (adtemplate_mr))./maskernumber;
   templatesq_mr = ((templatesq_mr*(maskernumber-(intnum-1)))+...
       (adtemplatesq_mr).^2)./maskernumber; 
end

decision = response;


