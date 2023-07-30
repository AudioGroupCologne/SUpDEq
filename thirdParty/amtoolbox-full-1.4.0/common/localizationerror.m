function [varargout] = localizationerror(m,varargin)
%LOCALIZATIONERROR Compute psychoacoustic performance parameters for sound localization experiments
%   Usage: [accL, precL, accP, precP, querr] = localizationerror(m)
%          [res, meta, par] = localizationerror(m,errorflag)
%
%   Input parameters:
%
%     m       : item list from a localization experiment. 
%
%   LOCALIZATIONERROR(m,errorflag) returns psychoacoustic performance 
%   parameters for experimental response patterns. 
%   m is a matrix. In each row of m, the information about the target 
%   and response in the columns must be provided.
%
%   [res, meta, par] = LOCALIZATIONERROR(m,errorflag) calculates error 
%   metric res given by the string errorflag. meta contains additional 
%   metadata, par contains additional calculation results.
%
%   The columns of 'm' comprise: 
% 
%   - 1:4 azi_target,ele_target,azi_response,ele_response
%
%   - 5:8 lat_target,pol_target,lat_response,pol_response
%
%   - 9   F/B-C resolved pol_response
%
%
%   The errorflag may be one of:
%
%     'accL'               accuracy (average) in the lateral dimension
%
%     'precL'              precision (std. dev.) in the lateral dimension
%
%     'precLcentral'       'precL' derived from responses close to horizontal meridian (+-30 deg elevation)
%
%     'accP'               accuracy (average) in the polar dimension
%
%     'precP'              precision (circular std. dev.) in the polar dimension 
%
%     'querr'              percentage of weighted polar errors > 45deg:[querr number_of_corrects total_number] = CalcAccPrec...
%
%     'accE'               accuracy in the elevation
%
%     'absaccE'            absolute accuracy in the elevation
%
%     'absaccL'            absolute lateral accuracy 
%
%     'accabsL'            mean absolute lateral error
%
%     'accabsP'            mean absolute polar error
%
%     'accPnoquerr'        accP with quadrant errors removed
%
%     'precE'              precision in the elevation
%
%     'querr90'            precentage of weighted polar errors > 90deg
%
%     'precPmedian'        unweighted polar precision error considering targets around the median plane only (+/-30deg lateral) 
%
%     'precPmedianlocal'   basing on precPmedian, quadrant errors (polar error >90 deg) are excluded. Similar to RMS local polar error (Middlebrooks, 1999)
%
%     'precPnoquerr'       precP with quadrant errors removed
%
%     'rmsL'               lateral RMS error according to Middlebrooks (1999)
%
%     'rmsPmedianlocal'    unweighted polar RMS considering targets around
%                          the median plane only (+/-30deg lateral) and excluding quadrant
%                          errors (polar error > 90deg). Identical to RMS local polar error in Middlebrooks (1999).
%
%     'rmsPmedian'         unweighted polar RMS considering targets around the median plane 
%                          only (+/-30deg lateral). Includes quadrant errors.
%
%     'querrMiddlebrooks'  quadrant errors as in Middlebrooks (1999). Use
%                          central positions within +/-30deg, querr is when polar
%                          error is >90deg [querr number_of_corrects total_number] = CalcAccPrec...
%
%     'corrcoefL'          lateral correlation coeffficient. Output parameter: [cc, p]
%
%     'corrcoefP'          polar correlation coeffficient, Only data within lateral +/-30deg are considered. Output parameter: [corr_coeff, p_of_significant_corr]
%
%     'SCC'                spherical/spatial correlation coefficient (Carlile et al., 1997)
%
%     'gainLstats'         performs linear regression to obtain lateral gain. See help of regress for a detailed description of the structure fields.
%
%     'gainL'              lateral gain only
%
%     'pVeridicalL'        Proportion of quasi-verdical responses in lateral
%                          dimension
%
%     'precLregress'       lateral scatter around linear regression line (only quasi-veridical responses included).
%
%     'sirpMacpherson2000' performs an ad-hoc selective, iterative
%                          regression procedure (SIRP) as proposed in
%                          Macpherson & Middlebrooks (2000 JASA) in order 
%                          to exclude outliers and reversals 
%                          and isolate the main concentration of responses
%                          in the computation of the linear fits. 
%                          As the original procedure does not necessarily
%                          converge, the implementation here runs up to 100
%                          iterations and selects the iteration that includes 
%                          the most responses or excludes the least outliers,
%                          Outlier distance criterion is 40deg, output parameters 
%                          [f,r] correspond to regression results for frontal and 
%                          rear hemisphere respectively, see help of regress for 
%                          a detailed description of the structure fields.
%
%     'gainPfront'         frontal polar gain only
%
%     'gainPrear'          rear polar gain only
%
%     'gainP'              polar gain averaged between front and back
%
%     'slopePfront'        slope in degrees of regression line (frontal only)
%
%     'slopePrear'         slope in degrees of regression line (rear only)
%
%     'slopeP'             slope in degrees of regression line (front and back)
%
%     'pVeridicalPfront'   proportion of quasi-verdical polar responses in the front
%
%     'pVeridicalPrear'    proportion of quasi-verdical polar responses in the back
%
%     'pVeridicalP'        proportion of quasi-verdical polar responses in total
%
%     'precPregressFront'  polar scatter around linear regression line for the front (only quasi-veridical responses included)
%
%     'precPregressRear'   same as precPregressFront but for the back
%
%     'precPregress'       average between precPregressFront and precPregressRear*
%
%     'perMacpherson2003'  polar error rate used in Macpherson & Middlebrooks
%                          (2003). They measured the deviation of responses from the
%                          linear predictors obtained by an ad-hoc SIRP. Polar errors  
%                          are defined by showing a deviation of >45deg 
%                          (i.e., not being quasi-veridical) with respect 
%                          to the linear flat stimulus prediction. Note that for this
%                          analysis the results from sirpMacpherson2000 are required 
%                          and handled as LOCALIZATIONERROR(m,f,r,'perMacpherson2003')
%
%
%   If no errorflag is provided, the function returns: accL, precL, precP, and querr
%
% 
%   References:
%     J. C. Middlebrooks. Virtual localization improved by scaling
%     nonindividualized external-ear transfer functions in frequency. J.
%     Acoust. Soc. Am., 106:1493--1510, 1999.
%     
%     E. A. Macpherson and J. C. Middlebrooks. Vertical-plane sound
%     localization probed with ripple-spectrum noise. J. Acoust. Soc. Am.,
%     114:430--445, 2003.
%     
%     S. Carlile, P. Leong, and S. Hyams. The nature and distribution of
%     errors in sound localization by human listeners. Hearing research,
%     114(1):179--196, 1997.
%     
%
%   See also: baumgartner2014
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/localizationerror.php


%   #Author : Piotr Majdak (2021) 
%   #Author : Robert Baumgartner (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Check input options 

definput.import={'localizationerror'};

[flags,kv]=ltfatarghelper({'f','r'},definput,varargin);


if size(m,2)==1,   m=m'; end

    
if isempty(flags.errorflag)
  nargout=6;
  thr=45;

  accL=mean(m(:,7)-m(:,5));
  varargout{1}=accL;
  precL=sqrt(sum((m(:,7)-m(:,5)-accL).^2)/size(m,1));  
  varargout{2}=precL;
  w=0.5*cos(deg2rad(m(:,7))*2)+0.5;  
  varargout{3}=mynpi2pi(circmean(m(:,8)-m(:,6),w));
  [x,mag]=circmean(m(:,8)-m(:,6));
  precP=rad2deg(acos(mag)*2)/2;        
  varargout{4}=precP;
  idx=find((abs(mynpi2pi(m(:,8)-m(:,6))).*w)<=thr);  
  querr=100-size(idx,1)/size(m,1)*100;
  varargout{5}=querr;
else
%   if iscell(errorflag), errorflag=errorflag{1}; end
%   if ~ischar(errorflag)
%     error('ERR must be a string');
%   end
%   if size(errorflag,1)~=1, errorflag=errorflag'; end
  nargout=3;
  metadata=[];
  par=[];
  switch flags.errorflag
    case 'accL'
      accL=mean(m(:,7)-m(:,5));  
      varargout{1}=accL;
      meta.ylabel='Lateral bias (deg)';
    case 'absaccL'
      accL=mean(m(:,7)-m(:,5));  
      varargout{1}=abs(accL);
      meta.ylabel='Lateral absolute bias (deg)';
    case 'accabsL'
      accL=mean(abs(m(:,7)-m(:,5)));  
      varargout{1}=accL;
      meta.ylabel='Lateral absolute error (deg)';
    case 'corrcoefL'
      [r,p] = corrcoef([m(:,5) m(:,7)]);
      varargout{1}=r(1,2);
      par.par1=p(1,2);
      meta.ylabel='Lateral correlation coefficient';
      meta.par1 = 'p for testing the hypothesis of no correlation';
    case 'corrcoefP'
      idx=find(abs(m(:,7))<=30);
      [r,p] = corrcoef([m(idx,6) m(idx,8)]);
      varargout{1}=r(1,2);
      par.par1=p(1,2);
      meta.ylabel='Polar correlation coefficient';
      meta.par1 = 'p for testing the hypothesis of no correlation';
    case 'precL'
      accL=mean(m(:,7)-m(:,5));  
      precL=sqrt(sum((m(:,7)-m(:,5)-accL).^2)/size(m,1));  
      varargout{1}=precL;
      meta.ylabel='Lateral precision error (deg)';
    case 'precLcentral'
      idx=find(abs(m(:,4))<=30); % responses close to horizontal meridian
      accL=mean(m(idx,7)-m(idx,5));  
      precL=sqrt(sum((m(idx,7)-m(idx,5)-accL).^2)/length(idx));  
      varargout{1}=precL;
      meta.ylabel='Central lateral precision error (deg)';
    case 'rmsL'
      idx=find(abs(m(:,7))<=60);
      m=m(idx,:);        
      rmsL=sqrt(sum((mynpi2pi(m(:,7)-m(:,5))).^2)/size(m,1));
      varargout{1}=rmsL;
      meta.ylabel='Lateral RMS error (deg)';
    case 'accP'
      varargout{1}=mynpi2pi(circmean(m(:,8)-m(:,6)));
      meta.ylabel='Polar bias (deg)';
    case 'accabsP'
      varargout{1}=mynpi2pi(circmean(abs(m(:,8)-m(:,6))));
      meta.ylabel='Polar absolute bias (deg)';
    case 'accPnoquerr'
      w=0.5*cos(deg2rad(m(:,7))*2)+0.5;
      idx=find(w>0.5);
      m=m(idx,:); w=w(idx,:);
      idx=find((abs(mynpi2pi(m(:,8)-m(:,6))).*w)<=45);                  
      varargout{1}=mynpi2pi(circmean(m(idx,8)-m(idx,6)));     
      par.par1=length(idx);
      meta.ylabel='Local polar bias (deg)';
      meta.par1='Number of considered responses';
    case 'precP'
      w=0.5*cos(deg2rad(m(:,7))*2)+0.5;
      precP=circstd(m(:,8)-m(:,6),w);
      varargout{1}=precP;
      meta.ylabel='Polar precision (deg)';
    case 'precPnoquerr'
%         w=0.5*cos(deg2rad(m(:,7))*2)+0.5;
%         idx=find(w>0.5);
%         thr=45; % threshold for the choice: quadrant error or not        
%         m=m(idx,:); w=w(idx,:);
      range=360; % polar range of targets
      thr=45; % threshold for the choice: quadrant error or not
      latmax=0.5*180*(acos(4*thr/range-1))/pi;
      idx=find(abs(m(:,7))<latmax); % targets outside not considered
      m=m(idx,:);
      w=0.5*cos((m(:,7))*2*pi/180)+0.5;
      idx=find((abs(mynpi2pi(m(:,8)-m(:,6))).*w)<=thr);                  
      m=m(idx,:);
      w=0.5*cos((m(:,7))*pi/180*2)+0.5;
      precP=circstd(m(:,8)-m(:,6),w);
      varargout{1}=real(precP);        
      meta.ylabel='Local polar precision error (deg)';
    case 'querr'
      range=360; % polar range of targets
      thr=45; % threshold for the choice: quadrant error or not
      latmax=0.5*180*(acos(4*thr/range-1))/pi;
      idx=find(abs(m(:,7))<latmax); % targets outside not considered
      m=m(idx,:);
      w=0.5*cos(pi*(m(:,7))/180*2)+0.5; % weigthing of the errors
      corrects=size(find((abs(mynpi2pi(m(:,8)-m(:,6))).*w)<=thr),1);
      totals=size(m,1);
      chancerate=2*thr/range./w;
      wrongrate=100-mean(chancerate)*100;
      querr=100-corrects/totals*100;
      varargout{1}=querr;
      meta.ylabel='Weighted quadrant errors (%)';
      par.par1=corrects;
      meta.par1='Number of correct responses';
      par.par2=totals;
      meta.par2='Number of responses within the lateral range'; 
      par.par3=wrongrate; % valid only for targets between 0 and 180°
      meta.par3='Wrong rate, valid only for targets between 0 and 180°';
      
    case 'querrMiddlebrooks' % as in Middlebrooks (1999); comparable to COC from Best et al. (2005)
      idx=abs(m(:,7))<=30;
      m=m(idx,:);
      idx=find((abs(mynpi2pi(m(:,8)-m(:,6))))<=90);  
      querr=100-size(idx,1)/size(m,1)*100;
      % chance performance (FIXME: not working with exp_baumgartner2014!)
%       tang = m(:,6); 
%       rang = min(tang):max(tang); % assume listener knows target range
%       p = ones(length(rang),length(tang)); % equally distributed response probabilities
%       chance = baumgartner2014_pmv2ppp(p,tang,rang,'QE');      
      varargout{1}=querr;
      meta.ylabel='Quadrant errors (%)';
      par.par1=size(idx,1);
      meta.par1='Number of confusions';
      par.par2=size(m,1);
      meta.par2='Number of responses within the lateral range';
%       par.chance=chance;
%       meta.chance='Chance performance';
      
    case 'querrMiddlebrooksRAU' % as in Middlebrooks (1999) but calculated in rau
      idx=find(abs(m(:,7))<=30);
      m=m(idx,:);
      idx=find((abs(mynpi2pi(m(:,8)-m(:,6))))<=90);  
      querr=100-size(idx,1)/size(m,1)*100;
      varargout{1}=rau(querr,size(m,1));
      meta.ylabel='RAU quadrant errors';
      par.par1=size(idx,1);
      meta.par1='Number of confusions';
      par.par2=size(m,1);
      meta.par2='Number of responses within the lateral range';
    case 'accE'
      varargout{1}=mynpi2pi(circmean(m(:,4)-m(:,2)));
      meta.ylabel='Elevation bias (deg)';
    case 'absaccE'
      varargout{1}=abs(mynpi2pi(circmean(m(:,4)-m(:,2))));
      meta.ylabel='Absolute elevation bias (deg)';
    case 'precE'
      precE=circstd(m(:,4)-m(:,2));      
      varargout{1}=precE;
      meta.ylabel='Elevation precision error (deg)';
    case 'lateral'        % obsolete
      nargout=2;
      varargout{1}=accL;
      varargout{2}=precL;
    case 'polar'          % obsolete
      nargout=3;
      varargout{1}=accP;
      varargout{2}=precP;        
      varargout{3}=querr;
    case 'precPmedianlocal'
      idx=find(abs(m(:,7))<=30);
      m=m(idx,:);
      idx=find(abs(mynpi2pi(m(:,8)-m(:,6)))<90);
      if isempty(idx), error('Quadrant errors of 100%, local bias is NaN'); end
      precPM=circstd(m(idx,8)-m(idx,6));        
      varargout{1}=precPM;
      meta.ylabel='Local central precision error (deg)';
    case 'precPmedian'
      idx=find(abs(m(:,7))<=30);
      m=m(idx,:);
      precPM=circstd(m(:,8)-m(:,6));        
      varargout{1}=precPM;
      meta.ylabel='Central precision error (deg)';
    case 'rmsPmedianlocal'
      idx=find(abs(m(:,7))<=30);
      m=m(idx,:);
      idx=find(abs(mynpi2pi(m(:,8)-m(:,6)))<90);        
      rmsPlocal=sqrt(sum((mynpi2pi(m(idx,8)-m(idx,6))).^2)/size(idx,1));
      varargout{1}=rmsPlocal;
      meta.ylabel='Local central RMS error (deg)';
      % chance performance (FIXME: not working with exp_baumgartner2014!)
%       tang = m(:,6); 
%       rang = min(tang):max(tang); % assume listener knows target range
%       p = ones(length(rang),length(tang)); % equally distributed response probabilities
%       chance = baumgartner2014_pmv2ppp(p,tang,rang,'PE');
%       par.chance = chance;
%       meta.chance = 'Chance performance';
    case 'rmsPmedian'
      idx=find(abs(m(:,7))<=30);
      m=m(idx,:);
      rmsP=sqrt(sum((mynpi2pi(m(:,8)-m(:,6))).^2)/size(m,1));
      varargout{1}=rmsP;
      meta.ylabel='Central RMS error (deg)';
    case 'SCC'
      trgt = m(:,[1 2]);
      resp = m(:,[3 4]);        
      T = pol2dcos(trgt);
      R = pol2dcos(resp);        
      sxy = T' * R;         
      sxx = T' * T;       
      syy = R' * R;               
      scc = det(sxy) / sqrt(det(sxx) * det(syy));
      varargout{1}=scc;
      meta.ylabel='Spatial correlation coefficient';
      
    case 'gainLstats'
      
      [l.b,l.bint,l.r,l.rint,l.stats]=regress(m(:,7),[ones(length(m),1) m(:,5)]);
      varargout{1}=l;
      meta.ylabel='Lateral gain';
      
    case 'gainL'
      
      l = localizationerror(m,'gainLstats');
      varargout{1}=l.b(2);
      meta.ylabel='Lateral gain';
      
    case 'precLregress'
      
      l = localizationerror(m,'gainLstats');
      x = -90:90;
      y = l.b(2)*m(:,5) + l.b(1); % linear regression
      dev = m(:,7) - y; % deviation
      dev = dev(abs(dev)<=45); % include only quasi-veridical responses
      prec = rms(dev);
      varargout{1} = prec;
      meta.ylabel = 'Lateral scatter (deg)';
     
    case 'pVeridicalL'
      tol = 45; % tolerance of lateral angle error for quasi-veridical responses
      l = localizationerror(m,'gainLstats');
      x = -90:90;
      y = l.b(2)*m(:,5) + l.b(1); % linear regression
      dev = m(:,7) - y; % deviation wrapped to +-180 deg
      prop = sum(abs(dev)<=tol)/length(dev);
      varargout{1} = prop*100;
      meta.ylabel = '% lateral quasi-veridical';
     
    case 'sirpMacpherson2000'
      delta = 40; % outlier tolerance in degrees
      Nmin = 5; % minimum number of responses used for regression
      maxiter = 100; % max number of iterations

      % Only responses for targets located within 30 degrees of the median plane are included 
      m = m( abs(m(:,5))<=30,: );
      
      for frontback = 1:2 % Analyze front and back hemispheres separately
        
        if frontback == 1 % Front
          mhemi = m( m(:,6)<=90 ,:); % frontal target locations
          idinit = mhemi(:,8)<=90; % initialization with correct responses
        else % Back
          mhemi = m( m(:,6)>=90 ,:); % rear targets
          idinit = mhemi(:,8)>=90;  % initialization with correct responses
        end

        if length(idinit)<Nmin
          r.b = nan(1,2);
          r.stats = nan(1,4);
        else
          iditer = false(length(idinit),maxiter);
          y = mhemi(:,8); % response
          x = mhemi(:,6); % target
          for iter = 1:maxiter
            B=regress(y(idinit),[ones(sum(idinit),1) mhemi(idinit,6)]);
            yhat=B(2)*x+B(1);
            dev = mod( (y-yhat) +180,360) -180; % deviation wrapped to +-180 deg
            iditer(:,iter)= abs(dev) < delta; 
            if isequal(iditer(:,iter),idinit)
              break % because it converged
            else
              idinit = iditer(:,iter);
            end
          end
          [nselect,iopt] = max(sum(iditer));
          idselect = iditer(:,iopt);

          [r.b,r.bint,r.r,r.rint,r.stats]=regress(...
            mhemi(idselect,8),[ones(nselect,1) mhemi(idselect,6)]);

        end
        varargout{frontback}=r;
      end

      par='SIRP results';
      
    case 'gainPfront'
      
      f = localizationerror(m,'sirpMacpherson2000');
      varargout{1}=f.b(2);
      meta.ylabel='Frontal polar gain';
      
    case 'gainPrear'
      
      [tmp,r] = localizationerror(m,'sirpMacpherson2000');
      varargout{1}=r.b(2);
      meta.ylabel='Rear polar gain';
      
    case 'gainP'
      
      [f,r] = localizationerror(m,'sirpMacpherson2000');
      varargout{1}=(f.b(2)+r.b(2))/2;
      meta.ylabel='Polar gain';
      
    case 'slopePfront'
      
      g = localizationerror(m,'gainPfront');
      varargout{1} = rad2deg(acos(1./sqrt(g.^2+1)));
      meta.ylabel='Frontal polar regression slope';
      
    case 'slopePrear'
      
      g = localizationerror(m,'gainPrear');
      varargout{1} = rad2deg(acos(1./sqrt(g.^2+1)));
      meta.ylabel='Rear polar regression slope';
      
    case 'slopeP'
      
      g = localizationerror(m,'gainP');
      varargout{1} = rad2deg(acos(1./sqrt(g.^2+1)));
      meta.ylabel='Polar regression slope';
      
    case 'pVeridicalPfront'
      
      tol = 45; % tolerance of polar angle error for quasi-veridical responses
      
      f = localizationerror(m,'sirpMacpherson2000');
      
      mf = m( abs(m(:,5))<=30 & m(:,6)< 90 ,:); % frontal central data only
      
      x = -90:270;
      yf=f.b(2)*mf(:,6)+f.b(1); % linear prediction for flat, frontal targets
      
      devf = mod( (mf(:,8)-yf) +180,360) -180; % deviation wrapped to +-180 deg
      
      pVer = sum(abs(devf) <= tol) / length(devf);
      
      varargout{1}=pVer*100;
      meta.ylabel = '% quasi-veridical (frontal)';
     
    case 'pVeridicalPrear'
      
      tol = 45; % tolerance of polar angle error for quasi-veridical responses
      
      [~,r] = localizationerror(m,'sirpMacpherson2000');
      
      mr = m( abs(m(:,5))<=30 & m(:,6)>=90 ,:); % rear central data only
      
      x = -90:270;
      yr = r.b(2)*mr(:,6) + r.b(1); % linear prediction for flat, rear targets
      
      devr = mod( (mr(:,8)-yr) +180,360) -180; % deviation wrapped to +-180 deg
      
      pVer = sum(abs(devr) <= tol) / length(devr);
      
      varargout{1}=pVer*100;
      meta.ylabel = '% quasi-veridical (rear)';
     
    case 'pVeridicalP'
      
      tol = 45; % tolerance of polar angle error for quasi-veridical responses
      [f,r] = localizationerror(m,'sirpMacpherson2000');
      
      mf = m( abs(m(:,5))<=30 & m(:,6)< 90 ,:); % frontal central targets
      mr = m( abs(m(:,5))<=30 & m(:,6)>=90 ,:); % rear central targets
      
      x = -90:270;
      yf = f.b(2)*mf(:,6) + f.b(1); % linear prediction for flat, frontal targets
      yr = r.b(2)*mr(:,6) + r.b(1); % linear prediction for flat, rear targets
      
      devf = mod( (mf(:,8)-yf) +180,360) -180; % deviation wrapped to +-180 deg
      devr = mod( (mr(:,8)-yr) +180,360) -180; % deviation wrapped to +-180 deg
      
      pVerf = sum(abs(devf) <= tol) / length(devf);
      pVerr = sum(abs(devr) <= tol) / length(devr);
      
      pVer = (pVerf + pVerr)/2;
      
      varargout{1}=pVer*100;
      meta.ylabel = '% quasi-veridical';
     
    case 'precPregressFront'
      
      tol = 45; % tolerance of polar angle error for quasi-veridical responses
      
      f = localizationerror(m,'sirpMacpherson2000');
      
      mf = m( abs(m(:,5))<=30 & m(:,6)< 90 ,:); % frontal central targets only
      
      x = -90:270;
      yf = f.b(2)*mf(:,6) + f.b(1); % linear prediction for flat, frontal targets
      
      devf = mod( (mf(:,8)-yf) +180,360) -180; % deviation wrapped to +-180 deg
      
      dev = devf(abs(devf)<=tol); % include only quasi-veridical responses
      prec = rms(dev);
      varargout{1} = prec;
      meta.ylabel = 'Frontal polar scatter (deg)';
     
    case 'precPregressRear'
      
      tol = 45; % tolerance of polar angle error for quasi-veridical responses
      
      [~,r] = localizationerror(m,'sirpMacpherson2000');
      
      mr = m( abs(m(:,5))<=30 & m(:,6)>=90 ,:); % rear central targets only
      
      x = -90:270;
      yr = r.b(2)*mr(:,6) + r.b(1); % linear prediction for flat, rear targets
      
      devr = mod( (mr(:,8)-yr) +180,360) -180; % deviation wrapped to +-180 deg
      
      dev = devr(abs(devr)<=tol); % include only quasi-veridical responses
      prec = rms(dev);
      varargout{1} = prec;
      meta.ylabel = 'Rear polar scatter (deg)';
      
    case 'precPregress'
      
      precF = localizationerror(m,'precPregressFront');
      precR = localizationerror(m,'precPregressRear');
      prec = (precF+precR)/2;
      varargout{1} = prec;
      meta.ylabel = 'Polar scatter (deg)';
     
    case 'perMacpherson2003'
      
      if isempty(kv.f) || isempty(kv.r)
        amt_disp('Regression coefficients missing! Input is interpreted as baseline data!','documentation');
        amt_disp('For this analysis the results from `sirpMacpherson2000()` for baseline data','documentation');
        amt_disp('are required and must be handled as `localizationerror(m,f,r,...)`.','documentation');
        [kv.f,kv.r] = localizationerror(m,'sirpMacpherson2000');
      end
      
      plotflag = false;
      tol = 45; % tolerance of polar angle to be not counted as polar error

      mf = m( abs(m(:,5))<=30 & m(:,6)< 90 ,:); % frontal central targets only
      mr = m( abs(m(:,5))<=30 & m(:,6)>=90 ,:); % rear central targets only
      
      x = -90:270;
      yf=kv.f.b(2)*mf(:,6)+kv.f.b(1); % linear prediction for flat, frontal targets
      yr=kv.r.b(2)*mr(:,6)+kv.r.b(1); % linear prediction for flat, rear targets
      
      devf = mod( (mf(:,8)-yf) +180,360) -180; % deviation wrapped to +-180 deg
      devr = mod( (mr(:,8)-yr) +180,360) -180;
      
      f.pe = sum(abs(devf) > tol);
      r.pe = sum(abs(devr) > tol);

      per = (f.pe + r.pe) / size(m,1) * 100; % polar error rate
      varargout{1} = per;
      meta.ylabel = 'Polar error rate (%)';
      
      if plotflag 
        plot(m(:,6),m(:,8),'ko');
        axis equal
        axis tight
        hold on
        plot(mf(:,6),yf,'LineWidth',2);
        plot(mr(:,6),yr,'LineWidth',2);
        plot(mf(:,6),yf+45,'r','LineWidth',2);
        plot(mf(:,6),yf-45,'r','LineWidth',2);
        plot(mr(:,6),yr+45,'r','LineWidth',2);
        plot(mr(:,6),yr-45,'r','LineWidth',2);
        title(['PER: ' num2str(per,2) '%'])
      end
       
    otherwise
      error(['Unknown ERR: ' err]);
  end
  if length(varargout) < 2
    varargout{2}=meta;
  end
  if length(varargout) < 3
    varargout{3}=par; 
  end
end

end

function out_deg=mynpi2pi(ang_deg)
ang=ang_deg/180*pi;
out_rad=sign(ang).*((abs(ang)/pi)-2*ceil(((abs(ang)/pi)-1)/2))*pi;
out_deg=out_rad*180/pi;
end

function [phi,mag]=circmean(azi,dim)
% [PHI,MAG]=circmean(AZI,DIM)
%
% Calculate the average angle of AZI according
% to the circular statistics (Batschelet, 1981).
% 
% For vectors, CIRCMEAN(AZI) is the mean value of the elements in AZI. For
% matrices, CIRCMEAN(AZI) is a row vector containing the mean value of
% each column.  For N-D arrays, CIRCMEAN(azi) is the mean value of the
% elements along the first non-singleton dimension of AZI.
%
% circmean(AZI,DIM) takes the mean along the dimension DIM of AZI. 
%
% PHI is the average angle in deg, 0..360�
% MAG is the magnitude showing the significance of PHI (kind of dispersion)
%   MAG of zero shows that there is no average of AZI.
%   MAG of one shows a highest possible significance
%   MAG can be represented in terms of standard deviation in deg by using
%       rad2deg(acos(MAG)*2)/2
%
% [PHI,MAG]=circmean(AZI,W) calculates a weighted average where each sample
% in AZI is weighted by the corresponding entry in W.
%
% Piotr Majdak, 6.3.2007

if isempty(azi) % no values
  phi=NaN;
  mag=NaN;
  return;
end

if prod(size(azi))==1 % one value only
  phi=azi;
  mag=1;
  return;
end

if ~exist('dim')  % find the dimension
  dim = min(find(size(azi)~=1));
  w=ones(size(azi));
elseif prod(size(dim))~=1
    % dim vector represents weights
  w=dim;
  dim = min(find(size(azi)~=1));
end


x=deg2rad(azi);

s=mean(sin(x).*w,dim)/mean(w);
c=mean(cos(x).*w,dim)/mean(w);

y=zeros(size(s));
for ii=1:length(s)
  if s(ii)>0 & c(ii)>0
    y(ii)=atan(s(ii)/c(ii));
  elseif c(ii)<0
    y(ii)=atan(s(ii)/c(ii))+pi;
  elseif s(ii)<0 & c(ii)>0
    y(ii)=atan(s(ii)/c(ii))+2*pi;
  elseif s(ii)==0
    y(ii)=0;
  else
    y(ii)=NaN;
  end
end
mag=sqrt(s.*s+c.*c);
phi=wrapTo360(rad2deg(y));
end

function mag=circstd(azi,dim)

% MAG=circstd(AZI,DIM)
%
% Calculate the standard deviation angle (in deg) of AZI (in deg) according
% to the circular statistics (Batschelet, 1981).
% 
% For vectors, CIRCSTD(AZI) is the std of the elements in AZI. For
% matrices, CIRCSTD(AZI) is a row vector containing the std of
% each column.  For N-D arrays, CIRCSTD(azi) is the std value of the
% elements along the first non-singleton dimension of AZI.
%
% CIRCSTD(AZI,DIM) takes the std along the dimension DIM of AZI. 
%
%
% MAG=CIRCSTD(AZI,W) calculates a weighted std where each vector strength
% of AZI is weighted by the corresponding entry in W.
%
% Piotr Majdak, 21.08.2008

if isempty(azi) % no values
  mag=NaN;
  return;
end

if prod(size(azi))==1 % one value only
  mag=0;
  return;
end

if ~exist('dim')  % find the dimension
  [phi,mag]=circmean(azi);
elseif prod(size(dim))~=1
  [phi,mag]=circmean(azi,dim);
end

mag=rad2deg(sqrt(2*(1-mag)));

end


