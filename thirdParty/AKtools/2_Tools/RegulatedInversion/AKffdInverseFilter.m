% s = AKffdInverseFilter(s)
% to be called from AKregulatedInversion.m
% see AKregulatedInversionDemo.m for documentation
%
% Calculate regularized LMS-inverse of a multichannel system response
% with regard to target functions of given magnitude and selectable phase
% (min, lin, zero) using fast frequncy deconvolution (see Norcross 2006, AES)
%
% Regularization and target functions have to be delivered as a,b-filter 
% coeffs, and with beta-weigths, respectively.
%
% Alexander Lindau, Fabian Brinkmann
% TU Berlin, audio communication group, 2008


% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

% for documentation only (old, as struct' s' is used for parameter input):
%
% function  [h,a_eq_bp_phase] = getFFDInverseFilter.m(x,fs,type,margin,b_REG,a_REG,beta,APPLY_REG,a_BP,b_BP,APPLY_BANDPASS,TITLE,do_plot)
%
% required args:
%     x              - input vector to be inverted
%     fs             - samplingrate, default: 1
%     type           - 'min_phase','lin_phase','zero_phase', default: 'min_phase'
%     margin         - to protect preringing of min_phase inverse filter [samples], default: 0
%     b_REG,a_REG    - coeffs for spectral weigthing of regularization,default: linear, i.e.  b_REG,a_REG mimic a length-N dirac pulse FIR Filter
%     beta           - overall gain of regularization, default: 0.01
%     APPLY_REG      - [0/1], use regularization, default: 1
%     a_BP,b_BP      - coeffs for band limitation of inversion target (assumes min phase!), default: linear
%     APPLY_BANDPASS - [0/1], apply bandpass on inversion target(select min-/lin-/zerophase bandpass as necessary), default = 0
%     TITLE          - title of plot, default: 'Norcross06 inverse filter'
%     do_plot        - [0/1], simple results plot, default = 0,
%     fade_in        - apply fade-in on norcross-filer, in samples, default = 0
%     fade_out       - apply fade-out on norcross-filer in samples, default = 0
%
%
% issues:
% TBD: to be changed to attribute-value-style in final version
% TBD: length of target function limited to s.Nsamples (windowing artifacts?)
% TBD: zero phase not tested again
% TBD: min_phase inversion is approx. 7 dB louder than lin_phase
%-------------------------------------------------------------------------%
function s = AKffdInverseFilter(s)

for ch = 1:s.numChannel
    
    %% 1.1 create target function magnitude response
    % decide about applying freq.-dependent regularization
    if ~s.APPLY_REG
        s.reg_beta=zeros(1,s.numChannel);
    end
       
    % get spectrum of regularization function
    A_r=fft(s.regularization(:,ch));
    
    % precalculations for Norcross 2006, formula 2.14
    % |A_r|^2
    A_r_abs_sq=A_r.*conj(A_r);
    % |X|^2
    C=fft(s.avg(:,ch));
        
    %% 1.2 limit dynamic range of system function below pass band 
    if s.dyn % bypass with s.dyn  = 0
        % limit dynamic range of target function
        abs_C = abs(C);
        abs_C_dyn = max(abs_C,10^(-abs(s.dyn)/20));
        angle_C = angle(C);
        C = abs_C_dyn.*exp(1i*angle_C);
        fprintf('\ngetFFDInverseFilter: Dynamic range of system function was limited to -%4.1f dB.\n',abs(s.dyn))
    end
   
    C_abs_sq=C.*conj(C);
    
    % calculate A_eq acc. to Norcross 2006 formula 2.14
    A_eq=1./( 1 + s.reg_beta(ch) * ( A_r_abs_sq ./ C_abs_sq ));
    
    % apply target function on target A_eq (or not)
    if s.target
        % apply target function by magnitude multiplication
        A_eq_bp = A_eq.*abs(fft(s.target));
    else
        A_eq_bp = A_eq;
        s.target    = zeros(s.Nsamples,1);
        s.target(1) =1;
    end
    
    %% 2.1 synthesize zero/lin/min-phase for target function magnitude response (for realizing formula 2.13 in Norcross 2006)
    disp(['Applying ' s.phase_type ' to target function'])
    s.a_eq_bp_phase(:,ch) = AKphaseManipulation(ifft(A_eq_bp, 'symmetric'), s.fs, s.phase_type, s.Nfft_double);
    A = fft(s.a_eq_bp_phase(:,ch));
    
    %% 2.2 limit dynamic range of target function below pass band 
    if s.dyn % bypass with s.dyn  = 0
        % limit dynamic range of target function
        abs_A = abs(A);
        abs_A_dyn = max(abs_A,10^(-abs(s.dyn)/20));
        angle_A = angle(A);
        A = abs_A_dyn.*exp(1i*angle_A);
        fprintf('\ngetFFDInverseFilter: Dynamic range of target function was limited to -%4.1f dB.\n',abs(s.dyn))
    end
   
    
    %% 3. calculate inverse filter acc. to Norcross 2006 formula 2.13
    H=( conj(C).*A)./C_abs_sq;
    % apply filters for return/plotting
    h=real(ifft(H));
    % compensate preringign by circshifting
    if strcmpi(s.phase_type, 'min_phase') && s.margin
        h=circshift(h,s.margin);
    end
    % fade in/out
    h = AKfade(h, [], s.fade_in, s.fade_out);
   
    s.compensation_filter(:,ch) = h;
       
end