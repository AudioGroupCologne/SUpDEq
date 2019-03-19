% radial = AKshRadial(kr, type, kind, n, derived)
% calulates spherical radial functions
%
% e.g.
% AKshRadial(0:.1:10, 'bessel', 1, 0:6, false)
% computes the spherical bessel function of first kind and order zero to
% six between kr=0 and kr=10
%
%
% I N P U T:
% kr      - x values where radial functions are calculated at. This is
%           usually the product of the wavenumber k and the radius r in
%           meter. Can be a scalar or vector
% type    - 'hankel', or 'bessel' so calculate spherical Hankel or Bessel
%           fucntions
% kind    - 1 or 2 to calculate functions of first or second kind
% n       - order of spherical radial function. Can be a scalar or vector
%           of positive integers including zero.
% derived - calculate the derivate of the radial function. true or false
%           (default = false)
%
% O U T P U T:
% radial  - radial function of size [numel(kr) numel(n)]
%
% R E F E R E N C E S
% [1] http://keisan.casio.com/has10/SpecExec.cgi?lang=en&id=system/2006/1222521303
% [2] Lawrence Ziomek: Fundamentals of acoustic field theory and space-time
%     signal processing. Frist edition. CRC Press, Boca Raton, Ann Arbor,
%     London, Tokyo. 1995.
% [3] P M Morse,  K U Ingard: Theoretical acoustics. Princeton University
%     Press, Princeton. First Edition with errata page. 1986.
% [4] Boaz Rafaely: Fundamentals of sperical array processing. First
%     edition. Springer, Berlin, Heidelberg, Germany. 2015.
%
% 07/2010  frank.schultz@tu-berlin.de (initial dev.)
% 12/2016  fabian.brinkmann@tu-berlin.de (vectorization)

function radial = AKshRadial(kr, type, kind, n, derived)

% set default parameter
if ~exist('derived', 'var')
    derived = false;
end

% format input
type    = lower(type);
type(1) = upper(type(1));
kind    = num2str(kind);
if derived
    derived = 'Derived';
else
    derived = '';
end

% calculate the radial functions
if numel(kr)>1 && numel(n)>1
    
    radial = zeros(numel(kr), numel(n));
    
    for nn = 1:numel(n)
        radial(:,nn) = eval(['AKsh' type kind derived '(kr, n(nn))']);
    end
    
else
    radial = eval(['AKsh' type kind derived '(kr, n)']);
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% radial functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Frank Schultz, FG Audiokommunikation, TU Berlin
%frank.schultz@tu-berlin.de, +49 175 15 49 763, Skype: j0shiiv
% 1.00 07/28/2010 ok

% spherical bessel function of first kind x = j_Ord(kr)
function [x] = AKshBessel1(kr,Ord)
    %(4.1-36) in Ziomek1995, 7.2.12 in Morse1986:
    x = sqrt(pi./(2.*kr)).*besselj(Ord+1/2,kr); %for kr!=0
end


% spherical bessel function of second kind (i.e. neumann) x = y_Ord(kr)
function [x] = AKshBessel2(kr,Ord)
    %(4.1-37) in Ziomek1995, 7.2.12 in Morse1986:
    x = sqrt(pi./(2.*kr)).*bessely(Ord+1/2,kr); %for kr!=0
end

% derivative of spherical bessel function of first kind x = j'_Ord(kr)
function [x] = AKshBessel1Derived(kr,Ord)
    %(4.1-51) in Ziomek1995:
%      x = sph_bessel_1st(kr,Ord-1) - ...
%          sph_bessel_1st(kr,Ord)*(Ord+1)./kr;  
    %7.2.13 in Morse1987:      
    x = (...
        Ord     .* AKshBessel1(kr,Ord-1) - ...
        (Ord+1) .* AKshBessel1(kr,Ord+1) ...
        )...
        ./(2*Ord+1); %for kr!=0              
end

% derivative of spherical bessel function of second kind (i.e. neumann) x = y'_Ord(kr)
function [x] = AKshBessel2Derived(kr,Ord)
    %(4.1-51) in Ziomek1995
%      x = sph_bessel_2nd(kr,Ord-1) - ...
%          sph_bessel_2nd(kr,Ord)*(Ord+1)./kr;
    %7.2.13 in Morse1987:      
    x = (...
        Ord     .* AKshBessel2(kr,Ord-1) - ...
        (Ord+1) .* AKshBessel2(kr,Ord+1) ...
        )...
        ./(2*Ord+1); %for kr!=0       
end

% spherical hankel function of first kind x = h^1_Ord(kr)
function [x] = AKshHankel1(kr,Ord)
    %(4.1-38) in Ziomek1995, 7.2.10 in Morse1986:
    x = AKshBessel1(kr,Ord) +...
        AKshBessel2(kr,Ord) * 1i; %for kr!=0    
end

% spherical hankel function of second kind x = h^2_Ord(kr)
function [x] = AKshHankel2(kr,Ord) %#ok<*DEFNU>
    %(4.1-39) in Ziomek1995:
    x = conj(AKshHankel1(kr,Ord)); %for kr!=0
end

% derivative of spherical hankel function of first kind x = h'^1_Ord(kr)
function [x] = AKshHankel1Derived(kr,Ord)
    %(4.1-38) in Ziomek1995, 7.2.10 in Morse1986:
    x = AKshBessel1Derived(kr,Ord) +...
        AKshBessel2Derived(kr,Ord) * 1i; %for kr!=0                            
end

% derivative of spherical hankel function of second kind x = h'^2_Ord(kr)
function [x] = AKshHankel2Derived(kr,Ord)
    %(4.1-39) in Ziomek1995:
    x = conj(AKshHankel1Derived(kr,Ord)); %for kr!=0
end
