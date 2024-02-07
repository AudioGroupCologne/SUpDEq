function y = smalt2014_efferentinnervation(x,CF)
%   Usage: y = smalt2014_efferentinnervation(x,CF) efferent innervation function
%
%   Input parameters:
%     x     : gamma function parameters; x(1) = k, x(2) = theta
%     CF    : frequency in kHz
%
%   Output parameters:
%     y     : efferent innervation function
%
%
%   for the model by Smalt, Heinz and Strickland (2014)
%   based on (Zilany and Bruce (JASA 2006, 2007) Auditory Nerve Model)
%
%   See also: smalt2014
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/smalt2014_efferentinnervation.php



k = x(1);
theta = x(2);

% define gamma pdf function
gammapdffunc = @(x)-1/(gamma(k)*theta.^(k))*x.^(k-1).*exp(-x./theta);


% Now get CF of peak in gamma function by solving derivivate function
% syms x;
% CF_peak = solve(-theta^(-k-1).*x.^(k-2).*exp(-x/theta).*(theta-theta*k+x)./gamma(k));

% or instead find peak by itertive method (faster) between 0 and 50 kHz
options = optimset('Algorithm','active-set','Display','off');
CF_peak = fmincon(gammapdffunc,8,[],[],[],[],0, 50,[],options);

% get value of function
gamma_peak = gammapdffunc(CF_peak);

% return normalized gamma function value
y = gammapdffunc(CF)./gamma_peak;
