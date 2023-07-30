function y = randvmf(kappa, mu)
% randvmf Generates a random sample of size 1 from the 
%                von Mises-Fisher distribution
%
% Required Inputs
% kappa       distribution parameter (kappa must be positive)
%
% Optional Inputs
% mu          distribution parameter, (for rotating the North pole to mu)
%
% 2022-04-22 Roberto Barumerli, Acoustic Research Institute, Vienna
%
% This function has been derived from the toolbox of Gy.Terdik, B.Wainwright
% available at https://github.com/TerdikGyorgy/3D-Simulation-Visualization
% That toolbox is licensed under 
% the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/thirdparty/randvmf.php

 
    % Verify parameter values and set defaults
    assert(~isempty(mu))
    assert(kappa > 0)
    
    Np=[0,0,1]; % z-axis (North Pole)
    
    %% density
    % Rubinstein 81, p.39, Fisher 87, p.59
    kappaS=sign(kappa);
    kappa=abs(kappa);
    U = rand(1,1); %random n-by-1 vector from uniform(0,1)
    x=log(2*U*sinh(kappa)+exp(-kappa))/kappa;
    x=kappaS*x;
    
    %%
    psi = 2*pi*rand(1,1);
    s_x = sqrt(1-x.^2); 
    y = [cos(psi).*s_x,sin(psi).*s_x,x];

    %% rotation (orient mean vector (mu) with the North Pole)
    mu=mu/norm(mu); % should be unit vector
    
    if norm(mu-Np) > eps % !!!!!modified!!!!!!!!!!!
        %% rotation (orient mean vector (mu) with the North Pole)
        if mu(3)~= 1
            Ux=cross(Np,mu);%axis of rotation
            Ux= Ux/norm(Ux);
            thetaX=acos(mu(3));
            Rg= rodrigues_rotn(Ux*thetaX);
            y=y*Rg;
        end
    end

end

function M = rodrigues_rotn(vector)
    % Rodrigues' formula for matrix rotation
    % Reference: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    % check Matrix notation section in that page
    % 
    % robaru, Wien
    theta = norm(vector);
    M =  eye(length(vector));

    if (theta < 1e-9)
        return
    end

    v = vector ./ theta;
    v = v(:);

    % anti-asymmetric matrix 
    % https://mathworld.wolfram.com/RodriguesRotationFormula.html
    K = [ 0    -v(3)  v(2); ...
          v(3)  0    -v(1);...
         -v(2)  v(1)  0];

    M = eye(length(v)) + sin(theta)*K + (1-cos(theta))*K^2; % Rodrigues' formula

    M = M';

end    


