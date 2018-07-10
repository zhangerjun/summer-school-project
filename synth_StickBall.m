function E = synth_StickBall(params,grads)

%
% Description: synthesise the dMRI signal for stick-ball model
%
% Output:   
% E - diffusion signal
%
% Input:
% params - size 8 vector of model parameters for stick-zeppelin model:
%       params(1) - S0
%       params(2) - fv: perfusion fraction
%       params(3) - dv: pseudo diffusivity 
%       params(4) - di: diffusivity 
%       params(5) - theta: angle from the z direction of the main direction of the tensor 
%       params(6) - phi: azimuthal angle of the main direction of the tensor 
%
% grads - gradient table (in form [gx gy gz b]) to synthesise the signal at 
%
%
% Author:
%   Paddy Slator (p.slator@ucl.ac.uk)
%
%

S0 = params(1);
f=params(2);%stick volume fraction 
dp=params(3);%stick diffusivity
d=params(4);%diffusivity
theta=params(5);
phi=params(6);


%get the b-values and gradient directions
bvals = grads(:,4);
bvecs = grads(:,1:3)';

%Synthesize the signals

%get the fibre direction
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(bvecs.*repmat(fibdir, [length(bvecs) 1])')';

E = S0*(f*exp(-bvals*dp.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*d));







