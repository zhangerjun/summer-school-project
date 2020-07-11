function E = synth_quadexp(params,grads)
%
% Description: synthesise the dMRI signal for quad-exponential model
%
% Output:   
% E - diffusion signal
%
% Input:
% params - size 2 vector of model parameters in SI units for IVIM model:
%       params(1) - S0: proton density
%       params(2) - f1: first volume fraction
%       params(3) - f2: second volume fraction
%       params(4) - f3: third volume fraction
%       params(5) - d1: first diffusion coefficient
%       params(6) - d2: second diffusion coefficient
%       params(7) - d3: third diffusion coefficient
%       params(8) - d4: fourth diffusion coefficient
%
% bvals - b-values to synthesise the signal at 
%
%
% Author:
%   Paddy Slator (p.slator@ucl.ac.uk)

%get the b-values
bvals = grads(:,4);


S0=params(1);
f1=params(2);
f2=params(3);
f3=params(4);
d1=params(5);
d2=params(6);
d3=params(7);
d4=params(8);

E = S0*( f1 * exp(-bvals.*d1) ...
    + f2 * exp(-bvals.*d2) ...
    + f3 * exp(-bvals.*d3) ...
    + (1 - f1 - f2 - f3) * exp(-bvals.*d4));  





end