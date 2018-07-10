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







end