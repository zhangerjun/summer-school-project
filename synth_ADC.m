function E = synth_ADC(params,grads)
%
% Description: synthesise the dMRI signal for ADC model
%
% Output:   
% E - diffusion signal
%
% Input:
% params - size 2 vector of model parameters in SI units for ADC model:
%       params(1) - S0: proton density
%       params (2) - d: apparent diffusion coefficient
% bvals - b-values to synthesise the signal at 
%
%
% Author:
%   Paddy Slator (p.slator@ucl.ac.uk)

%only need the b-values
bvals = grads(:,4);








end