function E = synth_IVIM(params,grads)
%
% Description: synthesise the dMRI signal for IVIM model
%
% Output:   
% E - diffusion signal
%
% Input:
% params - size 2 vector of model parameters for IVIM model:
%       params(1) - S0: proton density
%       params(2) - f: perfusion fraction
%       params(3) - dp: pseudo diffusion coefficient
%       params(4) - d: apparent diffusion coefficient
%
% grads - gradient table (in form [gx gy gz b]) to synthesise the signal at 
%
%
% Author:
%   Paddy Slator (p.slator@ucl.ac.uk)

%only need the b-values
bvals = grads(:,4);





end