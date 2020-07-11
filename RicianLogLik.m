function loglik = RicianLogLik(meas, signals, sig)
% 
% camino.m--------------------------------------------------------------
% Log likelihood with Rician noise
%
% rat = RicianLikeRat(y1, y2, A, sig1, sig2)
% 
% Description: Computes the log likelihood of the measurements given the model signals
% for the Rician noise model and the noise standard deviation.
%  
% Paramaters:
% 
% loglik - log likelihood
% meas - measured data
% signals - signal estimated from the model 
% sig -  standard deviation of the Gaussian distributions underlying
% the Rician distribution.
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Author:
%   Daniel C Alexander (d.alexander@ucl.ac.uk)
%

sumsqsc = (signals.^2 + meas.^2)./(2*sig.^2);
scp = meas.*signals./(sig.^2);
lb0 = logbesseli0(scp);
logliks = - 2*log(sig) - sumsqsc + log(meas) + lb0; % log(signals) changed to log(meas)
loglik = sum(logliks);
