function loglik = RicianLogLik(meas, signals, sig)
% 
% camino.m--------------------------------------------------------------
% Log likelihood with Rician noise
%% 
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
%
%

%note - you can use logbesseli0 to calculate the Bessel function term