% function [mu_bar,sigma_bar] = update(mu_bar,sigma_bar,H_bar,S_bar,nu_bar)
% This function should perform the update process(sequential update).
% You need to make sure that the output sigma_bar is symmetric.
% The last line makes sure that ouput sigma_bar is always symmetric.
% Inputs:
%           mu_bar(t)       3X1
%           sigma_bar(t)    3X3
%           H_bar(t)        2X3
%           S_bar(t)        2X2
%           nu_bar(t)       2X1
% Outputs:
%           mu_bar(t)       3X1
%           sigma_bar(t)    3X3
function [mu_bar,sigma_bar] = update(mu_bar,sigma_bar,H_bar,S_bar,nu_bar)
%Alg 3
%Compute Kalman Gain
K = (sigma_bar*H_bar')/S_bar;
%Update of the state
mu_bar = mu_bar+K*nu_bar;
%We need that K*H and I have the same sizes. 
I = eye(size(K*H_bar));
sigma_bar = (I-K*H_bar)*sigma_bar;
%Update of the covariance matrix
sigma_bar = (sigma_bar + sigma_bar')/2;
end