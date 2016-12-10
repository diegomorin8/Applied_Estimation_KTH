% function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           3X1
%           sigma_bar(t)        3X3
%           Q                   2X2
%           z_i(t)              2X1
%           M                   2XN
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1X1
%           outlier             1X1
%           nu^i(t)             2XN
%           S^i(t)              2X2XN
%           H^i(t)              2X3XN
function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)
%Alg 2
%For each observation
 for obs = 1:size(M,2)
        zbar = observation_model(mu_bar,M,obs);  
        H(:,:,obs)=jacobian_observation_model(mu_bar,M,obs,zbar,1);
        S(:,:,obs)=H(:,:,obs)*sigma_bar*H(:,:,obs)'+Q;
        nu(:,obs)=z_i-zbar;
        %As implied in the lab notes
        nu(2,obs)=mod(nu(2,obs)+pi,2*pi)-pi;
        %As it is defined in the lab notes
        D(obs)=(nu(:,obs)'/S(:,:,obs))*nu(:,obs);
        fi(obs)=det(2*pi*S(:,:,obs)).^(-1/2)*exp(-1/2*D(obs));
 end
    %Get the index for the max value of fi.
    c = find(fi == max(fi));
    %Outlier identification
    outlier = (D(c)>=Lambda_m);
end