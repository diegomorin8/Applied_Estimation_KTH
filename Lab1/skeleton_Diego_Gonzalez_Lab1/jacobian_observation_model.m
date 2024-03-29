% function H = jacobian_observation_model(mu_bar,M,j,z_hat,i)
% This function is the implementation of the H function
% Inputs:
%           mu_bar(t)   3X1
%           M           2XN
%           j           1X1 which M column
%           z_hat       2Xn
%           i           1X1  which z_hat column
% Outputs:  
%           H           2X3
function H = jacobian_observation_model(mu_bar,M,j,z_hat,i)
%Eq (6)
H = [ (mu_bar(1)-M(1,j))/z_hat(1,i) (mu_bar(2)-M(2,j))/z_hat(1,i) 0;
    -(mu_bar(2)-M(2,j))/(z_hat(1,i)*z_hat(1,i)) (mu_bar(1)-M(1,j))/(z_hat(1,i)*z_hat(1,i)) -1 ];

end
