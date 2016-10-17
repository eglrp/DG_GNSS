function [X, P] = kalman_P(X_p, P_p, Z, Q, XS, A)
% SYNTAX: [X, P] = kalman_P(X_p, P_p)
%
% INPUT:
%     X_p = apriori state estimate
%     P_p = apriori covariance matrix
%     Z   = measurement input
%     Q   = measurement covariance matrix
%     XS  = satellite location matrix
%     A   = distance geometry reference matrix
%
% OUTPUT:
%     X   = state estimate
%     P   = state covariance matrix
%
% DESCRIPTION:
%     Performs Kalman Filter calculations using a position only state
%     transition model for state predictions and the distance geometry
%     method for integrating measurements.
%
% -------------------------------------------------------------------------
%
%   Copyright 2016, Hadi Tabatabaee, All rights reserved.
%
% -------------------------------------------------------------------------

% State transition model covariance matrix:
e_x = 5;
e_y = 5;
e_z = 5;
epsilon = [e_x; e_y; e_z];

Q = epsilon*epsilon';

% State transition:
X_m = X_p;
P_m = P_p + Q;


% Measurement incorporation:
H = A*pinv([XS; zeros(1,3)]);

% Kalman gain calculation:
K = P_m*H'*inv(H*P_m*H' + Q);

% Measurement update stage:
X = X_m + K*(Z - H*X_m);
P = (I - K*H)*P_m;

end

