function [X, P] = kalman_PV(X_p, P_p, Z, Q, dt)
% SYNTAX: [X, P] = kalman_PV(X_p, P_p)
%
% INPUT:
%     X_p = apriori state estimate
%     P_p = apriori covariance matrix
%     Z   = measurement input
%     Q   = measurement covariance matrix
%
% OUTPUT:
%     X   = state estimate
%     P   = state covariance matrix
%
% DESCRIPTION:
%     Performs Kalman Filter calculations using a position-velocity state
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
e_vx = 5;
e_vy = 5;
e_vz = 5;
epsilon = [e_x; e_y; e_z; e_vx; e_vy; e_vz];

Q = epsilon*epsilon';

% State transition:

A = [1, 0, 0, dt, 0, 0;
     0, 1, 0, 0, dt, 0;
     0, 0, 1, 0, 0, dt;
     0, 0, 0, 1, 0,  0;
     0, 0, 0, 0, 1,  0;
     0, 0, 0, 0, 0,  1];
 
X_m = A*X_p;
P_m = P_p*A*P_p' + Q;


% Measurement incorporation:
H = A*pinv([XS; zeros(1,3)]);

% Kalman gain calculation:
K = P_m*H'*inv(H*P_m*H' + Q);

% Measurement update stage:
X = X_m + K*(Z - H*X_m);
P = (I - K*H)*P_m;