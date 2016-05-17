function [X, P] = kalman_PV(X_p, P_p, Z, Q)
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

% PLACEHOLDER