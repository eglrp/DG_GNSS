function [X, Sigma, dtR] = KF_PV_simple(XS, pr)

global state1 state_cov1

if ~any(state1)
    [XR, XR_cov, dtR] = DG_sol(XS, pr);
    state = [XR;0;0;0];
    %state_cov = XR_cov;
    state_cov1 = eye(6); % whatever!!
        X = state;
    Sigma = state_cov1;
    
else
    % new measurements:
    [z, R, dtR] = DG_sol(XS, pr);
    %R = 50*eye(3);
    [V, V_cov] = velocitySim(z, R, state1(1:3), state_cov1(1:3, 1:3), 1);
    
    z = [z;V];
    R = [R, zeros(3,3); zeros(3,3), V_cov];
    
    % add process noise to state covariance estimate:
    
    A = [eye(3), eye(3); zeros(3,3), eye(3)];
    state1 = A*state1; 
    state_cov1 = A*state_cov1*A' + 4*eye(6);
    
    % Kalman gain:
    K = state_cov1*pinv(state_cov1 + R);
    
    % state update:
    state1 = state1 + K*(z - state1);
    state_cov1 = (eye(6) - K)*state_cov1;
    X = state1;
    Sigma = state_cov1;
    
end
