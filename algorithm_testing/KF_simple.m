function [XR, XR_cov, dtR] = KF_simple(XS, pr)

global state state_cov

if ~any(state)
    [XR, XR_cov, dtR] = DG_sol(XS, pr);
    state = XR;
    %state_cov = XR_cov;
    state_cov = zeros(3,3);
else
    % new measurements:
    [z, R, dtR] = DG_sol(XS, pr);
    %R = 50*eye(3);
    % add process noise to state covariance estimate:
    state_cov = state_cov;
    
    % Kalman gain:
    K = state_cov*pinv(state_cov + R);
    
    % state update:
    state = state + K*(z - state);
    state_cov = (eye(3) - K)*state_cov;
    XR = state;
    XR_cov = state_cov;
    
end
