function [XR, XR_cov, dtR] = KF_simple2(XS, pr)

global state state_cov

if ~any(state)
    %[XR, XR_cov, dtR] = DG_sol(XS, pr);
    [XR, XR_cov, dtR] = leastSquare2(XS, pr);
    state = XR;
    state_cov = XR_cov;
else
    % new measurements:
    %[z, R, dtR] = DG_sol(XS, pr);
    [z, R, dtR] = leastSquare2(XS, pr)
    
    % add process noise to state covariance estimate:
    state_cov = state_cov + 20*eye(3);
    
    % Kalman gain:
    K = state_cov*pinv(state_cov + R);
    
    % state update:
    state = state + K*(z - state);
    state_cov = (eye(3) - K)*state_cov;

end
