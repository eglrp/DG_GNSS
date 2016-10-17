function [XR, XR_cov, dtR] = kalmanFilter_P(XS, pr)
global state state_cov

% State Projection using P-model: Just use the state variable from last
% time

% Defining needed variables:
n = length(pr);
Ar = zeros(n);
hr = ones(n,1);
c = 299792458; % speed of light.
P = eye(n) - 1/n*(hr*hr');

% Creating distance geometry reference matrix
for i = 1:n
    for j = 1:n
        Ar(i,j) = (XS(:,i) - XS(:,j))'*(XS(:,i) - XS(:,j));
    end
end
    
% form u & r vectors
u = diag(pr(1:10)*pr(1:10)');
r = pr(1:10);
    
% weighting constants
tol = 25;
k_const = 100;

if ~any(state)
    % Just run a simple GPS algorithm to get first state estimate.
    % Ar*x = u solution
    x_ur = pinv(P*Ar + k_const*(hr*hr'), tol)*(P*u + k_const*hr);
    x_ur_last = (1/n)*hr'*(u - Ar*x_ur);

    x_u = [x_ur; x_ur_last];
    
    % Ar*x = r solution
    x_rr = pinv(P*Ar + k_const*(hr*hr'), tol)*(P*r+ k_const*hr);
    x_rr_last = (1/n)*hr'*(r - Ar*x_rr);

    x_r = [x_rr; x_rr_last];
    
    % clock bias calculation
    if x_r'*[r;0] < 0
        cdt = (x_u'*[r;0] + x_r'*[u;1] + sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));
    else
        cdt = (x_u'*[r;0] + x_r'*[u;1] - sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));
    end
    % calculate final x vector
    x_r = x_ur - 2*cdt*x_rr;
    
    XR = XS*x_r;
    dtR = cdt/c;
    
   %calculate DOP and covariance 
   %figure out how to calculate cov of state from previous cov
   [PDOP, HDOP, VDOP, cov_XYZ, cov_ENU] = DOP(XR, XS);
    G = XS'*pinv(XS*XS'); %least-norm solution
    state = x_r;
    state_cov = G*6*cov_XYZ*G';
    XR_cov = 6*cov_XYZ;
else
    % Measurement update:
    z = [u; r];
    
    % Measurement covariance mess:
    [~, ~, ~, cov_XYZ, ~] = DOP(XS*state, XS);
    G = XS'*pinv(XS*XS'); %least-norm solution
    
    % State projection -- add process noise
    state_cov = state_cov + G*36*eye(3)*G';
    
    % M matrix
    M = pinv(P*Ar, tol)*P;
    
    % dtR calculations:
    dtR = clock_bias(M,Ar,u,r)/c;
    
    % H matrix:
    H = [M, -2*c*dtR*M];
    
    % Measurement covariance matrix calculation:
    J = (XS*H)'*pinv((XS*H)*(XS*H)');
    R = J*6*cov_XYZ*J'; % Measurement covariance matrix
         
    % Kalman gain calculation:
    K = state_cov*pinv(state_cov + H*R*H');
    
    % State estimate:
    state = state + K*(H*z - state);
    state_cov = (eye(length(K))-K)*state_cov;
    
    % Position estimate and covariance
    XR = XS*state;
    XR_cov = XS*state_cov*XS';
end






