function [V, V_cov] = velocitySim(XR, XR_cov, XR0, XR0_cov, dt)
% simulates a velocity for the Kalman filter input based on the previous
% state and the current measurement.

V = (XR - XR0)./dt;
for i = 1:3
    dV(i,1) = sqrt(XR_cov(i,i) + XR0_cov(i,i));
end
V_cov = dV*dV';

end




