function [XR, XR_cov, dtR] = LS_dtR(XS, PS, dtR, XR0)
% SYNTAX: [XR, dtR] = leastSquare(sat, PS)

% Description:
%     Leat-Square solution to the linearized GPS observation equations.
%     Provides estimate of receiver position and clock bias.

%  INPUTS:
%     XS  = satellite coordinates
%     PS  = pseudorange (corrected for satellite clock error)
%     dtR = satellite clock error initial estimate
%     XR0 = reciever position initial estimate
  
%  OUTPUTS:
%     XR  = receiver ECEF coordinates (vector)
%     XR_cov = covariance matrix
%     dtR = receiver clock bias (in meters)
 
% ------------------------------------------------------------------------
%  This program uses some functions from GoGPS MATLAB.
% ------------------------------------------------------------------------
%  Copyright 2016, Hadi Tabatabaee, All rights reserved.
% ------------------------------------------------------------------------


nsat = length(PS); % total number of satellites
if ~any(XR0)
    XR0 = [0; 0; 0]; % initial guess for receiver position
end
if ~any(dtR)
    dtR = 0; % receiver clock bias initialization
end
b = [0; 0; 0; 0]; % residual initialization
PS0 = zeros(nsat, 1);
A = zeros(nsat, 4);

iter = 0;
maxiter = 10;

while iter <= maxiter
    
    for i = 1:nsat
        PS0(i) = (sqrt((XS(:,i) - XR0)'*(XS(:,i) - XR0))) - dtR; % the computed pseudorange, note that dtR has been subtracted.
        b(i) = PS(i) - PS0(i);
        A(i,:) = [(XR0(1) - XS(1,i)) ./ (PS(i)-dtR), ... %column for X coordinate
                  (XR0(2) - XS(2,i)) ./ (PS(i)-dtR), ... %column for Y coordinate
                  (XR0(3) - XS(3,i)) ./ (PS(i)-dtR), ... %column for Z coordinate
                  1];
    end
    
    XR = A\b;
    %XR = inv(A'*A)*A'*b;
    %XR = (A'*A)\(A'*b);
    XR0 = XR(1:3) + XR0;
    dtR = XR(4);
    iter = iter + 1;
end
XR = XR0;
dtR = dtR/goGNSS.V_LIGHT;
cov_XYZ = (A'*A)^-1;
cov_XYZ = cov_XYZ(1:3,1:3);
%cov_ENU = global2localCov(cov_XYZ, XR);
XR_cov = 36*cov_XYZ;
PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
%HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
%VDOP = sqrt(cov_ENU(3,3));
end
