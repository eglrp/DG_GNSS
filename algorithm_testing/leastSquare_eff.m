function [XR, iter, dtR] = leastSquare_eff(XS, PS)
% SYNTAX: [XR, dtR] = leastSquare(sat, PS)

% Description:
%     Leat-Square solution to the linearized GPS observation equations.
%     Provides estimate of receiver position and clock bias.

%  INPUTS:
%     XS  = satellite coordinates
%     PS  = pseudorange (corrected for satellite clock error)
  
%  OUTPUTS:
%     XR  = receiver ECEF coordinates (vector)
%     XR_cov = covariance matrix
%     dtR = receiver clock bias (in meters)
 
% ------------------------------------------------------------------------
% 
%  Copyright 2016, Hadi Tabatabaee, All rights reserved.
%  
% ------------------------------------------------------------------------


nsat = length(PS); % total number of satellites

XR0 = [0; 0; 0]; % initial guess for receiver position
XR = XR0;
dtR = 0; % receiver clock bias initialization
b = [1; 1; 1; 1]; 
PS0 = zeros(nsat, 1);
A = zeros(nsat, 4);
c = 299792458;

iter = 0;
maxiter = 7;

while (iter <= maxiter)
    
    for i = 1:nsat
        PS0(i) = (sqrt((XS(:,i) - XR0)'*(XS(:,i) - XR0))) - dtR; % the computed pseudorange, note that dtR has been subtracted.
        b(i) = PS(i) - PS0(i);
        A(i,:) = [(XR0(1) - XS(1,i)) ./ (PS(i)-dtR), ... %column for X coordinate
                  (XR0(2) - XS(2,i)) ./ (PS(i)-dtR), ... %column for Y coordinate
                  (XR0(3) - XS(3,i)) ./ (PS(i)-dtR), ... %column for Z coordinate
                  -1];
    end
    
    XR = A\b;
    %XR = inv(A'*A)*A'*b;
    %XR = (A'*A)\(A'*b);
    XR0 = XR(1:3) + XR0;
    dtR = XR(4);
    iter = iter + 1;
end
XR = XR0;
dtR = dtR/c;

end
