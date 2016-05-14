function [XR, dtR, A, b] = leastSquare(sat, XS, PS)
% SYNTAX: [XR, dtR] = leastSquare(sat, PS)

% Description:
%     Leat-Square solution to the linearized GPS observation equations.
%     Provides estimate of receiver position and clock bias.

%  INPUTS:
%     sat = available satellite indexes
%     XS  = satellite coordinates
%     PS  = pseudorange (corrected for satellite clock error)
  
%  OUTPUTS:
%     XR  = receiver ECEF coordinates (vector)
%     dtR = receiver clock bias (in meters)
 
% ------------------------------------------------------------------------
% 
%  Copyright 2016, Hadi Tabatabaee, All rights reserved.
%  
% ------------------------------------------------------------------------


nsat = length(sat); % total number of satellites

k = find(~(PS(:)==0)); % find pseudorange indexes that are non-zero
PS1 = PS(k);
PS = PS1;

XR0 = [0; 0; 0]; % initial guess for receiver position
dtR = 0; % receiver clock bias initialization
b = [0; 0; 0; 0]; % residual initialization
PS0 = zeros(nsat, 1);
A = zeros(nsat, 4);

iter = 0;
maxiter = 10;

while iter <= maxiter
    
    for i = 1:nsat
        PS0(i) = (sqrt((XS(i,:)' - XR0)'*(XS(i,:)' - XR0))) - dtR; % the computed pseudorange, note that dtR has been subtracted.
        b(i) = PS(i) - PS0(i);
        for j = 1:3
            A(i,j) = -(XS(i,j) - XR0(j))/(PS(i)-dtR);
        end
        A(i,4) = 1;
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

end
