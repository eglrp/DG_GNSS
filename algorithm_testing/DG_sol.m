function [XR, XR_cov, dtR] = DG_sol(XS, pr)
% This function provides GPS reciever location using the distance geometry
% method for over-determined problems (n>4).
% ------------------------------------------------------------------------
% Input:
%       XS  =  3xn matrix containing satellite coordinates
%       pr  =  nx1 vector containing observed pseudoranges
% Output:
%       XR  =  3x1 vector of receiver location
%       XR_cov = 3x3 covariance matrix for receiver location
%       dtR = receiver clock bias in seconds (scalar)
% ------------------------------------------------------------------------
%   Copyright 2016, Hadi Tabatabaee. All rights reserved.
% ------------------------------------------------------------------------

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
u = diag(pr*pr');
r = pr;
    
% weighting constants
tol = 25;
k_const = 10000;

% Ar*x = u solution
x_ur = pinv(P*Ar + k_const*(hr*hr'), tol)*(P*u + k_const*hr);
x_ur_last = (1/n)*hr'*(u - Ar*x_ur);

x_u = [x_ur; x_ur_last];

% Ar*x = r solution
x_rr = pinv(P*Ar + k_const*(hr*hr'), tol)*(P*r+ k_const*hr);
x_rr_last = (1/n)*hr'*(r - Ar*x_rr);

x_r = [x_rr; x_rr_last];
    
% clock bias calculation
% Check if there is a real solution, if there is no real solution, the best
% result would be the local minimum or maximum (x = -b/(2a))

delta = (x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1];

if delta <= 0
    cdt = (x_u'*[r;0] + x_r'*[u;1])/(2*(1 + 2*x_r'*[r;0]));

elseif x_r'*[r;0] < 0
    cdt = (x_u'*[r;0] + x_r'*[u;1] + sqrt(delta))/(2*(1 + 2*x_r'*[r;0]));
else
    cdt = (x_u'*[r;0] + x_r'*[u;1] - sqrt(delta))/(2*(1 + 2*x_r'*[r;0]));
end
% calculate final x vector
x_r = x_ur - 2*cdt*x_rr;

XR = XS*x_r;
dtR = cdt/c;

[PDOP, HDOP, VDOP, cov_XYZ, cov_ENU] = DOP(XR, XS);
G = XS'*pinv(XS*XS'); %least-norm solution

XR_cov = 36*cov_XYZ;

end


