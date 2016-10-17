function [XR, XR_cov, dtR] = DG_sol_genInv(XS, pr)

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
x_ur = genInv(P*Ar + k_const*(hr*hr'))*(P*u + k_const*hr);
x_ur_last = (1/n)*hr'*(u - Ar*x_ur);

x_u = [x_ur; x_ur_last];

% Ar*x = r solution
x_rr = genInv(P*Ar + k_const*(hr*hr'))*(P*r+ k_const*hr);
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
% 
% [PDOP, HDOP, VDOP, cov_XYZ, cov_ENU] = DOP(XR, XS);
% G = XS'*pinv(XS*XS'); %least-norm solution
% 
% XR_cov = 36*cov_XYZ;
XR_cov = 0;
end


