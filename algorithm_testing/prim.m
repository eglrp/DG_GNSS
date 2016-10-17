function [XR, dtR] = prim(XS, r)

% Defining needed variables:
n = length(r);
Ar = zeros(n);
hr = ones(n,1);
c = 299792458; % speed of light.
P = eye(n) - 1/n*ones(n,n);

% Creating distance geometry reference matrix
for i = 1:n
    for j = 1:n
        tmp = XS(:,i) - XS(:,j);
        Ar(i,j) = tmp'*tmp;
    end
end

u = diag(r*r');
X = pinv(P*Ar)*P*u;

XR = XS*X;
dtR = 0;
end


