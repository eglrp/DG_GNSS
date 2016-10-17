function [XR, dtR] = DG4sat(XS, pr)
% This function performs the distance geometery solution for only four
% satellites.
% Defining needed variables:
n = length(pr);
Ar = zeros(n);
hr = (10^14)*ones(n,1);
c = 299792458; % speed of light.
for i = 1:n
    for j = 1:n
        tmp = XS(:,i) - XS(:,j);
        Ar(i,j) = tmp'*tmp;
    end
end
A = [Ar, hr;hr',0];    
% form u & r vectors
u = [diag(pr*pr');(10^14)];
r = [pr;0];

Ainv = pinv(A);
b = r'*Ainv*u;
b1 = r'*Ainv*r;

dtR_dist = (b -sign(b)*sqrt(b^2 -u'*Ainv*u*(0.5 + b1)))/(2*(0.5 + b1));
dtR = dtR_dist/c;

PR = pr - dtR_dist*ones(4,1);
u = [diag(PR*PR'); (10^14)];
XS = [XS, zeros(3,1)];

XR = XS*Ainv*u;
% condA = cond(A);
% condAr = cond(Ar);
end