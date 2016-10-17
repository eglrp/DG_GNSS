function [XR, dtR] = distG(XS, pr)

[~, dtR] = DG4sat(XS(:,1:4,:), pr(1:4));
c = 299792458;
n = numel(pr);
% u corrected for clock bias
u = (pr - c*dtR*ones(n,1)).^2;

for i = 1:n
    for j = 1:n
        tmp = XS(:,i) - XS(:,j);
        Ar(i,j) = tmp'*tmp;
    end
end
P = eye(n) - 1/n*ones(n,n);

% Solution:
x = pinv(P*Ar, 25)*P*u;

XR = XS*x;