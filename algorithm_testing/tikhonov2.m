function [XR_T, dtR_T] = tikhonov2(XS, pr)
% Tikhonov Regularization:
n = numel(pr);
sigma = 2.5*10^13;
c = 299792458;
Ar = zeros(n);
hr = ones(n,1);
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

Xr_r = pinv((P*Ar)'*(P*Ar) + (sigma*eye(n))'*(sigma*eye(n)))*(P*Ar)'*P*r;
Xr_l = (1/n)*hr'*(u - Ar*Xr_r);
Xr = [Xr_r; Xr_l];

Xu_r = pinv((P*Ar)'*(P*Ar) + (sigma*eye(10))'*(sigma*eye(10)))*(P*Ar)'*P*u;
Xu_l = (1/n)*hr'*(u - Ar*Xu_r);
Xu = [Xu_r; Xu_l];

% Clock bias
if Xr'*[r;0] < 0
    cdt_T = (Xu'*[r;0] + Xr'*[u;1] + sqrt((Xu'*[r;0] + Xr'*[u;1])^2 - 2*(1 + 2*Xr'*[r;0])*Xu'*[u;1]))/(2*(1 + 2*Xr'*[r;0]));
else
    cdt_T = (Xu'*[r;0] + Xr'*[u;1] - sqrt((Xu'*[r;0] + Xr'*[u;1])^2 - 2*(1 + 2*Xr'*[r;0])*Xu'*[u;1]))/(2*(1 + 2*Xr'*[r;0]));
end

X = Xu_r - 2*cdt_T*Xr_r;
XR_T = XS*X;
dtR_T = cdt_T/c;

end