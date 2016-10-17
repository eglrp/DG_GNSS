function [XR_T, dtR_T] = tikhonov(XS, r)
% Tikhonov Regularization:
n = numel(r);
sigma = 5*10^10;
c = 299792458;
Ar = zeros(n);
hr = ones(n,1);
P = eye(n) - 1/n*ones(n,n);

% Creating distance geometry reference matrix
for i = 1:n
    for j = 1:n
        tmp = (XS(:,i) - XS(:,j));
        Ar(i,j) = tmp'*tmp;
    end
end
    
% form u & r vectors
u = diag(r*r');

PAr = P*Ar;
GenInv = pinv(PAr'*PAr + (sigma*eye(n))'*(sigma*eye(n)))*PAr'; 
Xr_r = GenInv*P*r;
Xr_l = (1/n)*hr'*(u - Ar*Xr_r);
Xr = [Xr_r; Xr_l];

Xu_r = GenInv*P*u;
Xu_l = (1/n)*hr'*(u - Ar*Xu_r);
Xu = [Xu_r; Xu_l];

% Clock bias
mrr = Xr_r'*r; 
mru = Xr'*[u;1];
mur = Xu_r'*r;
muu = Xu'*[u;1];

if mrr < 0
    cdt_T = (mur + mru + sqrt((mur + mru)^2 - 2*(1 + 2*mrr)*muu))/(2*(1 + 2*mrr));
else
    cdt_T = (mur + mru - sqrt((mur + mru)^2 - 2*(1 + 2*mrr)*muu))/(2*(1 + 2*mrr));
end

X = Xu_r - 2*cdt_T*Xr_r;
XR_T = XS*X;
dtR_T = cdt_T/c;

end