

A = P*Ar*10^-5;
rk = 4; % rank of matrix
n = size(A,1); % no. of satellites
B = [A,[eye(rk); zeros(n-rk,rk)]];

Br = rref(B);

A11_inv = Br(1:rk, n+1:n+rk);
Q_m = Br(1:rk, rk+1:n);
P_m = -Br(rk+1:n, n+1:n+rk); % probably always zero

% Now the generalized inverse:
Ginv = [eye(rk); Q_m']*inv(eye(rk) + Q_m*Q_m')*A11_inv*inv(eye(rk) + P_m'*P_m)*[eye(rk) P_m'];

% Now x solutions:
u_xr = Ginv*(P*u);
u_xl = (1/n)*hr'*(u - Ar*u_xr);
u_xl2 = (1/n)*hr'*(u - P*Ar*u_xr);

u_x = [u_xr; u_xl];
u_x2 = [u_xr; u_xl2];

r_xr = Ginv*(P*r);
r_xl = (1/n)*hr'*(r - Ar*r_xr);
r_xl2 = (1/n)*hr'*(r - P*Ar*r_xr);

r_x = [r_xr; r_xl];
r_x2 = [r_xr; r_xl2];


if r_x'*[r;0] < 0
    cdt_G = (u_x'*[r;0] + r_x'*[u;1] + sqrt((u_x'*[r;0] + r_x'*[u;1])^2 - 2*(1 + 2*r_x'*[r;0])*u_x'*[u;1]))/(2*(1 + 2*r_x'*[r;0]));
else
    cdt_G = (u_x'*[r;0] + r_x'*[u;1] - sqrt((u_x'*[r;0] + r_x'*[u;1])^2 - 2*(1 + 2*r_x'*[r;0])*u_x'*[u;1]))/(2*(1 + 2*r_x'*[r;0]));
end


if r_x2'*[r;0] < 0
    cdt_G2 = (u_x2'*[r;0] + r_x2'*[u;1] + sqrt((u_x2'*[r;0] + r_x2'*[u;1])^2 - 2*(1 + 2*r_x2'*[r;0])*u_x2'*[u;1]))/(2*(1 + 2*r_x2'*[r;0]));
else
    cdt_G2 = (u_x2'*[r;0] + r_x2'*[u;1] - sqrt((u_x2'*[r;0] + r_x2'*[u;1])^2 - 2*(1 + 2*r_x2'*[r;0])*u_x2'*[u;1]))/(2*(1 + 2*r_x2'*[r;0]));
end

% calculate final x vector
xr_G = u_xr - 2*cdt_G*r_xr;
xr_G2 = u_xr - 2*cdt_G2*r_xr;

XR_G = XS(:,:,400)*xr_G;
XR_G2 = XS(:,:,400)*xr_G2;
