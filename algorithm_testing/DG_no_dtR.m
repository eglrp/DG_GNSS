function XR = DG_no_dtR(XS, r, dtR)
c = 299792458;
u = (r-c*dtR).^2;
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

x = pinv(P*Ar,25)*(P*u);
XR = XS*x;

end
