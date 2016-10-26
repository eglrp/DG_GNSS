function [XR, dtR] = Geo_buildup3(XS, pr, indx, El)
% This function uses the localization method from the geometric build-up
% method for positioning.
% The receiver clock bias is still found using the distance geometry
% method.
% Chooses highest elevation satellite as base satellite.

[~, ~, dtR] = DG4sat(XS, pr, indx);
c = 299792458;
n = numel(pr);

% u corrected for clock bias
u = (pr - c*dtR*ones(n,1)).^2;

% Selecting base satellite:
[~,k] = max(El);
% k=3;
% Form the normal equation:
A = 2*(ones(n-1,1)*XS(:,k)' - XS(:,[1:k-1,k+1:n])');
b = XS(:,k)'*XS(:,k)*ones(n-1,1) - diag(XS(:,[1:k-1,k+1:n])'*XS(:,[1:k-1,k+1:n])) - u(k)*ones(n-1,1) + u([1:k-1,k+1:n]);

% Solving the LS using normal equation:
XR = pinv(A'*A)*A'*b;

end