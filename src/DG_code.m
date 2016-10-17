function [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, is_bias] = DG_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold)

% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, is_bias] = DG_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold);
%
% INPUT:
%   XR_approx     = receiver approximate position (X,Y,Z)
%   XS            = satellite position (X,Y,Z)
%   pr_R          = code observations
%   snr_R         = signal-to-noise ratio
%   elR           = satellite elevation (vector)
%   distR_approx  = approximate receiver-satellite distance (vector)
%   dtS           = satellite clock error (vector)
%   err_tropo_RS  = tropospheric error
%   err_iono_RS   = ionospheric error
%   sys           = array with different values for different systems
%   SPP_threshold = maximum RMS of code single point positioning to accept current epoch
%
% OUTPUT:
%   XR   = estimated position (X,Y,Z)
%   dtR  = receiver clock error (scalar)
%   cov_XR  = estimated position error covariance matrix
%   var_dtR = estimated clock error variance
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   cond_num = condition number on the eigenvalues of the N matrix
%   bad_obs = vector with ids of observations found as outlier
%   bad_epoch = 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   sigma02_hat = [a posteriori sigma (SPP sigma), v_hat'*(invQ)*v_hat), n-m] 
%   residuals_obs = vector with residuals of all input observation, computed from the final estimates
%   is_bias = inter-systems bias (vector with all possibile systems)
%
% DESCRIPTION:
%   Calculates receiver position and clock bias using the distance geometry
%   method.
%--------------------------------------------------------------------------
%
% Copyright (C) 2016, Hadi Tabatabaee. All rights reserved.
%
%--------------------------------------------------------------------------

is_bias = [1; 0; 0; 0; 0; 0];
v_light = goGNSS.V_LIGHT;
nSats = nnz(pr_R);
index = find(pr_R);
Ar = zeros(nSats);
hr = ones(nSats,1);
h = [hr;0];
m = 4; % number of unkowns

% reduced r and u vectors (not corrected for clock bias)
r = pr_R(index);
u = diag(r*r');

% Form the reference matrix using the satellite positions.
for i = 1:nSats
    for j = 1:nSats
        Ar(i,j) = (XS(i,:) - XS(j,:))*(XS(i,:) - XS(j,:))';
    end
end
A = [Ar, hr; hr', 0];

P = eye(nSats) - 1/(nSats)*hr*hr';

% Condition number of A (may change later)
cond_num = cond(A);

% Weighting constants (modify maybe...)
tol = 25;
k_const = 100;

% Ar*x = u solution
x_ur = pinv(P*Ar + k_const*hr*hr', tol)*(P*u + k_const*hr);
x_ur_last = (1/nSats)*hr'*(u - Ar*x_ur);

x_u = [x_ur; x_ur_last];
    
% Ar*x = r solution
x_rr = pinv(P*Ar + k_const*hr*hr', tol)*(P*r+ k_const*hr);
x_rr_last = (1/nSats)*hr'*(r - Ar*x_rr);

x_r = [x_rr; x_rr_last];

% clock bias calculation
cdt_plus = (x_u'*[r;0] + x_r'*[u;1] + sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));

cdt_minus = (x_u'*[r;0] + x_r'*[u;1] - sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));

if x_r'*[r;0] < 0
    dtR = cdt_plus;
else
    dtR = cdt_minus;
end

% calculate final x vector
x_r = x_ur - 2*dtR*x_rr;

% clock bias
dtR = dtR/v_light;

% position estimate
XR = XS'*x_r;

% covariance matrix computation (skipped for now)
%cov_XR = [];
var_dtR = [];

%sigma02_hat
sigma02_hat = NaN(1,3);
bad_obs = [];
bad_epoch = [];
residuals_obs = zeros(nSats, 1);

if (~any(distR_approx))
    %approximate receiver-satellite distance
    XR_mat = XR_approx(:,ones(nSats,1))';
    distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
end

%design matrix
Ad = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
    (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
    (XR_approx(3) - XS(:,3)) ./ distR_approx, ... %column for Z coordinate
     ones(nSats,1)];        %column for receiver clock delay (multiplied by c)
  
%DOP computation (copied from goGPS)
if (nargout > 4)
   cov_XYZ = (Ad'*Ad)^-1;
   cov_XYZ = cov_XYZ(1:3,1:3);
   cov_ENU = global2localCov(cov_XYZ, XR);
   
   cov_XR = cov_XYZ;

   PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
   HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
   VDOP = sqrt(cov_ENU(3,3));
end

end


