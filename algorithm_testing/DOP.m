function [PDOP, HDOP, VDOP, cov_XYZ, cov_ENU] = DOP(XR, XS)
% Dilution of Precision and covariance matrix calculations for GPS
% estimates.
% Uses some functions from GoGPS MATLAB.
%--------------------------------------------------------------------------
% Copyright 2016, Hadi Tabatabaee, all rights reserved.
%--------------------------------------------------------------------------

n = size(XS,2);
% Note thatt XS needs to be an nx3 matrix (n = no. of sats)
for i = 1:n
    distR_approx = sqrt((XR - XS(:,i))'*((XR - XS(:,i))));
% Design matrix:
    A(i,:) = [(XR(1) - XS(1,i)) ./ distR_approx, ... %column for X coordinate
        (XR(2) - XS(2,i)) ./ distR_approx, ... %column for Y coordinate
        (XR(3) - XS(3,i)) ./ distR_approx, ... %column for Z coordinate
        1];
end
 % DOP computation:
   cov_XYZ = pinv(A'*A);
   cov_XYZ = cov_XYZ(1:3,1:3);
   cov_ENU = global2localCov(cov_XYZ, XR);
   
   PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
   HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
   VDOP = sqrt(cov_ENU(3,3));
   
end
   
