function [XR] = DG_SA_code(XS, pr_R, dtR, dtS, err_tropo_RS, err_iono_RS, A, Ainv, sat)

% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, is_bias] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold);
%
% INPUT:
%   XS            = satellite position (X,Y,Z)
%   pr_R          = code observations
%   dtS           = satellite clock error (vector)
%   err_tropo_RS  = tropospheric error
%   err_iono_RS   = ionospheric error
%   A             = satellite reference matrix
%   Ainv          = inverse of reference matrix
%   sat           = satellite index
%
% OUTPUT:
%   XR   = estimated position (X,Y,Z)

% DESCRIPTION:
%   Absolute positioning by means of using the algebraic distance geometry solution on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           
%  This program was written by Hadi Tabatabaee. It uses functions and
%  variables from the open-source goGPS project.
%
%----------------------------------------------------------------------------------------------

if isempty(sat)
    k = find(~(pr_R(:)==0)); % find pseudorange indexes that are non-zero
    PR_RAW = pr_R(k);
else
    PR_RAW = pr_R(sat);
end

PR = PR_RAW - dtR*goGNSS.V_LIGHT*[1;1;1;1];

% Add an if statement to calculate A and Ainv if they have not already been
% calculated

u = diag(PR*PR');
u = [u; 1];

XS = [XS; [0,0,0]]';


XR = XS*Ainv*u;


end
