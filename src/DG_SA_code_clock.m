function [dtR, A, Ainv] = DG_SA_code_clock(XS, dtS, err_iono, err_tropo, pr_R, sat)
% Clock bias claculations using distance geometry for the case of exactly four satellites 
%   Tabib, S.: A coordinate free distance geometry approach to the GPS and kinematical
% computations for general localization. Ph.D. thesis, University of California, Davis
% (2007)

% input:
%   XS - satellite positions (each satellite on a row)
%   pr_R - code observations
%   sat - satellite index
   
% output:
%   dtR - receiver clock error (in m)
%   A - reference matrix
%   Ainv - inverse of reference matrix

if isempty(sat)
    k = find(~(pr_R(:)==0)); % find pseudorange indexes that are non-zero
    PR = pr_R(k);
else
    PR = pr_R(sat);
end

%PR = PR - goGNSS.V_LIGHT*dtS; % + err_iono + err_tropo;


r = [PR; 0];
u = diag(PR*PR');
u = [u; 1];
QR2 = (XS(1,:) - XS(2,:))*(XS(1,:) - XS(2,:))';
QS2 = (XS(1,:) - XS(3,:))*(XS(1,:) - XS(3,:))';
QT2 = (XS(1,:) - XS(4,:))*(XS(1,:) - XS(4,:))';
RS2 = (XS(2,:) - XS(3,:))*(XS(2,:) - XS(3,:))';
RT2 = (XS(2,:) - XS(4,:))*(XS(2,:) - XS(4,:))';
ST2 = (XS(3,:) - XS(4,:))*(XS(3,:) - XS(4,:))';

A = [0, QR2, QS2, QT2, 1;
     QR2, 0, RS2, RT2, 1;
     QS2, RS2, 0, ST2, 1;
     QT2, RT2, ST2, 0, 1;
     1, 1, 1, 1, 0];
 
% Consider scaling the matrix to avoid singularity

% supress MATLAB warning on singularity
id = 'MATLAB:nearlySingularMatrix';
warning('off',id); 
 
%Ainv = 1\A;  
Ainv = inv(A);
 
b = r'*Ainv*u;



dtR_dist = (b -sign(b)*sqrt((r'*Ainv*u)^2 -u'*Ainv*u*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
dtR = dtR_dist/goGNSS.V_LIGHT;

end