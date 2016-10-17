    % This program tests the distance geometry algorithm for the over defined
% case.
% It uses 400 data points for a static receiver locked onto 10 satellites.
% Position estimates are using the pseudorange broadcasted on the L1
% frequency band.
% No attempts have been made to remove ionospheric and tropospheric delays.
% The statellite positions were calculated using a piecewise smooth cubic
% spline on the satellie positions provided by the precise ephemerides.

% THIS ONE GIVES GOOD CLOCK BIAS ESTIMATE !!!

% COMPARING WITH THE TWO LESS ACCURATE METHODS:
% METHOD 3 IS SIGNIFICANTLY LESS ACCURATE
% METHOD 2 GIVES ALMOST IDENTICAL RESULTS TWO THE MAIN METHOD (UP TO (
% DECIMAL PLACES --> ALL SIGNIFICANT DIGITS)
%-------------------------------------------------------------------------
% Hadi Tabatabaee, all rights reserved.
%-------------------------------------------------------------------------

load('sat_final.mat');
load('signal.mat');

time_steps = numel(time);

pr = PR(k,:);
nSats = size(pr,1);
ph_l1 = PH_L1(k,:);
ph_l2 = PH_L2(k,:);
snr = SNR(k,:);
Ar = zeros(nSats);
hr = ones(nSats,1);
h = [hr; 0];

for t = 1:time_steps
    % form reference matrix
    Ar = zeros(nSats);
 
    for i = 1:nSats
        for j = 1:nSats
            Ar(i,j) = (XS(:,i,t) - XS(:,j,t))'*(XS(:,i,t) - XS(:,j,t));
        end
    end

    A = [Ar, hr; hr', 0];
    
    P = eye(nSats) - 1/(nSats)*hr*hr';
    
    % form u & r vectors
    u = diag(pr(1:10,t)*pr(1:10,t)');
    r = pr(1:10,t);
    
    % weighting constants
    tol = 25;
    k_const = 100;
    
    % Ar*x = u solution
    x_ur = pinv(P*Ar + k_const*hr*hr', tol)*(P*u + k_const*hr);
    x_ur_last = (1/nSats)*hr'*(u - Ar*x_ur);

    x_u = [x_ur; x_ur_last];
    
%     % Less accurate estimates:
%     x_ur2 = pinv(P*Ar, tol)*(P*u); 
%     x_ur3 = inv(P*Ar)*P*u;
%     
%     x_ur_last2 = (1/nSats)*hr'*(u - Ar*x_ur2);
%     x_ur_last3 = (1/nSats)*hr'*(u - Ar*x_ur3);
%     
%     x_u2 = [x_ur2; x_ur_last2];
%     x_u3 = [x_ur3; x_ur_last3];
%     
    % Ar*x = r solution
    x_rr = pinv(P*Ar + k_const*hr*hr', tol)*(P*r+ k_const*hr);
    x_rr_last = (1/nSats)*hr'*(r - Ar*x_rr);

    x_r = [x_rr; x_rr_last];
    
%     % Less accurate estimates:
%     x_rr2 = pinv(P*Ar, tol)*(P*r);
%     x_rr3 = inv(P*Ar)*P*r;
%     
%     x_rr_last2 = (1/nSats)*hr'*(r - Ar*x_rr2);
%     x_rr_last3 = (1/nSats)*hr'*(r - Ar*x_rr3);
%     
%     x_r2 = [x_rr2; x_rr_last2];
%     x_r3 = [x_rr3; x_rr_last3];
%     
    % clock bias calculation
    cdt_plus = (x_u'*[r;0] + x_r'*[u;1] + sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));

    cdt_minus = (x_u'*[r;0] + x_r'*[u;1] - sqrt((x_u'*[r;0] + x_r'*[u;1])^2 - 2*(1 + 2*x_r'*[r;0])*x_u'*[u;1]))/(2*(1 + 2*x_r'*[r;0]));

    if x_r'*[r;0] < 0
        cdt = cdt_plus;
    else
        cdt = cdt_minus;
    end
    
%     % Less accurate estimates:
%     % clock bias calculation
%     if x_r2'*[r;0] < 0
%         cdt2 = (x_u2'*[r;0] + x_r2'*[u;1] + sqrt((x_u2'*[r;0] + x_r2'*[u;1])^2 - 2*(1 + 2*x_r2'*[r;0])*x_u2'*[u;1]))/(2*(1 + 2*x_r2'*[r;0]));
%     else
%         cdt2 = (x_u2'*[r;0] + x_r2'*[u;1] - sqrt((x_u2'*[r;0] + x_r2'*[u;1])^2 - 2*(1 + 2*x_r2'*[r;0])*x_u2'*[u;1]))/(2*(1 + 2*x_r2'*[r;0]));
%     end
%     
%     if x_r3'*[r;0] < 0
%         cdt3 = (x_u3'*[r;0] + x_r3'*[u;1] + sqrt((x_u3'*[r;0] + x_r3'*[u;1])^2 - 2*(1 + 2*x_r3'*[r;0])*x_u3'*[u;1]))/(2*(1 + 2*x_r3'*[r;0]));
%     else
%         cdt3 = (x_u3'*[r;0] + x_r3'*[u;1] - sqrt((x_u3'*[r;0] + x_r3'*[u;1])^2 - 2*(1 + 2*x_r3'*[r;0])*x_u3'*[u;1]))/(2*(1 + 2*x_r3'*[r;0]));
%     end
%     
    
    % calculate final x vector
    x_r = x_ur - 2*cdt*x_rr;

%     % Less accurate estimate:
%     x_r2 = x_ur2 - 2*cdt2*x_rr2;
%     x_r3 = x_ur3 - 2*cdt3*x_rr3;

    % position estimate
    XR(:,t) = XS(:,:,t)*x_r;
    
%     % Less accurate esimates
%     XR2(:,t) = XS(:,:,t)*x_r2;
%     XR3(:,t) = XS(:,:,t)*x_r3;
end