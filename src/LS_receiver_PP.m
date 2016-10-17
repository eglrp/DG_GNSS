
% Least squares software GPS receiver
%------------------------------------------------------------------------ 
% 
% This program provides GPS receiver functionality for post-processing
% data recorded in Receiver Independant Exchange (RINEX) files. The
% positioning algorithm is based on a least square solution for the
% receiver position.
%
% Input:
% Navigation RINEX file
% Observation RINEX file
%
% Output:
% GPS position estimates in an ECEF reference frame (X, Y, Z)
% GPS position estimates in a geographic reference frame (lat., long., h)
%
%------------------------------------------------------------------------
%
% * Tabib, S., A coordinate free distance geometry approach to the GPS and 
% kinematical computations for general localization, Ph.D. thesis,
% University of California, Davis (2007).
%
%------------------------------------------------------------------------
%
% Copyright 2016, Hadi Tabatabaee, All rights reserved.
%
%------------------------------------------------------------------------

% Loading RINEX observation file
[pr1, ~, ~, ~, ~, ~, ~, ~, ~, time, ~, date, pos, interval, antoff, antmod, codeC1] = load_RINEX_obs('capo083v_cut.16o', []);

nEpochs = length(time);      


nSatTot = 31;
err_iono = zeros(nSatTot,nEpochs);
err_tropo = zeros(nSatTot,nEpochs);            
dtR = zeros(length(time),1);
dtR_dot = zeros(length(time),1);
XR = zeros(3, length(time));
%XS = zeros(7, 3, length(time));

% Locate satellites using Eph (this will need receiver clock bias as input)
% For now assume dtR = 0 for first epoch.

% Estimate drR for next epochs by using dtRdot.
for i = 1 : length(time)
%     sat0 = find(pr1(:,i) ~= 0);sat0 = find(pr1(:,i) ~= 0);
    sat0 = [7; 9; 23; 30; 16; 8; 3; 5];
    [XS, dtS] = sp3_lookup('igu18893_18.sp3', sat0, date(i,:));
    [XR(:,i), dtR(i), A, b] = leastSquare(sat0, XS, pr1(:,i));
    %[XR_geo(1,i), XR_geo(2,i), XR_geo(3,i)] = llh(XR(1,i), XR(2,i), XR(3,i));
end      

% Saving outputs
time_stamp = datestr(now, 'mmddyyHHMMSS');
mkdir(strcat('./data/', 'PP_LS_', time_stamp));
pathname = strcat('./data/', 'PP_LS_', time_stamp, '/');sat0 = find(pr1(:,i) ~= 0);
%save(strcat(pathname,'PP_LS_XS_', time_stamp), 'XS');
save(strcat(pathname,'PP_LS_XR_', time_stamp), 'XR');
save(strcat(pathname,'PP_LS_time_', time_stamp), 'time');
save(strcat(pathname,'PP_LS_dtR_', time_stamp), 'dtR');
%save(strcat(pathname,'PP_LS_dtS_', time_stamp), 'dtS');
save(strcat(pathname,'PP_LS_pr1_', time_stamp), 'pr1');

