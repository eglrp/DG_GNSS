% Distance geometry software GPS receiver
%------------------------------------------------------------------------ 
% 
% This program provides GPS receiver functionality for post-processing
% data recorded in Receiver Independant Exchange (RINEX) files. The
% positioning algorithm is based on a distance geoemetry approach
% introduced by Shahrdad Tabib* and is part of an ongoing research project
% at the University of California, Davis.
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

% Loading RINEX navigation file
[Eph, iono, constellations] = load_RINEX_nav('brdm1660.13p', [], 0);
%[Eph, iono, constellations] = load_RINEX_nav('site247j.01n', [], 0);

% Loading RINEX observation file
[pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, time_ref, time, week, date, pos, interval, antoff, antmod, codeC1] = load_RINEX_obs('gmsd1660_cut.13o', []);
%[pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, time_ref, time, week, date, pos, interval, antoff, antmod, codeC1] = load_RINEX_obs('site247j.01o', []);

nEpochs = length(time);      


nSatTot = constellations.nEnabledSat;
err_iono = zeros(nSatTot,nEpochs);
err_tropo = zeros(nSatTot,nEpochs);            
dtR = zeros(length(time),1);
dtR_dot = zeros(length(time),1);
XR = zeros(3, length(time));
XS = zeros(4, 3, length(time));
time_interval = interval; %initialization of time interval

% Locate satellites using Eph (this will need receiver clock bias as input)
% For now assume dtR = 0 for first epoch.


for i = 1 : length(time)
    sat0 = find(pr1(:,i) ~= 0);

    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys, traveltime] = satellite_positions(time(i), pr1(:,i), sat0, Eph, [], [], err_tropo, err_iono, dtR(i,1));

    r = [pr1(sat0, i);0] ;
    u_bar_r = diag(pr1(sat0, i)*pr1(sat0, i)');
    u_bar = [u_bar_r; 1];
    h_r = ones(numel(pr1(sat0, i)), 1);

    P = eye(numel(pr1(sat0, i))) - 1/(numel(pr1(sat0, i)))*(h_r*h_r');

    A = zeros(numel(pr1(sat0, i))+1);
 
    for k = 1:numel(pr1(sat0, i))
        for j = 1:numel(pr1(sat0, i))
             A(k,j) = (XS(k,:) - XS(j,:))*(XS(k,:) - XS(j,:))';
        end
        A(k,numel(pr1(sat0, k))+1) = 1;
        A(numel(pr1(sat0, k))+1, k) = 1;
    end
 
    A_r = A(1:numel(pr1(sat0, i)), 1:numel(pr1(sat0, i)));
    % X matrix
    x_r_bar = pinv(P*A_r)*P*u_bar_r;
    XR_m(:,:,i) = XS'*x_r_bar;
    
end      




% Saving outputs
time_stamp = datestr(now, 'mmddyyHHMMSS');
mkdir(strcat('./data/', 'DG_', time_stamp));
pathname = strcat('./data/', 'DG_', time_stamp, '/');
save(strcat(pathname,'DG_XR_', time_stamp), 'XR');
save(strcat(pathname,'DG_time_', time_stamp), 'time');
save(strcat(pathname,'DG_dtR_', time_stamp), 'dtR');
save(strcat(pathname,'DG_pr1_', time_stamp), 'pr1');

