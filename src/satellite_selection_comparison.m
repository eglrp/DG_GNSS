% Distance geometry software GPS receiver -- overdefined
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

% Loading RINEX observation file
[pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, time_ref, time, week, date, pos, interval, antoff, antmod, codeC1] = load_RINEX_obs('gmsd1660_cut.13o', []);

nEpochs = length(time);      

nSatTot = constellations.nEnabledSat;
err_iono = zeros(nSatTot,nEpochs);
err_tropo = zeros(nSatTot,nEpochs);            
dtR_min = zeros(length(time),1);
dtR_max = zeros(length(time),1);
dtR_max2 = zeros(length(time),1);
dtR_dot_min = zeros(length(time),1);
dtR_dot_max = zeros(length(time),1);
XR_min = zeros(3, length(time));
XR_max = zeros(3, length(time));
XR_max2 = zeros(3, length(time));
XS_min = zeros(4, 3, length(time));
XS_max = zeros(4, 3, length(time));
XS_max2 = zeros(4, 3, length(time));




% Locate satellites using Eph (this will need receiver clock bias as input)
% For now assume dtR = 0 for first epoch.

% Estimate drR for next epochs by using dtRdot.
for i = 1 : length(time)
    
    k_min = min_snr(pr1(:,i), snr1(:,i));
    k_max = max_snr(pr1(:,i), snr1(:,i));
    
    % CALCULATIONS FOR K_MIN
    % clock bias estimation for satellite positioning purposes.
%     if i > 2
%         dtR_dot_min(i-1) = (dtR_min(i-1,1) - dtR_min(i-2,1))/(time(i-1) - time(i-2));
%     end
%     if i > 1
%         dtR_min(i,1) = dtR_min(i-1,1) + (time(i) - time(i-1))*dtR_dot_min(i-1);
%     end

%     if i > 1
%         dtR_min(i, 1) = dtR_min(i-1,1);
%     end
    [XS_min(:,:,i), dtS_min(:,i), XS_tx_min, VS_tx_min, time_tx_min, no_eph_min, sys_min] = satellite_positions(time(i), pr1(:,i), k_min, Eph, [], [], err_iono, err_tropo, dtR_min(i,1));
    [dtR_min(i,1), A_min, Ainv_min] = DG_SA_code_clock(XS_min(:,:,i), dtS_min(:,i), err_iono, err_tropo, pr1(:,i), k_min);
    [XR_min(:,i)] = DG_SA_code(XS_min(:,:,i), pr1(:,i), dtR_min(i,1), dtS_min(:,i), err_iono, err_tropo, A_min, Ainv_min, k_min);
    %[XR_geo_min(1,i), XR_geo_min(2,i), XR_geo_min(3,i)] = llh(XR_min(1,i), XR_min(2,i), XR_min(3,i));
    
    % CALCULATIONS FOR K_MAX
    % clock bias estimation for satellite positioning purposes.
%     if i > 2
%         dtR_dot_max(i-1) = (dtR_max(i-1,1) - dtR_max(i-2,1))/(time(i-1) - time(i-2));
%     end
%     if i > 1
%         dtR_max(i,1) = dtR_max(i-1,1) + (time(i) - time(i-1))*dtR_dot_max(i-1);
%     end
%     if i > 1
%         dtR_min(i, 1) = dtR_min(i-1,1);
%     end
    [XS_max(:,:,i), dtS_max(:,i), XS_tx_max, VS_tx_max, time_tx_max, no_eph_max, sys_max] = satellite_positions(time(i), pr1(:,i), k_max, Eph, [], [], err_iono, err_tropo, dtR_max(i,1));
    [dtR_max(i,1), A_max, Ainv_max] = DG_SA_code_clock(XS_max(:,:,i), dtS_max(:,i), err_iono, err_tropo, pr1(:,i), k_max);
    [XR_max(:,i)] = DG_SA_code(XS_max(:,:,i), pr1(:,i), dtR_max(i,1), dtS_max(:,i), err_iono, err_tropo, A_max, Ainv_max, k_max);
    %[XR_geo_max(1,i), XR_geo_max(2,i), XR_geo_max(3,i)] = llh(XR_max(1,i), XR_max(2,i), XR_max(3,i));
    
    % LS CALCULATIONS FOR K_MAX
    [XS_max2(:,:,i), dtS_max2(:,i), XS_tx_max2, VS_tx_max2, time_tx_max2, no_eph_max2, sys_max2] = satellite_positions(time(i), pr1(:,i), k_max, Eph, [], [], err_iono, err_tropo, dtR_max2(i,1));
    [XR_max2(:,i), dtR_max2(i), A, b] = leastSquare(k_max, XS_max2(:,:,i), pr1(:,i));
    
end      

% Saving outputs
time_stamp = datestr(now, 'mmddyyHHMMSS');
mkdir(strcat('./data/', 'DG_compare_', time_stamp));
pathname = strcat('./data/', 'DG_compare_', time_stamp, '/');
save(strcat(pathname,'DG_XS_min_', time_stamp), 'XS_min');
save(strcat(pathname,'DG_XS_max_', time_stamp), 'XS_max');
save(strcat(pathname,'DG_XR_min_', time_stamp), 'XR_min');
save(strcat(pathname,'DG_XR_max_', time_stamp), 'XR_max');
%save(strcat(pathname,'DG_XR_geo_min_', time_stamp), 'XR_geo_min');
%save(strcat(pathname,'DG_XR_geo_max_', time_stamp), 'XR_geo_max');
save(strcat(pathname,'DG_time_', time_stamp), 'time');
save(strcat(pathname,'DG_dtR_min_', time_stamp), 'dtR_min');
save(strcat(pathname,'DG_dtR_max_', time_stamp), 'dtR_max');
save(strcat(pathname,'DG_dtS_min_', time_stamp), 'dtS_min');
save(strcat(pathname,'DG_dtS_max_', time_stamp), 'dtS_max');
save(strcat(pathname,'DG_pr1_', time_stamp), 'pr1');

% plot dtR comparison
figure, hold on
plot(dtR_min, 'blue')
plot(dtR_max, 'red')
hold off

% subplot XYZ comparison
figure, hold on
subplot(3,1,1);
plot(XR_min(1,:), 'blue')
plot(XR_max(1,:), 'red')

subplot(3,1,2);
plot(XR_min(2,:), 'blue')
plot(XR_max(2,:), 'red')

subplot(3,1,3);
plot(XR_min(3,:), 'blue')
plot(XR_max(3,:), 'red')
hold off