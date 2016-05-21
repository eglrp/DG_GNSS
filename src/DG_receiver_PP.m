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

% Loading RINEX observation file
[pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, time_ref, time, week, date, pos, interval, antoff, antmod, codeC1] = load_RINEX_obs('capo083v_cut.16o', []);

nEpochs = length(time);      

nSatTot = 31;
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
dtS_min = zeros(4, length(time));
dtS_max = zeros(4, length(time));
sat_max = [];
PR = zeros(4, length(time));

% Option 1: Estimate XS for next epochs by using satellite velocity.
for i = 1 : length(time)
        
    k_min = min_snr(pr1(:,i), snr1(:,i));
    k_max = max_snr(pr1(:,i), snr1(:,i));
    
    [XS_min(:,:,i), dtS_min(:,i)] = sp3_lookup('igu18894_18_cut.sp3', k_min, date(i,:)); 
    [dtR_min(i,1), A_min, Ainv_min] = DG_SA_code_clock(XS_min(:,:,i), dtS_min(:,i), err_iono, err_tropo, pr1(:,i), k_min);
    [XR_min(:,i)] = DG_SA_code(XS_min(:,:,i), pr1(:,i), dtR_min(i,1), dtS_min(:,i), err_iono, err_tropo, A_min, Ainv_min, k_min);
    
    [XS_max(:,:,i), dtS_max(:,i)] = sp3_lookup('igu18894_18_cut.sp3', k_max, date(i,:));
    [dtR_max(i,1), A_max, Ainv_max] = DG_SA_code_clock(XS_max(:,:,i), dtS_max(:,i), err_iono, err_tropo, pr1(:,i), k_max);
    [XR_max(:,i)] = DG_SA_code(XS_max(:,:,i), pr1(:,i), dtR_max(i,1), dtS_max(:,i), err_iono, err_tropo, A_max, Ainv_max, k_max);
    
    sat_max = [sat_max; k_max];
end      

% Saving outputs
time_stamp = datestr(now, 'mmddyyHHMMSS');
mkdir(strcat('./data/', 'DG_PP_', time_stamp));
pathname = strcat('./data/', 'DG_PP_', time_stamp, '/');
save(strcat(pathname,'PP_DG_XS_min_', time_stamp), 'XS_min');
save(strcat(pathname,'PP_DG_XS_max_', time_stamp), 'XS_max');
save(strcat(pathname,'PP_DG_XR_min_', time_stamp), 'XR_min');
save(strcat(pathname,'PP_DG_XR_max_', time_stamp), 'XR_max');
save(strcat(pathname,'PP_DG_time_', time_stamp), 'time');
save(strcat(pathname,'PP_DG_dtR_min_', time_stamp), 'dtR_min');
save(strcat(pathname,'PP_DG_dtR_max_', time_stamp), 'dtR_max');
save(strcat(pathname,'PP_DG_dtS_min_', time_stamp), 'dtS_min');
save(strcat(pathname,'PP_DG_dtS_max_', time_stamp), 'dtS_max');
save(strcat(pathname,'PP_DG_pr1_', time_stamp), 'pr1');


