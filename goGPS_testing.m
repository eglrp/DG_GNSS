
% Readin RINEX observation file
[pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, time_ref, time, week, date, pos, interval, antoff, antmod, codeC1] = load_RINEX_obs('capo083v_cut.16o', []);

bad_epoch=NaN;
var_SPP = NaN(1,3);
is_bias=NaN(6,1);

SPP_threshold=4; %meters 
bad_sat=[];

XR0 = [];
if (isempty(XR0))
    XR0 = zeros(3,1);
end
Omegae_dot = goGNSS.OMEGAE_DOT_GPS;
POS = zeros(3, length(time));

for i = 1 : length(time)
    
    index  = find(pr1(:,i) ~= 0);
    % Get sat positions from SP3:
    [XS1, dtS] = sp3_linear('igu18893_18.sp3', index, date(i,:));
    
    for j = 1:length(index)
        [XS2] = earth_rotation_correction(pr1(index(j))/goGNSS.V_LIGHT, XS1(j,:), Omegae_dot);
        XS(j,:) = XS2';
    end
    
    nsat_avail = length(index);

    if (nsat_avail < 4) %if available observations are not enough, return empty variables
        XR   = [];
        dtR  = [];
        az   = [];
        el   = [];
        dist = [];
        sat  = [];
        err_tropo = [];
        err_iono  = [];
    end

    %iterative least-squares from XR0,i.e. given coordinates or the center of the Earth (i.e. [0; 0; 0])
    n_iter_max = 5;
    n_iter = 0;
    var_SPP(1) = Inf;
    while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
        [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP] = LS_SA_code2([], XS, pr1(index,i), zeros(nsat_avail,1), zeros(nsat_avail,1), zeros(nsat_avail,1), dtS, zeros(nsat_avail,1), zeros(nsat_avail,1), [], SPP_threshold);       
        %bad_sat(sat(bad_obs))=1;
        XR0 = XR;
        n_iter = n_iter + 1;
    end
    POS(:,i) = XR;
    CBIAS(i) = dtR;
end