% compare all combinations of four satellites.
load('sat_final.mat');
load('signal.mat');

% All possible choices of satellite combinations with just four satellites.
indx = pick(1:10, 7, '');
m = size(indx,1);

for i = 1:size(indx,1)
    [PDOP(i), HDOP(i), VDOP(i)] = DOP(pos, XS(:,indx(i,:),1));
end

pr = PR(k,:);
nSats = size(pr,1);
ph_l1 = PH_L1(k,:);
ph_l2 = PH_L2(k,:);
snr = SNR(k,:);
Ar = zeros(nSats);
hr = ones(nSats,1);
h = [hr; 0];

for i = 1:size(indx,1)
    [XR(:,i), ~,dtR(i)] = DG_sol(XS(:,indx(i,:),1), pr(indx(i,:),1));
    [XR_LS(:,i), ~, dtR_LS(i)] = LS_dtR(XS(:,indx(i,:),1), pr(indx(i,:),1), [], []);
end



%% Cartesian coordinates plotting
figure
title('ECEF')
subplot(4,1,1)
xlim([0 m]);
hold on
plot(1:m, XR(1,:)-pos(1));
plot(1:m, XR_LS(1,:)-pos(1));
plot(1:m, 0*pos(1)*ones(m,1));
hold off

subplot(4,1,2)
xlim([0 m]);
hold on
plot(1:m, XR(2,:)-pos(2));
plot(1:m, XR_LS(2,:)-pos(2));
plot(1:m, 0*pos(2)*ones(m,1));
hold off

subplot(4,1,3)
xlim([0 m]);
hold on
plot(1:m, XR(3,:) - pos(3));
plot(1:m, XR_LS(3,:)-pos(3));
plot(1:m, 0*pos(3)*ones(m,1));
hold off

subplot(4,1,4)
xlim([0 m]);
hold on
plot(PDOP)
plot(HDOP)
plot(VDOP)
hold off

%% Geodetic Coordinates plotting
phi = zeros(1,m);
phi_LS = zeros(1,m);
lam = zeros(1,m);
lam_LS = zeros(1,m);
h = zeros(1,m);
h_LS = zeros(1,m);


[phi(:), lam(:), h(:)] = ECEF2GPS(XR(1,:), XR(2,:), XR(3,:));
[phi_LS(:), lam_LS(:), h_LS(:)] = ECEF2GPS(XR_LS(1,:), XR_LS(2,:), XR_LS(3,:));
[phi_pos, lam_pos, h_pos] = ECEF2GPS(pos(1), pos(2), pos(3));

% Convert to degrees
phi_pos = (phi_pos./pi)*180;
lam_pos = (lam_pos./pi)*180;

% Convert and normalize
phi = (phi./pi)*180 - phi_pos;
lam = (lam./pi)*180 - lam_pos;
phi_LS = (phi_LS./pi)*180 - phi_pos;
lam_LS = (lam_LS./pi)*180 - lam_pos;



figure
title('Geographic coordinates')
subplot(4,1,1); title('Longitude');
xlim([0 m]);
hold on
plot(phi)
plot(phi_LS)
plot(1:m, 0*phi_pos*ones(m,1))
hold off

subplot(4,1,2); title('Latitude');
xlim([0 m]);
hold on
plot(lam)
plot(lam_LS)
plot(1:m, 0*lam_pos*ones(m,1))
hold off

subplot(4,1,3); title('Elevation (m)');
xlim([0 m]);
hold on
plot(h)
plot(h_LS)
plot(1:m, h_pos*ones(m,1))
hold off

subplot(4,1,4); title('Dilution of Precision (DOP)');
xlim([0 m]);
hold on
plot(PDOP)
plot(HDOP)
plot(VDOP)
hold off

%% Planar coordinates plotting

EAST = zeros(1,m);
EAST_LS = zeros(1,m);
NORTH = zeros(1,m);
NORTH_LS = zeros(1,m);
ht = zeros(1,m);
ht_LS = zeros(1,m);


[EAST(:), NORTH(:), ht(:), utm_zone] = cart2plan(XR(1,:), XR(2,:), XR(3,:));
[EAST_LS(:), NORTH_LS(:), ht_LS(:), utm_zone_LS] = cart2plan(XR_LS(1,:), XR_LS(2,:), XR_LS(3,:));
[EAST_pos, NORTH_pos, ht_pos, utm_zone_pos] = cart2plan(pos(1), pos(2), pos(3));

figure
title('Planimetric coordinates')
subplot(4,1,1); title('East coordinates (m)');
hold on
plot(EAST(:))
plot(EAST_LS(:))
plot(1:m, EAST_pos*ones(m,1))
hold off

subplot(4,1,2); title('North coordinates (m)');
hold on
plot(NORTH(:))
plot(NORTH_LS(:))
plot(1:m, NORTH_pos*ones(m,1))
hold off

subplot(4,1,3); title('Elevation (m)');
hold on
plot(ht(:))
plot(ht_LS(:))
plot(1:m, ht_pos*ones(m,1))
hold off

subplot(4,1,4); title('Dilution of precision (DOP)');
hold on
plot(PDOP)
plot(HDOP)
plot(VDOP)
hold off

