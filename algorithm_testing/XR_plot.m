function XR_plot(XR1, XR2, XR3, pos, le)
% Compares different position estimates from three methods with a known
% actual position.
% ---------------------------------------------------------------------
% Inputs:
%       pos = actual ECEF coordinates of receiver (if known)
%       XR1 = ECEF coordinates achieved from first method
%             XR1 is a 3xm matrix where m is the number of epochs
%       le = structure containing legend names under field 'legends'
% ---------------------------------------------------------------------
% Copyright 2016, Hadi Tabatabaee. All rights reserved.
% ---------------------------------------------------------------------

if isempty(le)
    legend1 = 'distance geometry';
    legend2 = 'DG - KF';
    legend3 = 'Least Squares';
    legend4 = 'Real Position';
else
    legend1 = le(1).legends;
    legend2 = le(2).legends;
    legend3 = le(3).legends;
    legend4 = le(4).legends;
end

n = size(XR1, 2);


%% ECEF PLOTS

maxX = max(max([XR1(1,:)'; XR2(1,:)'; XR3(1,:)']) - pos(1), 0);
minX = min(min([XR1(1,:)'; XR2(1,:)'; XR3(1,:)']) - pos(1), 0);
ranX = maxX - minX;
maxY = max(max([XR1(2,:)'; XR2(2,:)'; XR3(2,:)']) - pos(2), 0);
minY = min(min([XR1(2,:)'; XR2(2,:)'; XR3(2,:)']) - pos(2), 0);
ranY = maxY - minY;
maxZ = max(max([XR1(3,:)'; XR2(3,:)'; XR3(3,:)']) - pos(3), 0);
minZ = min(min([XR1(3,:)'; XR2(3,:)'; XR3(3,:)']) - pos(3), 0);
ranZ = maxZ - minZ;
p = 0.05;

figure
subplot(3,1,1); title('X - coordinate (ECEF)');
ylabel('X (m)')
xlim([0 n+1])
ylim([minX-p*ranX maxX+p*ranX])
hold on
plot(1:n, XR1(1,:) - pos(1), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(1:n, XR2(1,:) - pos(1), 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(1:n, XR3(1,:) - pos(1), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*pos(1)*ones(n,1), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)

subplot(3,1,2); title('Y - coordinate (ECEF)');
ylabel('Y (m)')
xlim([0 n+1])
ylim([minY-p*ranY maxY+p*ranY])
hold on
plot(1:n, XR1(2,:) - pos(2), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(1:n, XR2(2,:) - pos(2), 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(1:n, XR3(2,:) - pos(2), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*pos(2)*ones(n,1), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)

subplot(3,1,3); title('Z - coordinate (ECEF)');
ylabel('Z (m)')
xlabel('Time (s)')
xlim([0 n+1])
ylim([minZ-p*ranZ maxZ+p*ranZ])
hold on
plot(1:n, XR1(3,:) - pos(3), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(1:n, XR2(3,:) - pos(3), 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(1:n, XR3(3,:) - pos(3), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*pos(3)*ones(n,1), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)

%% LLH PLOTS
clear h


[phi1(:), lam1(:), h1(:)] = cart2geod(XR1(1,:), XR1(2,:), XR1(3,:));
[phi2(:), lam2(:), h2(:)] = cart2geod(XR2(1,:), XR2(2,:), XR2(3,:));
[phi3(:), lam3(:), h3(:)] = cart2geod(XR3(1,:), XR3(2,:), XR3(3,:));
[phi_pos, lam_pos, h_pos] = cart2geod(pos(1), pos(2), pos(3));

% Convert to degrees
phi_pos = (phi_pos./pi)*180;
lam_pos = (lam_pos./pi)*180;

% Convert and normalize
phi1 = (phi1./pi)*180 - phi_pos;
lam1 = (lam1./pi)*180 - lam_pos;
phi2 = (phi2./pi)*180 - phi_pos;
lam2 = (lam2./pi)*180 - lam_pos;
phi3 = (phi3./pi)*180 - phi_pos;
lam3 = (lam3./pi)*180 - lam_pos;

maxphi = max(max([phi1'; phi2'; phi3']), 0);
minphi = min(min([phi1'; phi2'; phi3']), 0);
ranphi = maxphi - minphi;
maxlam = max(max([lam1'; lam2'; lam3']), 0);
minlam = min(min([lam1'; lam2'; lam3']), 0);
ranlam = maxlam - minlam;
maxh = max(max([h1'; h2'; h3']) - h_pos, 0);
minh = min(min([h1'; h2'; h3']) - h_pos, 0);
ranh = maxh - minh;

figure
title('Geographic coordinates')
subplot(3,1,1); title('Longitude');
xlim([0 n+1]);
ylim([minphi-p*ranphi maxphi+p*ranphi])
hold on
plot(phi1, 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255] , 'LineWidth', 1.0)
plot(phi2, 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(phi3, '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*phi_pos*ones(n,1), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)
hold off

subplot(3,1,2); title('Latitude');
xlim([0 n+1]);
ylim([minlam-p*ranlam maxlam+p*ranlam])
hold on
plot(lam1, 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(lam2, 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(lam3, '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*lam_pos*ones(n,1), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)
hold off

subplot(3,1,3); title('Elevation (m)');
xlim([0 n+1]);
ylim([minh-p*ranh maxh+p*ranh])
hold on
plot(h1 - h_pos*ones(1,n), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(h2 - h_pos*ones(1,n), 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(h3 - h_pos*ones(1,n), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*h_pos*ones(n,1), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)
hold off

%% 3D POSITIONING ERROR
for i = 1:n
    norm1(i) = norm(XR1(:,i) - pos);
    norm2(i) = norm(XR2(:,i) - pos);
    norm3(i) = norm(XR3(:,i) - pos);
end

figure
title('3D positioning error (m)');
ylabel('Error (m)')
xlabel('Time (s)')
hold on

plot(1:n, norm1, 'LineWidth', 2.0)
plot(1:n, norm2, 'LineWidth', 2.0)
plot(1:n, norm3, 'LineWidth', 2.0)

legend(legend1, legend2, legend3)


%% ENU PLOTS

XR_ENU1 = ecef2enu(XR1, pos, [phi_pos; lam_pos; h_pos]);
XR_ENU2 = ecef2enu(XR2, pos, [phi_pos; lam_pos; h_pos]);
XR_ENU3 = ecef2enu(XR3, pos, [phi_pos; lam_pos; h_pos]);

maxE = max(max([XR_ENU1(1,:)'; XR_ENU2(1,:)'; XR_ENU3(1,:)']), 0);
minE = min(min([XR_ENU1(1,:)'; XR_ENU2(1,:)'; XR_ENU3(1,:)']), 0);
ranE = maxX - minX;
maxN = max(max([XR_ENU1(2,:)'; XR_ENU2(2,:)'; XR_ENU3(2,:)']), 0);
minN = min(min([XR_ENU1(2,:)'; XR_ENU2(2,:)'; XR_ENU3(2,:)']), 0);
ranN = maxY - minY;
maxU = max(max([XR_ENU1(3,:)'; XR_ENU2(3,:)'; XR_ENU3(3,:)']), 0);
minU = min(min([XR_ENU1(3,:)'; XR_ENU2(3,:)'; XR_ENU3(3,:)']), 0);
ranU = maxZ - minZ;




figure
title('Geographic coordinates')
subplot(3,1,1); title('East (m)');
ylabel('East (m)');
xlim([0 n+1]);
ylim([minE-p*ranE maxE+p*ranE])
hold on
plot(1:n, XR_ENU1(1,:), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(1:n, XR_ENU2(1,:), 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(1:n, XR_ENU3(1,:), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*ones(1,n), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)
hold off

subplot(3,1,2); title('North (m)');
ylabel('North (m)');
xlim([0 n+1]);
ylim([minN-p*ranN maxN+p*ranN])
hold on
plot(1:n, XR_ENU1(2,:), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(1:n, XR_ENU2(2,:), 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(1:n, XR_ENU3(2,:), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*ones(1,n), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)
hold off

subplot(3,1,3); title('Up (m)');
xlim([0 n+1]);
ylim([minU-p*ranU maxU+p*ranU])
ylabel('Elevation (m)')
hold on
plot(1:n, XR_ENU1(3,:), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(1:n, XR_ENU2(3,:), 'color', [216/255 82/255 24/255], 'LineWidth', 2.0)
plot(1:n ,XR_ENU3(3,:), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(1:n, 0*ones(1,n), '--', 'color', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)
hold off

%% PLANAR 2D PLOT
figure
title('Geographic Coordinates - 2D')
xlabel('East (m)');
ylabel('North (m)');
hold on
plot(XR_ENU1(1,:), XR_ENU1(2,:), 'gs', 'MarkerSize',5,'MarkerEdgeColor',[0,113/255,188/255], 'LineWidth', 1.0)
plot(XR_ENU2(1,:), XR_ENU2(2,:), 'rx' ,'MarkerEdgeColor', [216/255 82/255 24/255])
plot(XR_ENU3(1,:), XR_ENU3(2,:), '^', 'MarkerSize',4,'MarkerEdgeColor',[236/255, 176/255, 31/255], 'MarkerFaceColor',[1, 215/255, 0])
plot(0,0, '*','MarkerEdgeColor', [125/255 45/255 141/255], 'LineWidth', 1.5)
legend(legend1, legend2, legend3, legend4)
hold off




end
