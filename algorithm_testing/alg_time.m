load('sat_final.mat');
load('signal.mat');

time_steps = numel(time);
pr = PR(k,:);
pr = pr(1:4,:);
nSats = size(pr,1);
ph_l1 = PH_L1(k,:);
ph_l2 = PH_L2(k,:);
snr = SNR(k,:);
Ar = zeros(nSats);
hr = ones(nSats,1);
h = [hr; 0];
iter = zeros(time_steps,1);
XS = XS(:,1:4,:);
for j = 1:100
    for i = 1:time_steps
        tic;
        [XR1(:,i), dtR1(i)] = DG_sol_eff(XS(:,:,i), pr(:,i));
        dg_time(i,j) = toc;

        tic;
        [XR2(:,i), dtR2(i)] = tikhonov(XS(:,:,i), pr(:,i));
        ti_time(i,j) = toc;

        tic;
        [XR3(:,i), iter(i), dtR3(i)] = leastSquare_eff(XS(:,:,i), pr(:,i));
        ls_time(i,j) = toc;
        
        tic;
        [XR4(:,i), dtR4(i)] = DG4sat(XS(:,:,i), pr(:,i));
        dg4_time(i,j) = toc;
    end
end
% 
% 
% legend1 = 'DG - with Tikhonov';
% legend2 = 'DG - without Tikhonov';
% legend3 = 'Least Squares';
% legend4 = 'Real Position';
% 
% % figure
% % plot(1:400, d1, 1:400, d2, 1:400, d3)
% 
% figure
% subplot(3,1,1); title('X - coordinate (ECEF)');
% ylabel('X (m)')
% hold on
% %plot(1:400, XR(1,:))
% plot(1:400, XR1(1,:) - pos(1))
% plot(1:400, XR2(1,:) - pos(1))
% plot(1:400, XR3(1,:) - pos(1))
% plot(1:400, 0*pos(1)*ones(400,1))
% legend(legend1, legend2, legend3, legend4)
% 
% subplot(3,1,2); title('Y - coordinate (ECEF)');
% ylabel('Y (m)')
% hold on
% %plot(1:400, XR(1,:))
% plot(1:400, XR1(2,:) - pos(2))
% plot(1:400, XR2(2,:) - pos(2))
% plot(1:400, XR3(2,:) - pos(2))
% plot(1:400, 0*pos(2)*ones(400,1))
% legend(legend1, legend2, legend3, legend4)
% 
% subplot(3,1,3); title('Z - coordinate (ECEF)');
% ylabel('Z (m)')
% xlabel('Time (s)')
% hold on
% %plot(1:400, XR(1,:))
% plot(1:400, XR1(3,:) - pos(3))
% plot(1:400, XR2(3,:) - pos(3))
% plot(1:400, XR3(3,:) - pos(3))
% plot(1:400, 0*pos(3)*ones(400,1))
% legend(legend1, legend2, legend3, legend4)
% 
% %% Convert from ECEF to geographic coordinates
% clear h
% 
% 
% [phi(:), lam(:), h(:)] = cart2geod(XR1(1,:), XR1(2,:), XR1(3,:));
% [phi2(:), lam2(:), h2(:)] = cart2geod(XR2(1,:), XR2(2,:), XR2(3,:));
% [phi_LS(:), lam_LS(:), h_LS(:)] = cart2geod(XR3(1,:), XR3(2,:), XR3(3,:));
% [phi_pos, lam_pos, h_pos] = cart2geod(pos(1), pos(2), pos(3));
% 
% % Convert to degrees
% phi_pos = (phi_pos./pi)*180;
% lam_pos = (lam_pos./pi)*180;
% 
% % Convert and normalize
% phi = (phi./pi)*180 - phi_pos;
% lam = (lam./pi)*180 - lam_pos;
% phi2 = (phi2./pi)*180 - phi_pos;
% lam2 = (lam2./pi)*180 - lam_pos;
% phi_LS = (phi_LS./pi)*180 - phi_pos;
% lam_LS = (lam_LS./pi)*180 - lam_pos;
% 
% 
% 
% figure
% title('Geographic coordinates')
% subplot(3,1,1); title('Longitude');
% xlim([0 400]);
% hold on
% plot(phi)
% plot(phi2)
% plot(phi_LS)
% plot(1:400, 0*phi_pos*ones(400,1))
% legend(legend1, legend2, legend3, legend4)
% hold off
% 
% subplot(3,1,2); title('Latitude');
% xlim([0 400]);
% hold on
% plot(lam)
% plot(lam2)
% plot(lam_LS)
% plot(1:400, 0*lam_pos*ones(400,1))
% legend(legend1, legend2, legend3, legend4)
% hold off
% 
% subplot(3,1,3); title('Elevation (m)');
% xlim([0 400]);
% hold on
% plot(h)
% plot(h2)
% plot(h_LS)
% plot(1:400, h_pos*ones(400,1))
% legend(legend1, legend2, legend3, legend4)
% hold off
% 
% %% ENU coordinates
% 
% XR_ENU1 = ecef2enu(XR1, pos, [phi_pos; lam_pos; h_pos]);
% XR_ENU2 = ecef2enu(XR2, pos, [phi_pos; lam_pos; h_pos]);
% XR_ENU3 = ecef2enu(XR3, pos, [phi_pos; lam_pos; h_pos]);
% 
% figure
% title('Geographic coordinates')
% subplot(3,1,1); title('East');
% xlim([0 400]);
% hold on
% plot(1:400,XR_ENU1(1,:))
% plot(1:400,XR_ENU2(1,:))
% plot(1:400,XR_ENU3(1,:))
% legend(legend1, legend2, legend3)
% hold off
% 
% subplot(3,1,2); title('North');
% xlim([0 400]);
% hold on
% plot(1:400,XR_ENU1(2,:))
% plot(1:400,XR_ENU2(2,:))
% plot(1:400,XR_ENU3(2,:))
% legend(legend1, legend2, legend3)
% hold off
% 
% subplot(3,1,3); title('Up (m)');
% xlim([0 400]);
% hold on
% plot(1:400,XR_ENU1(3,:))
% plot(1:400,XR_ENU2(3,:))
% plot(1:400,XR_ENU3(3,:))
% legend(legend1, legend2, legend3)
% hold off
% 
% 
% 
