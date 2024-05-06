
% This file contains new SPL plots
close all
clear
clc

%% Setup of speaker system
% % NHK 22 channel system
% r_spk = [0.97; 1.55*ones(8, 1); 1.21*ones(10, 1); 1.34*ones(3, 1)]; % in metres
% theta_deg_spk = [0; 51.3*ones(8, 1); 90*ones(10, 1); 115*ones(3, 1)]; % in degrees
% phi_deg_spk = [0;0;45;90;135;180;225;270;315;0;22.5;45;90;135;180;225;270;315;337.5;0;45;315]; % in degrees
% 
% theta_spk = deg2rad(theta_deg_spk);
% phi_spk = deg2rad(phi_deg_spk);
% 
% [x_spk, y_spk, z_spk] = sph2cart(phi_spk, pi/2 - theta_spk, r_spk);
% 
% Cart_spk = [x_spk, y_spk, z_spk];

% % 5 channel system
% r_spk = 1.21*ones(5, 1);
% theta_deg_spk = [90;90;90;90;90];
% phi_deg_spk = [0;45;135;225;315];
% 
% theta_spk = deg2rad(theta_deg_spk);
% phi_spk = deg2rad(phi_deg_spk);
% 
% [x_spk, y_spk, z_spk] = sph2cart(phi_spk, pi/2 - theta_spk, r_spk);
% 
% Cart_spk = [x_spk, y_spk, z_spk];

% 8 channel system
r_spk = 1*ones(8, 1);
theta_deg_spk = [58.3;58.3;58.3;90;90;121.7;121.7;148.3];
phi_deg_spk = [288;216;72;18;126;324;180;72];

theta_spk = deg2rad(theta_deg_spk);
phi_spk = deg2rad(phi_deg_spk);

[x_spk, y_spk, z_spk] = sph2cart(phi_spk, pi/2 - theta_spk, r_spk);

Cart_spk = [x_spk, y_spk, z_spk];


%% Setup 
f = 200:100:8e3;
c = 343.21; % speed of sound at 20 degrees Celcius
rho = 1.2042; % density of air at 20 degrees Celcius

k = 2*pi*f/c;

theta_pw = pi/2;
phi_pw = deg2rad(160);

R = 0.5; % radius of the reproduction region, in metres
N = 4; % SH truncation order for Eigenmike
% N = ceil(exp(1)*max(k)*R/2);
R_sw = 0.15; % sweet spot radius

%%  Plot setup
% Contour lines
phi_contour = 0:0.01:2*pi;
x_contour = cos(phi_contour);
y_contour = sin(phi_contour);

z1 = cosd(58.3);
r1 = sind(58.3);

z2 = cosd(90);
r2 = sind(90);

z3 = cosd(121.7);
r3 = sind(121.7);

z4 = cosd(148.3);
r4 = sind(148.3);

[x_u, y_u, z_u] = sphere(256);

% Speaker labels
spk_labels = {'(58.3, 288)', '(58.3, 216)', '(58.3, 72)', '(90, 18)', '(90, 126)',...
    '(121.7, 324)', '(121.7, 180)', '(148.3, 72)'};

figure;
surf(x_u, y_u, z_u, 'FaceColor', 'y', 'EdgeColor', 'none', 'FaceAlpha',0.1); % spherical speaker array aperture
hold on;
surf(x_u*R, y_u*R, z_u*R, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha',0.1); % spherical listening region
hold on;
surf(x_u*R_sw, y_u*R_sw, z_u*R_sw, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha',0.1); % sweet spot
hold on;
plot3(x_contour * r1, y_contour * r1, z1*ones(size(x_contour)), 'k'); % first contour
hold on; 
plot3(x_contour * r2, y_contour * r2, z2*ones(size(x_contour)), 'k'); % second contour
hold on;
plot3(x_contour * r3, y_contour * r3, z3*ones(size(x_contour)), 'k'); % third contour
hold on;
plot3(x_contour * r4, y_contour * r4, z4*ones(size(x_contour)), 'k'); % fourth contour
hold on;
scatter3(x_spk, y_spk, z_spk, 'x', 'k', 'LineWidth',10); % loudspeakers
text(x_spk+0.05, y_spk+0.05, z_spk+0.05, spk_labels); % loudspeaker coordinate labels
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
axis equal;

clear phi_contour x_contour y_contour z1 r1 z2 r2 z3 r3 z4 r4 x_u y_u z_u spk_labels;

%% Calculate the desired velocity coefficients
[n, m, ~] = getDegreeOrderPairs(N);
[n_vel, m_vel, ~] = getDegreeOrderPairs(N-1);
Ynm_pw = sphHarm_mat(n, m, theta_pw, phi_pw);
scaling_pw = 4*pi*1i.^(n);

% size(alpha_nm_pw) = [numel(n), 1]
alpha_nm_pw = scaling_pw .* conj(Ynm_pw);

%% Calculate operator matrix for SH coefficients of velocity
% size(glbp2vel_x_mat) = [numel(n_vel), numel(n)]
[glbp2vel_x_mat, glbp2vel_y_mat, glbp2vel_z_mat]  = getLocalVelTranslationMats(N, rho, c);

%% Calculate SH coefficients of velocity
% size(SH_coeffs_vel_x_d) = [numel(n_vel), 1]
SH_coeffs_vel_x_d = glbp2vel_x_mat * alpha_nm_pw;
SH_coeffs_vel_y_d = glbp2vel_y_mat * alpha_nm_pw;
SH_coeffs_vel_z_d = glbp2vel_z_mat * alpha_nm_pw;

% size(SH_coeffs_vel_d_mat) = [3*numel(n_vel), 1]
SH_coeffs_vel_d_mat = [SH_coeffs_vel_x_d; SH_coeffs_vel_y_d; SH_coeffs_vel_z_d];

%% Calculate the SH coefficients of velocity due to each speaker
% size(Ynm_spk) = [numel(n), 1, numel(theta_spk)]
Ynm_spk = sphHarm_mat(n, m, theta_spk, phi_spk);

% size(hankel_kr_spk) = [numel(n), numel(k), numel(theta_spk)]
hankel_kr_spk = sph_Hankel_2_adapted(n, k, r_spk);

% size(alpha_nm_spk) = [numel(n), numel(k), numel(theta_spk)]
alpha_nm_spk = -1i*k .*hankel_kr_spk .* conj(Ynm_spk);

% size(SH_coeffs_vel_x_spk) = [numel(n_vel), numel(k), numel(theta_spk)]
SH_coeffs_vel_x_spk = pagemtimes(glbp2vel_x_mat, alpha_nm_spk);
SH_coeffs_vel_y_spk = pagemtimes(glbp2vel_y_mat, alpha_nm_spk);
SH_coeffs_vel_z_spk = pagemtimes(glbp2vel_z_mat, alpha_nm_spk);

%% Calculate the weight of each speaker using velocity matching
w_spk_VM = zeros(numel(theta_spk), numel(k));
cond_num_VM = zeros(1, numel(k));
for freq_idx = 1:numel(k)
    % size(SH_ceoffs_vel_spk_mat) = [3*numel(n_vel), numel(theta_spk)]
    SH_coeffs_vel_spk_mat = [squeeze(SH_coeffs_vel_x_spk(:, freq_idx, :));
                             squeeze(SH_coeffs_vel_y_spk(:, freq_idx, :));
                             squeeze(SH_coeffs_vel_z_spk(:, freq_idx, :))];
    cond_num_VM(freq_idx) = cond(SH_coeffs_vel_spk_mat);
    
    w_spk_VM(:, freq_idx) = pinv(SH_coeffs_vel_spk_mat) * SH_coeffs_vel_d_mat;
end
%% Calculate the weight of each speaker using pressure matching
w_spk_PM = zeros(numel(theta_spk), numel(k));
cond_num_PM = zeros(1, numel(k));
for freq_idx = 1:numel(k)
    % size(matching_mat_PM) = [numel(n), numel(theta_spk)]
    matching_mat_PM = squeeze(alpha_nm_spk(:, freq_idx, :));
    
    w_spk_PM(:, freq_idx) = pinv(matching_mat_PM) * alpha_nm_pw;
    
    cond_num_PM(freq_idx) = cond(matching_mat_PM);
end
%% Reconstruction grid
% Setup the grid point
dx = R/10;
x = -R:dx:R;
y = -R:dx:R;
z = -R:dx:R;

[x_mat, y_mat, z_mat] = meshgrid(x, y, z);

Cart_recon_col = [x_mat(:), y_mat(:), z_mat(:)];
[~, ~, r_recon] = cart2sph(x_mat(:), y_mat(:), z_mat(:));

%% Calculate reproduced pressure
% Unit pressure due to each loudspeaker
p_unit = pressure_tf_point_source(Cart_recon_col, Cart_spk, k);

% Add loudspeaker weights
p_recon_VM_col = zeros(numel(x_mat), numel(k));
p_recon_PM_col = zeros(numel(x_mat), numel(k));

for freq_idx = 1:numel(k)
    % velocity matching
    p_recon_VM_col(:, freq_idx) = p_unit(:, :, freq_idx) * w_spk_VM(:, freq_idx);

    % pressure matching
    p_recon_PM_col(:, freq_idx) = p_unit(:, :, freq_idx) * w_spk_PM(:, freq_idx);
end

% Rearrage reconstructed pressure
p_recon_VM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
p_recon_PM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
for freq_idx = 1:numel(k)
    % velocity matching
    p_recon_VM(:, :, :, freq_idx) = reshape(p_recon_VM_col(:, freq_idx), size(x_mat));

    % pressure matching
    p_recon_PM(:, :, :, freq_idx) = reshape(p_recon_PM_col(:, freq_idx), size(x_mat));
end

%% Calculate desired pressure
p_des_col = pressure_tf_pw(Cart_recon_col, theta_pw, phi_pw, k);

p_des = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
for freq_idx = 1:numel(k)
    p_des(:, :, :, freq_idx) = reshape(p_des_col(:, :, freq_idx), size(x_mat));
end

%% Calculate reproduced velocity
% Unit velocity due to each loudspeaker
[vel_x_unit, vel_y_unit, vel_z_unit] = velocity_tf_point_source(Cart_recon_col, Cart_spk, k, rho, c);

% Add loudspeaker weights
vel_x_recon_VM_col = zeros(numel(x_mat), numel(k));
vel_x_recon_PM_col = zeros(numel(x_mat), numel(k));

vel_y_recon_VM_col = zeros(numel(x_mat), numel(k));
vel_y_recon_PM_col = zeros(numel(x_mat), numel(k));

vel_z_recon_VM_col = zeros(numel(x_mat), numel(k));
vel_z_recon_PM_col = zeros(numel(x_mat), numel(k));

for freq_idx = 1:numel(k)
    % velocity matching
    vel_x_recon_VM_col(:, freq_idx) = vel_x_unit(:, :, freq_idx) * w_spk_VM(:, freq_idx);
    vel_y_recon_VM_col(:, freq_idx) = vel_y_unit(:, :, freq_idx) * w_spk_VM(:, freq_idx);
    vel_z_recon_VM_col(:, freq_idx) = vel_z_unit(:, :, freq_idx) * w_spk_VM(:, freq_idx);

    % pressure matching
    vel_x_recon_PM_col(:, freq_idx) = vel_x_unit(:, :, freq_idx) * w_spk_PM(:, freq_idx);
    vel_y_recon_PM_col(:, freq_idx) = vel_y_unit(:, :, freq_idx) * w_spk_PM(:, freq_idx);
    vel_z_recon_PM_col(:, freq_idx) = vel_z_unit(:, :, freq_idx) * w_spk_PM(:, freq_idx);
end

% Rearrage reconstructed velocity
vel_x_recon_VM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
vel_x_recon_PM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));

vel_y_recon_VM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
vel_y_recon_PM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));

vel_z_recon_VM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
vel_z_recon_PM = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));

for freq_idx = 1:numel(k)
    % velocity matching
    vel_x_recon_VM(:, :, :, freq_idx) = reshape(vel_x_recon_VM_col(:, freq_idx), size(x_mat));
    vel_y_recon_VM(:, :, :, freq_idx) = reshape(vel_y_recon_VM_col(:, freq_idx), size(x_mat));
    vel_z_recon_VM(:, :, :, freq_idx) = reshape(vel_z_recon_VM_col(:, freq_idx), size(x_mat));

    % pressure matching
    vel_x_recon_PM(:, :, :, freq_idx) = reshape(vel_x_recon_PM_col(:, freq_idx), size(x_mat));
    vel_y_recon_PM(:, :, :, freq_idx) = reshape(vel_y_recon_PM_col(:, freq_idx), size(x_mat));
    vel_z_recon_PM(:, :, :, freq_idx) = reshape(vel_z_recon_PM_col(:, freq_idx), size(x_mat));
end

%% Calculate desired velocity
[vel_x_des_col, vel_y_des_col, vel_z_des_col] = velocity_tf_pw(Cart_recon_col, theta_pw, phi_pw, k, rho, c);

vel_x_des = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
vel_y_des = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));
vel_z_des = zeros(size(x_mat, 1), size(x_mat, 2), size(x_mat, 3), numel(k));

for freq_idx = 1:numel(k)
    vel_x_des(:, :, :, freq_idx) = reshape(vel_x_des_col(:, :, freq_idx), size(x_mat));
    vel_y_des(:, :, :, freq_idx) = reshape(vel_y_des_col(:, :, freq_idx), size(x_mat));
    vel_z_des(:, :, :, freq_idx) = reshape(vel_z_des_col(:, :, freq_idx), size(x_mat));
end

%% Calculate angular error between the desired velocity and the reproduced velocity
ang_err_real_VM = zeros(numel(x_mat), numel(k));
ang_err_imag_VM = zeros(numel(x_mat), numel(k));
DOT_real_VM = zeros(numel(x_mat), numel(k));
DOT_imag_VM = zeros(numel(x_mat), numel(k));

ang_err_real_PM = zeros(numel(x_mat), numel(k));
ang_err_imag_PM = zeros(numel(x_mat), numel(k));
DOT_real_PM = zeros(numel(x_mat), numel(k));
DOT_imag_PM = zeros(numel(x_mat), numel(k));

for freq_idx = 1:numel(k)
    vel_des_col = [vel_x_des_col(:, freq_idx), vel_y_des_col(:, freq_idx), vel_z_des_col(:, freq_idx)];
    
    vel_recon_VM_col = [vel_x_recon_VM_col(:, freq_idx), vel_y_recon_VM_col(:, freq_idx), vel_z_recon_VM_col(:, freq_idx)];
    
    vel_recon_PM_col = [vel_x_recon_PM_col(:, freq_idx), vel_y_recon_PM_col(:, freq_idx), vel_z_recon_PM_col(:, freq_idx)];
    
    [ang_err_real_VM(:, freq_idx), ang_err_imag_VM(:, freq_idx), DOT_real_VM(:, freq_idx), DOT_imag_VM(:, freq_idx)] = ...
        get_ang_error(vel_des_col, vel_recon_VM_col);

    [ang_err_real_PM(:, freq_idx), ang_err_imag_PM(:, freq_idx), DOT_real_PM(:, freq_idx), DOT_imag_PM(:, freq_idx)] = ...
        get_ang_error(vel_des_col, vel_recon_PM_col);
end

ang_err_real_VM_sel_whole = ang_err_real_VM(r_recon<=R, :);
ang_err_imag_VM_sel_whole = ang_err_imag_VM(r_recon<=R, :);
ang_err_real_PM_sel_whole = ang_err_real_PM(r_recon<=R, :);
ang_err_imag_PM_sel_whole = ang_err_imag_PM(r_recon<=R, :);

ave_ang_err_real_VM_whole = mean(ang_err_real_VM_sel_whole, 'omitnan');
ave_ang_err_imag_VM_whole = mean(ang_err_imag_VM_sel_whole, 'omitnan');
ave_ang_err_real_PM_whole = mean(ang_err_real_PM_sel_whole, 'omitnan');
ave_ang_err_imag_PM_whole = mean(ang_err_imag_PM_sel_whole, 'omitnan');


ang_err_real_VM_sel_centre = ang_err_real_VM(r_recon<=R_sw, :);
ang_err_imag_VM_sel_centre = ang_err_imag_VM(r_recon<=R_sw, :);
ang_err_real_PM_sel_centre = ang_err_real_PM(r_recon<=R_sw, :);
ang_err_imag_PM_sel_centre = ang_err_imag_PM(r_recon<=R_sw, :);

ave_ang_err_real_VM_centre = mean(ang_err_real_VM_sel_centre, 'omitnan');
ave_ang_err_imag_VM_centre = mean(ang_err_imag_VM_sel_centre, 'omitnan');
ave_ang_err_real_PM_centre = mean(ang_err_real_PM_sel_centre, 'omitnan');
ave_ang_err_imag_PM_centre = mean(ang_err_imag_PM_sel_centre, 'omitnan');
%% Plot reconstructed pressure and velocity on 2D plane
equator_idx = round(numel(z)/2);
freq_idx = 6;
freq_plot = f(freq_idx);

figure;
tiledlayout(2, 2, 'TileSpacing','tight', 'padding', 'tight');
nexttile;
imagesc(x, y, real(p_recon_VM(:, :, equator_idx, freq_idx)));
hold on;
viscircles([0, 0], R);
hold on;
viscircles([0, 0], R_sw, 'Color', 'c');
ax = gca;
ax.YDir = 'normal';
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('VM');
clim([-1, 1]);

nexttile;
imagesc(x, y, real(p_recon_PM(:, :, equator_idx, freq_idx)));
hold on;
viscircles([0, 0], R);
hold on;
viscircles([0, 0], R_sw, 'Color', 'c');
ax = gca;
ax.YDir = 'normal';
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('PM');
cB = colorbar;
cB.Label.String = 'Amplitude';
clim([-1, 1]);

% figure;
% tiledlayout(1, 2, 'TileSpacing','tight', 'padding', 'tight');
nexttile;
quiver(x, y, real(vel_x_recon_VM(:, :, equator_idx, freq_idx)), real(vel_y_recon_VM(:, :, equator_idx, freq_idx)));
hold on;
viscircles([0, 0], R);
hold on;
viscircles([0, 0], R_sw, 'Color', 'c');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('VM');

nexttile;
quiver(x, y, real(vel_x_recon_PM(:, :, equator_idx, freq_idx)), real(vel_y_recon_PM(:, :, equator_idx, freq_idx)));
hold on;
viscircles([0, 0], R);
hold on;
viscircles([0, 0], R_sw, 'Color', 'c');
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('PM');

%% Plot desired pressure and velocity
figure;
imagesc(x, y, real(p_des(:, :, equator_idx, freq_idx)));
hold on;
viscircles([0, 0], R);
ax = gca;
ax.YDir = 'normal';
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
cB = colorbar;
cB.Label.String = 'Amplitude';
%caxis([-1, 1]);

figure;
quiver(x, y, real(vel_x_des(:, :, equator_idx, freq_idx)), real(vel_y_des(:, :, equator_idx, freq_idx)));
hold on;
viscircles([0, 0], R);
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);

%% Plot condition number
figure;
semilogx(f, cond_num_VM, '-*');
hold on;
semilogx(f, cond_num_PM, '-s');
xlabel('Frequency (Hz)');
ylabel('Condition number');
% ylim([1.3, 1.45]);
xlim([200, 8e3]);
legend('VM', 'PM');
grid on;

%% Plot velocity error
figure;
semilogx(f, mag2db(ave_ang_err_real_VM_whole), '-*');
hold on;
semilogx(f, mag2db(ave_ang_err_real_PM_whole), '-s');

hold on;
semilogx(f, mag2db(ave_ang_err_real_VM_centre), '-x');
hold on;
semilogx(f, mag2db(ave_ang_err_real_PM_centre), '-o');
xlim([200, 8e3]);
%ylim([-25, -4]);
ylabel('Direction error (dB)');
xlabel('Frequency (Hz)');
xticks([200, 500, 1000, 2000, 4000, 8000]);
grid on;
legend('VM, R = 0.5', 'PM, R = 0.5', 'VM, R = 0.15', 'PM, R = 0.15')
