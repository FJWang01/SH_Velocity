
% This file caculates and plots the velocity vectors due to a plane wave
% and a point source 

close all
clear
clc

%% Setup
f = 0.5e3;
c = 343.21; % speed of sound at 20 degrees Celcius
rho = 1.2042; % density of air at 20 degrees Celcius

k = 2*pi*f/c;

% Plane wave
theta_pw = pi/2;
phi_pw = deg2rad(160);

% point source
theta_s = theta_pw;
phi_s = phi_pw;
r_s = 0.7;
[x_s, y_s, z_s] = sph2cart(phi_s, pi/2 - theta_s, r_s);
Cart_s = [x_s, y_s, z_s];

R = 0.5; % radius of the reproduction region
N_theo = ceil(exp(1) * k * R/2); % SH truncation order
N = N_theo;
%% Calculate the PW pressure coefficients 
[n, m, ~] = getDegreeOrderPairs(N);
[n_vel, m_vel, ~] = getDegreeOrderPairs(N-1);
Ynm_pw = sphHarm_mat(n, m, theta_pw, phi_pw);
scaling_pw = 4*pi*1i.^(n);
alpha_nm_pw = scaling_pw .* conj(Ynm_pw);

%% Calculate the point source pressure coefficients
Ynm_s = sphHarm_mat(n, m, theta_s, phi_s);
hn_kr_s = sph_Hankel_2_adapted(n, k, r_s);
alpha_nm_ps = -1i*k*hn_kr_s .* conj(Ynm_s);

%% Calculate operator matrix for SH coefficients of velocity
[glbp2vel_x_mat, glbp2vel_y_mat, glbp2vel_z_mat]  = getLocalVelTranslationMats(N, rho, c);

%% Calculate SH coefficients of velocity
% PW
SH_coeffs_vel_x_pw = glbp2vel_x_mat * alpha_nm_pw;
SH_coeffs_vel_y_pw = glbp2vel_y_mat * alpha_nm_pw;
SH_coeffs_vel_z_pw = glbp2vel_z_mat * alpha_nm_pw;

% point source
SH_coeffs_vel_x_ps = glbp2vel_x_mat * alpha_nm_ps;
SH_coeffs_vel_y_ps = glbp2vel_y_mat * alpha_nm_ps;
SH_coeffs_vel_z_ps = glbp2vel_z_mat * alpha_nm_ps;

%% Find velocity in a region
% Setup the grid point
dx = R/10;
x = -R:dx:R;
y = -R:dx:R;
z = 0;

[x_mat, y_mat] = meshgrid(x, y);
z_mat = zeros(size(x_mat));
Cart_mat = [x_mat(:), y_mat(:), z_mat(:)];

[azim_mat, elev_mat, R_mat] = cart2sph(x_mat, y_mat, z_mat);
elev_mat = pi/2 - elev_mat;

% Calculate velocity in the region
jn_vel_kr_mat = sph_Bessel_1_adapted(n_vel, k, R_mat(:));
jn_vel_kr_mat = squeeze(jn_vel_kr_mat);

% Incorporate radial part to SH coefficients of velocity
% PW
SH_coeffs_vel_x_inc_rad_pw = SH_coeffs_vel_x_pw .* jn_vel_kr_mat;
SH_coeffs_vel_y_inc_rad_pw = SH_coeffs_vel_y_pw .* jn_vel_kr_mat;
SH_coeffs_vel_z_inc_rad_pw = SH_coeffs_vel_z_pw .* jn_vel_kr_mat;

% point source
SH_coeffs_vel_x_inc_rad_ps = SH_coeffs_vel_x_ps .* jn_vel_kr_mat;
SH_coeffs_vel_y_inc_rad_ps = SH_coeffs_vel_y_ps .* jn_vel_kr_mat;
SH_coeffs_vel_z_inc_rad_ps = SH_coeffs_vel_z_ps .* jn_vel_kr_mat;

% Calculate velocity
% PW
Ynm_vel_mat = sphHarm_mat(n_vel, m_vel, elev_mat(:), azim_mat(:));
Ynm_vel_mat = squeeze(Ynm_vel_mat);

vel_x_pw = sum(SH_coeffs_vel_x_inc_rad_pw .* Ynm_vel_mat);
vel_y_pw = sum(SH_coeffs_vel_y_inc_rad_pw .* Ynm_vel_mat);
vel_z_pw = sum(SH_coeffs_vel_z_inc_rad_pw .* Ynm_vel_mat);

vel_x_mat_pw = reshape(vel_x_pw, size(x_mat));
vel_y_mat_pw = reshape(vel_y_pw, size(x_mat));
vel_z_mat_pw = reshape(vel_z_pw, size(x_mat));

% point source
vel_x_ps = sum(SH_coeffs_vel_x_inc_rad_ps .* Ynm_vel_mat);
vel_y_ps = sum(SH_coeffs_vel_y_inc_rad_ps .* Ynm_vel_mat);
vel_z_ps = sum(SH_coeffs_vel_z_inc_rad_ps .* Ynm_vel_mat);

vel_x_mat_ps = reshape(vel_x_ps, size(x_mat));
vel_y_mat_ps = reshape(vel_y_ps, size(x_mat));
vel_z_mat_ps = reshape(vel_z_ps, size(x_mat));

%% Calculate theoretical velocity
% PW
[vel_x_pw_theo, vel_y_pw_theo, vel_z_pw_theo] = velocity_tf_pw(Cart_mat, theta_pw, phi_pw, k, rho, c);

vel_x_pw_theo_mat = reshape(vel_x_pw_theo, size(x_mat));
vel_y_pw_theo_mat = reshape(vel_y_pw_theo, size(x_mat));
vel_z_pw_theo_mat = reshape(vel_z_pw_theo, size(x_mat));

% Point source
[vel_x_ps_theo, vel_y_ps_theo, vel_z_ps_theo] = velocity_tf_point_source(Cart_mat, Cart_s, k, rho, c);
vel_x_ps_theo_mat = reshape(vel_x_ps_theo, size(x_mat));
vel_y_ps_theo_mat = reshape(vel_y_ps_theo, size(x_mat));
vel_z_ps_theo_mat = reshape(vel_z_ps_theo, size(x_mat));

%% Plot velocity 
figure;
tiledlayout(1, 2, 'TileSpacing','tight', 'Padding','compact');
nexttile;
quiver(x, y, real(vel_x_mat_pw), real(vel_y_mat_pw));
hold on;
quiver(x, y, real(vel_x_pw_theo_mat), real(vel_y_pw_theo_mat), '--');
hold on;
viscircles([0, 0], R);
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('Real part of velocity PW - SH');

nexttile;
quiver(x, y, real(vel_x_mat_ps), real(vel_y_mat_ps));
hold on;
quiver(x, y, real(vel_x_ps_theo_mat), real(vel_y_ps_theo_mat), '--');
hold on;
viscircles([0, 0], R);
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('Real part of velocity PS - SH');

% nexttile;
% quiver(x, y, real(vel_x_pw_theo_mat), real(vel_y_pw_theo_mat));
% hold on;
% viscircles([0, 0], R);
% xlabel('x (m)');
% ylabel('y (m)');
% axis equal;
% xlim([-R, R]);
% ylim([-R, R]);
% title('Real part of velocity PW - Theoretical');

% nexttile;
% quiver(x, y, real(vel_x_ps_theo_mat), real(vel_y_ps_theo_mat));
% hold on;
% viscircles([0, 0], R);
% xlabel('x (m)');
% ylabel('y (m)');
% axis equal;
% xlim([-R, R]);
% ylim([-R, R]);
% title('Real part of velocity PS - Theoretical');

figure;
tiledlayout(1, 2, 'TileSpacing','tight', 'Padding','compact');
nexttile;
quiver(x, y, imag(vel_x_mat_pw), imag(vel_y_mat_pw));
hold on;
quiver(x, y, imag(vel_x_pw_theo_mat), imag(vel_y_pw_theo_mat), '--');
hold on;
viscircles([0, 0], R);
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('Imag part of velocity PW - SH');

nexttile;
quiver(x, y, imag(vel_x_mat_ps), imag(vel_y_mat_ps));
hold on;
quiver(x, y, imag(vel_x_ps_theo_mat), imag(vel_y_ps_theo_mat), '--');
hold on;
viscircles([0, 0], R);
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('Imag part of velocity PS - SH');

% nexttile;
% quiver(x, y, imag(vel_x_pw_theo_mat), imag(vel_y_pw_theo_mat));
% hold on;
% viscircles([0, 0], R);
% xlabel('x (m)');
% ylabel('y (m)');
% axis equal;
% xlim([-R, R]);
% ylim([-R, R]);
% title('Imag part of velocity PW - Theoretical');

% nexttile;
% quiver(x, y, imag(vel_x_ps_theo_mat), imag(vel_y_ps_theo_mat));
% hold on;
% viscircles([0, 0], R);
% xlabel('x (m)');
% ylabel('y (m)');
% axis equal;
% xlim([-R, R]);
% ylim([-R, R]);
% title('Imag part of velocity PS - Theoretical');

%% Plot for SPL paper
figure;
tiledlayout(1, 2, 'TileSpacing','tight', 'Padding','compact');
nexttile;
quiver(x, y, real(vel_x_mat_pw), real(vel_y_mat_pw));
hold on;
viscircles([0, 0], R);
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('Real part of velocity PW - SH');

nexttile;
quiver(x, y, real(vel_x_mat_ps), real(vel_y_mat_ps));
hold on;
viscircles([0, 0], R);
xlabel('x (m)');
ylabel('y (m)');
axis equal;
xlim([-R, R]);
ylim([-R, R]);
title('Real part of velocity PS - SH');