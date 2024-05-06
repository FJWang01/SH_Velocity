
function p = pressure_tf_pw(R, theta_pw, phi_pw, k)
% p = pressure_tf_pw(R, theta_pw, phi_pw, k)
% 
% This function calculates the pressure at point receiver R due to plane
% wave with incidence direction (theta_pw, phi_pw)
%
% Inputs:
% R - locations of point receivers in Cartesian coordinates
% theta_pw, phi_pw - plane wave incidence directions
% k - wavenumbers, must be a row vector
%
% Outputs:
% p - pressure at point receiver R due to plane wave with incidence
%     direction (theta_pw, phi_pw)
%   size(p) = [size(R, 1), numel(theta_pw), numel(k)]

%% Check the dimensions of inputs
if ~isequal(size(R, 2), 3)
    error('@@ pressure_tf_pw: R must have three columns');
else
    % do nothing
end

if ~isequal(size(theta_pw), size(phi_pw))
    error('@@ pressure_tf_pw: theta_pw and phi_pw must be of the same size');
else
    % do nothing
end

validateattributes(k, {'double'}, {'row'});
validateattributes(theta_pw, {'double'}, {'column'});
validateattributes(phi_pw, {'double'}, {'column'});

%% Calculate pressure at point receiver R
% Find unit vector in (theta_pw, phi_pw) direction
[x_pw, y_pw, z_pw] = sph2cart(phi_pw, pi/2 - theta_pw, 1);

[x_pw_mat, x_R_mat] = meshgrid(x_pw, R(:, 1));
[y_pw_mat, y_R_mat] = meshgrid(y_pw, R(:, 2));
[z_pw_mat, z_R_mat] = meshgrid(z_pw, R(:, 3));

dot_prod = x_pw_mat .* x_R_mat + y_pw_mat .* y_R_mat + z_pw_mat .* z_R_mat;

k_3D = reshape(k, 1, 1, numel(k));
dot_prod_3D = repmat(dot_prod, 1, 1, numel(k));

p = exp(1i * k_3D .* dot_prod_3D);
end