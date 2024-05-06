
function [v_x, v_y, v_z] = velocity_tf_point_source(R, Rs, k, rho, c)
% [v_x, v_y, v_z] = velocity_tf_point_source(R, Rs, k, rho, c)
% 
% This function calculates the velocity vector at receiver R due to point 
% source at Rs.
%
% Inputs:
% Rs - locations of point sources in Cartesian coordinates
% R - locations of point receivers in Cartesian coordinates
% k - wavenumbers, must be a row vector
% rho - density of air, in kg/m^3, scalar
% c - speed of sound in metres per second, scalar
%
% Outputs:
% v_x, v_y, v_z - velocity vector at R due to point source at Rs
%   size(v_x) = size(v_y) = size(v_z) = [size(R, 1), size(Rs, 1), numel(k)]

%% Check the dimensions of inputs
if ~isequal(size(R, 2), size(Rs, 2), 3)
    error('@@ velocity_tf_point_source: R and Rs must have three columns');
else
    % do nothing
end

validateattributes(k, {'double'}, {'row'});
validateattributes(rho, {'double'}, {'scalar'});
validateattributes(c, {'double'}, {'scalar'});

%% Calculate the radial derivative of the Green's function
% The derivative points from Rs to R

% Distance between Rs and R
% size(dist_mat) = [size(R, 1), size(Rs, 1)]
dist_mat = pdist2(R, Rs);

% Radial derivative 
dist_mat_3D = repmat(dist_mat, 1, 1, numel(k));
k_3D = reshape(k, [1, 1, numel(k)]);

% size(rad_deriv_Rs2R) = [size(R, 1), size(Rs, 1), numel(k)]
rad_deriv_Rs2R = rad_derivative_point_source(k_3D, dist_mat_3D);

%% Calculate the velocity pointing from Rs to R
scaling_factor = 1i./rho./k_3D./c;

% size(v_Rs2R) = [size(R, 1), size(Rs, 1), numel(k)]
v_Rs2R = scaling_factor .* rad_deriv_Rs2R;

%% Calculate the velocity in x y z directions
% Calculate the unit vector in the Rs2R direction
R_expanded = repmat(R, size(Rs, 1), 1);
Rs_expanded = kron(Rs, ones(size(R, 1), 1));
Rs2R = R_expanded - Rs_expanded; % vector pointing from Rs to R
norm_Rs2R = vecnorm(Rs2R, 2, 2);
u_Rs2R = Rs2R./norm_Rs2R;

% Calculate the projection
x_u_Rs2R = reshape(u_Rs2R(:, 1), size(R, 1), size(Rs, 1));
y_u_Rs2R = reshape(u_Rs2R(:, 2), size(R, 1), size(Rs, 1));
z_u_Rs2R = reshape(u_Rs2R(:, 3), size(R, 1), size(Rs, 1));

v_x = v_Rs2R .* x_u_Rs2R;
v_y = v_Rs2R .* y_u_Rs2R;
v_z = v_Rs2R .* z_u_Rs2R;
end