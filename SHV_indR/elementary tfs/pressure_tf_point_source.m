
function p = pressure_tf_point_source(R, Rs, k)
% p = pressure_tf_point_source(R, Rs, k)
% 
% This function calculates the pressure at point receiver R due to point 
% source at Rs.
%
% Inputs:
% Rs - locations of point sources in Cartesian coordinates
% R - locations of point receivers in Cartesian coordinates
% k - wavenumbers, must be a row vector
%
% Outputs:
% p - pressure at point receiver R due to point source at Rs
%   size(p) = [size(R, 1), size(Rs, 1), numel(k)]

%% Check the dimensions of inputs
if ~isequal(size(R, 2), size(Rs, 2), 3)
    error('@@ pressure_tf_point_source: R and Rs must have three columns');
else
    % do nothing
end

validateattributes(k, {'double'}, {'row'});

%% Calculate pressure at point receiver R
% Distance between Rs and R
% size(dist_mat) = [size(R, 1), size(Rs, 1)]
dist_mat = pdist2(R, Rs);

k_3D = reshape(k, [1, 1, numel(k)]);
dist_mat_3D = repmat(dist_mat, 1, 1, numel(k));

p = exp(-1i * k_3D .* dist_mat_3D)/4/pi./dist_mat_3D;

end