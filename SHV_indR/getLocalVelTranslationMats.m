
function [glbp2vel_x_mat, glbp2vel_y_mat, glbp2vel_z_mat]  = getLocalVelTranslationMats(N, rho, c)
% This function finds the operator matrix that converts global SH coefficients of
% pressure to SH coefficients of local velocity 
%
% Inputs:
% N - SH truncation order of global SH coefficients
% rho - density of air in kg/m^3
% c - speed of sound in m/s
%
% Outputs:
% glbp2vel_x_mat - operator matrix to obtain SH coefficients of beta_10
% glbp2vel_y_mat - operator matrix to obtain SH coefficients of
%                      beta_1{-1}
% glbp2vel_z_mat - operator matrix to obtain SH coefficients of beta_11
%   size(glbp2vel_x_mat_mat) = size(glbp2vel_y_mat) =
%   size(glbp2vel_z_mat) = [N^2, (N+1)^2]

%% Setup 
[n, m, ~] = getDegreeOrderPairs(N);
[n_ddash, m_ddash, ~] = getDegreeOrderPairs(N-1);

[n_mat, n_ddash_mat] = meshgrid(n, n_ddash);
[m_mat, m_ddash_mat] = meshgrid(m, m_ddash);

ones_mat = ones(size(n_mat));
zeros_mat = zeros(size(n_mat));

%% Calculate common scaling factor
scaling_factor_1 = 4*pi*1i.^(1+n_ddash_mat-n_mat).*(-1).^(m_mat);
scaling_factor_2 = sqrt(3*(2*n_mat + 1) .* (2*n_ddash_mat + 1)/4/pi);

scaling_factor = scaling_factor_1 .* scaling_factor_2;

%% beta_10
W1_beta_10 = Wigner3j_Vector(n_mat, ones_mat, n_ddash_mat, zeros_mat, zeros_mat, zeros_mat);
W2_beta_10 = Wigner3j_Vector(n_mat, ones_mat, n_ddash_mat, -m_mat, zeros_mat, m_mat);

W_mult_beta_10 = scaling_factor .* W1_beta_10 .* W2_beta_10;

mask_beta_10 = m_mat == m_ddash_mat;

% calculate operator matrix
glbp2beta10_mat = zeros(size(W_mult_beta_10));
glbp2beta10_mat(mask_beta_10) = W_mult_beta_10(mask_beta_10); 

%% beta_1{-1}
W1_beta_1neg1 = Wigner3j_Vector(n_mat, ones_mat, n_ddash_mat, zeros_mat, zeros_mat, zeros_mat);
W2_beta_1neg1 = Wigner3j_Vector(n_mat, ones_mat, n_ddash_mat, -m_mat, -ones_mat, m_mat+ones_mat);

W_mult_beta_1neg1 = scaling_factor .* W1_beta_1neg1 .* W2_beta_1neg1;

mask_beta_1neg1 = (m_mat + 1) == m_ddash_mat;

% calculate operator matrix
glbp2beta1neg1_mat = zeros(size(W_mult_beta_1neg1));
glbp2beta1neg1_mat(mask_beta_1neg1) = W_mult_beta_1neg1(mask_beta_1neg1);

%% beta_11
W1_beta_11 = Wigner3j_Vector(n_mat, ones_mat, n_ddash_mat, zeros_mat, zeros_mat, zeros_mat);
W2_beta_11 = Wigner3j_Vector(n_mat, ones_mat, n_ddash_mat, -m_mat, ones_mat, m_mat-ones_mat);

W_mult_beta_11 = scaling_factor .* W1_beta_11 .* W2_beta_11;

mask_beta_11 = (m_mat - 1) == m_ddash_mat;

% calculate operator matrix
glbp2beta11_mat = zeros(size(W_mult_beta_11));
glbp2beta11_mat(mask_beta_11) = W_mult_beta_11(mask_beta_11);

%% Velocity in the x-direction
vel_scaling_factor = 1/3 *1i/rho/c;
glbp2vel_x_mat = vel_scaling_factor * (sqrt(3/8/pi) * glbp2beta1neg1_mat - sqrt(3/8/pi) * glbp2beta11_mat);
glbp2vel_y_mat = vel_scaling_factor * (-sqrt(3/8/pi) *1i * glbp2beta1neg1_mat - sqrt(3/8/pi) *1i * glbp2beta11_mat);
glbp2vel_z_mat = vel_scaling_factor * sqrt(3/4/pi) * glbp2beta10_mat;
end