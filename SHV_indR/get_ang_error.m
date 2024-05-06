
function [ang_err_real, ang_err_imag, DOT_real, DOT_imag] = get_ang_error(des, reprod)
% [ang_err_real, ang_err_imag, DOT_real, DOT_imag] = get_ang_error(des, reprod)
% 
% This function calculates the angular error between the desired vector and
% the reproduced vector

%% Check inputs
if ~isequal(size(des, 2), size(reprod, 2))
    error('@@ get_ang_error: des and reprod must have the same number of columns');
else 
    % do nothing 
end

if ~isequal(size(des, 1), size(reprod, 1))
    error('@@ get_ang_error: des and reprod must be of the same size');
else
    % do nothing
end

%% Calculate DOT
% Real part
num_real = dot(real(des), real(reprod), 2); 
denom_real = vecnorm(real(des), 2, 2) .* vecnorm(real(reprod), 2, 2);
DOT_real = num_real./denom_real;

% Imaginary part
num_imag = dot(imag(des), imag(reprod), 2); 
denom_imag = vecnorm(imag(des), 2, 2) .* vecnorm(imag(reprod), 2, 2);
DOT_imag = num_imag./denom_imag;

%% Calculate angular error
ang_err_real = acos(DOT_real)/pi;
ang_err_imag = acos(DOT_imag)/pi;
end