
function dG_dr = rad_derivative_point_source(k, r)
% dG_dr = rad_derivative_point_source(k, r)
%
% This function calculates the radial derivative of the Green's function
% exp(-ikr)/4/pi/r
%
% Inputs:
% k - wavenumber
% r - radius
%
% Output
% dG_dr - d/dr [exp(-ikr)/4/pi/r]

%% Calculate the derivative using quotient rule
term_1 = exp(-1i * k .* r)./r;
term_2 = (-1i * k .* r - 1)./r;

dG_dr = 1/4/pi * term_1 .* term_2;
end