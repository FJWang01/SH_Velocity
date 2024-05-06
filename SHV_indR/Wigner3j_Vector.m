function W = Wigner3j_Vector(j1, j2, j3, m1, m2, m3)
% Origin
%
% Project   : wigner3jvector
% File      : couplingcoefficients.m
% Author    : Lachlan Birnie
% Creation  : 27-09-18
% Version   : v7.20 (30-09-19)
%
% Description:
%
%   Compute Wigner 3j symbol using Racah formula on the elements of input
%   vector terms.
%
% Inputs: 
%
%   / j1 j2 j3 \
%   \ m1 m2 m3 /
%
%       - where all terms are matrices of the same size
%
% Outputs:
%
%   W   Wigner 3j symbols corresponding to the input vectors.
%
%       W(i) = / j1(i) j2(i) j3(i) \
%              \ m1(i) m2(i) m3(i) /
%
% Acknowledgements:
%
%   Inspired by Wigner3j.m by Kobi Kraus, Technion, 25-06-08
%   - avaliable at: https://au.mathworks.com/matlabcentral/fileexchange
%                   /20619-wigner3j-symbol
%

%% Check if all inputs are of the same size
if ~isequal(size(j1), size(j2), size(j3), size(m1), size(m2), size(m3))
    error('@@Wigner3j_Vector: all inputs must be of the same size');
else
    % do nothing
end

%% Convert inputs to column vectors
j1_col = j1(:);
j2_col = j2(:);
j3_col = j3(:);

m1_col = m1(:);
m2_col = m2(:);
m3_col = m3(:);

%% Mains
    % Find invalid inputs from the vectors:
    i_invalid = ( (j3_col > (j1_col + j2_col)) ...          % j3 out of interval
                | (j3_col < abs(j1_col - j2_col)) ...       % j3 out of interval
                | ((m1_col + m2_col + m3_col) ~= 0) ...     % non-conserving agular momentum 
                | (abs(m1_col) > j1_col) ...            % m is larger than j
                | (abs(m2_col) > j2_col) ...            % m is larger than j
                | (abs(m3_col) > j3_col) ...            % m is larger than j
                | ((~any([m1_col,m2_col,m3_col],2)) & (mod(j1_col+j2_col+j3_col,2)==1)) ... % m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
                );
    
    % find valid inputs:
    i_valid = ~i_invalid;
    n_valid = sum(i_valid);
    
    % return 0 if no valid inputs:
    if (n_valid == 0)
        W = zeros(size(j1));
        return;
    end
        
    % Select valid inputs:
    vj1 = j1_col(i_valid);
    vj2 = j2_col(i_valid);
    vj3 = j3_col(i_valid);
    vm1 = m1_col(i_valid);
    vm2 = m2_col(i_valid);
    vm3 = m3_col(i_valid);
    
    % compute terms:
    t1 = vj2 - vm1 - vj3;
    t2 = vj1 + vm2 - vj3;
    t3 = vj1 + vj2 - vj3;
    t4 = vj1 - vm1;
    t5 = vj2 + vm2;
    
    tmin = max( 0,  max( t1, t2 ) );
    tmax = min( t3, min( t4, t5 ) );
    
    % find largest summation:
    n_t =  max(tmax-tmin)+1;
    
    % fill in summation term matrix, pad with NaN.
    t = zeros(n_valid, n_t);
    for i = 1:n_valid
        t(i,:) = [(tmin(i):tmax(i)),nan(1,n_t-length(tmin(i):tmax(i)))];
    end
    
    % more terms?
    x = zeros(n_valid,n_t,6);
    x(:,:,1) = t;
    x(:,:,2) = t-t1;
    x(:,:,3) = t-t2;
    x(:,:,4) = t3-t;
    x(:,:,5) = t4-t;
    x(:,:,6) = t5-t;
    
    x2 = [ vj1+vj2+vj3+1 ...
           ,vj1+vj2-vj3 ...
           ,vj1-vj2+vj3 ...
           ,-vj1+vj2+vj3 ...
           ,vj1+vm1 ...
           ,vj1-vm1 ...
           ,vj2+vm2 ...
           ,vj2-vm2 ...
           ,vj3+vm3 ...
           ,vj3-vm3];
    
    % solve everything summation terms:
    sterm = (-1).^t .* exp( sum(-gammaln( x+1 ), 3) ...
                            + sum([-1,ones(1,9)].*gammaln( x2 +1 ), 2) ...
                            .* 0.5 );
                        
    % change NaNs to zero:
    sterm(isnan(sterm)) = 0;
    
    % sum over t:
    w_valid = sum(sterm, 2) .* (-1).^(vj1-vj2-vm3);
    
    % map back to full size:
    W = zeros(size(j1_col));
    W(i_valid) = w_valid;
    
    W = reshape(W, size(j1));
       
end
