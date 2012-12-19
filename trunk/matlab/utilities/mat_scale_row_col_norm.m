function [A_scaled, varargout] = mat_scale_row_col_norm(A, varargin)

% [A_scaled, d_left_opt, d_right_opt, iters_opt] =
%   mat_scale_row_col_norm(A [, p] [, tol] [, max_iters])
%
% Compute a left and right diagonal scaling of a matrix A so that
% after scaling each row's and each column's vector p-norms are almost
% 1 for square matrix and close to each other and 1 for rectangular
% matrix.  This is an iterative algorithm.
%
% varargin can be p, tol, max_iters.
% p must be is a real number in [1,inf].  Default is inf.
% tol must be a 'small' non-negative real number.  Default is 1E-2.
% tol measures how much are the row norms far from being equal to
% each other.  It is a relative tolerance, which means it is independent
% of the magnitude of matrix enrties.  The same applies to column norms.
% max_iters must be a non-negative integer or inf.  Default is 10.
%
% The output consists optionally of d_left_opt and d_right_opt that are
% diagonals of the left and right scaling matrices, respectively.  The
% scaled matrix can be computed using
% sparse(1:n1,1:n1,d_left_opt,n1,n1) * A * sparse(1:n2,1:n2,d_right_opt,n2,n2)
%
% Number of iterations taken are returned in a fourth optional output
% variable.
%
% Algorithm modified from
% http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.9.1143
% or ftp://ftp.numerical.rl.ac.uk/pub/reports/drRAL2001034.pdf
%
% (A Scaling Algorithm to Equilibrate Both Rows and Columns Norms
% in Matrices (2001) by Daniel Ruiz)
% Modifications deal with complex matrices, and convergence criteria for
% rectangular matrices when p ~= inf.
%
% If A is Hermitian or skew-Hermitian or complex-symmetric (symmetric or
% skew-symmetric when real), then the scaled matrix also satisfies the
% same properties and d_left_opt and d_right_opt are equal.  If A is
% multiplied by a non-zero number, then the scaled matrix for the
% multiple of A does not change.  This holds only for the exact solution
% and may not hold if the solution is computed iteratively as done below.

n_opt_out = max(nargout - 1, 0);
assert(n_opt_out <= 3, 'mat_scale_row_col_norm: Too many optional output arguments (%d)', n_opt_out);

optargin = size(varargin, 2);
assert(optargin <= 3, 'mat_scale_row_col_norm: Too many optional input arguments (%d)', optargin);

if(optargin == 0 || isempty(varargin{1}))
    p = inf;
else
    p = varargin{1};
    assert(numel(p) == 1, 'mat_scale_row_col_norm: second input must be a number');
    assert(isreal(p),     'mat_scale_row_col_norm: second input must be real');
    assert(1 <= p,        'mat_scale_row_col_norm: second input (p = %f) must not be less than 1', p);
end

if(optargin <= 1 || isempty(varargin{2}))
    tol = 1E-2;
else
    tol = varargin{2};
    assert(numel(tol) == 1, 'mat_scale_row_col_norm: third input must be a number');
    assert(isreal(tol),     'mat_scale_row_col_norm: third input must be real');
    assert(0 <= tol,        'mat_scale_row_col_norm: third input (tol = %f) must not be less than 0', tol);
end

if(optargin <= 2 || isempty(varargin{3}))
    max_iters = 10;
else
    max_iters = varargin{3};
    assert(numel(max_iters) == 1, 'mat_scale_row_col_norm: fourth input must be a number');
    assert(isreal(max_iters),     'mat_scale_row_col_norm: fourth input must be real');
    assert(0 <= max_iters,        'mat_scale_row_col_norm: fourth input (max_iters = %d) must not be less than 0', max_iters);
    assert(~(max_iters <= int32(inf) && max_iters ~= int32(max_iters)), 'mat_scale_row_col_norm: fourth input (max_iters = %d) must be an integral value or infinity', max_iters);
end

[n1 n2] = size(A);

d_left_opt = ones(n1, 1);
d_right_opt = ones(n2, 1);

if(n1 == 0 || n2 == 0)
    A_scaled = A;
    if(0 < n_opt_out)
        varargout(1) = {d_left_opt};
    end

    if(1 < n_opt_out)
        varargout(2) = {d_right_opt};
    end

    if(2 < n_opt_out)
        varargout(3) = {0};
    end
    
    return
end

DR = zeros(n1, 1);
DC = zeros(n2, 1);

if(p == inf)
    DR_lim = 1;
else
    DR_lim = (n2/n1)^(0.25/p);
end

DC_lim = 1/DR_lim;

iters = 0;

for i = 1:max_iters

    if(p == 1)
        for i1 = 1:n1
            DR(i1) = sqrt(sum(abs(A(i1,:))));
        end

        for i2 = 1:n2
            DC(i2) = sqrt(sum(abs(A(:,i2))));
        end
    elseif(p == inf)
        for i1 = 1:n1
            DR(i1) = sqrt(max(abs(A(i1,:))));
        end

        for i2 = 1:n2
            DC(i2) = sqrt(max(abs(A(:,i2))));
        end
    else
        for i1 = 1:n1
            DR(i1) = sqrt(sum(abs(A(i1,:)) .^ p)^(1/p));
        end

        for i2 = 1:n2
            DC(i2) = sqrt(sum(abs(A(:,i2)) .^ p)^(1/p));
        end
    end

    % check zero row / columns when coming in the loop
    % zero row/cols cause two problems -
    % 1. divide by zero
    % 2. need to change limiting values DR_lim and DC_lim
    % If we find zero row/cols, we remove those rows/cols
    % and call this function again with the remaining matrix.
    % After computing the scaling and scaled matrix, we
    % return the big-size versions.
    if i == 1
        if(~isempty(find(DR == 0,1)) || ~isempty(find(DC == 0,1)))

            % works for matrix with no entries too.

            [scaled_A_loc, d_left_loc, d_right_loc, iters] = ...
                mat_scale_row_col_norm(A(0 < DR, 0 < DC), p, tol, max_iters);

            % create appropriately sized output

            if(~isempty(d_left_loc))
                d_left_opt(DR == 0) = 1;
            end

            if(~isempty(d_right_loc))
                d_right_opt(DC == 0) = 1;
            end

            d_left_opt(0 < DR) = d_left_loc;
            d_right_opt(0 < DC) = d_right_loc;
            A(0 < DR, 0 < DC) = scaled_A_loc;
            break;
        end
    end

    % If we're here, all entries of DR and DC are positive.
    if(norm(DR - DR_lim, inf) <= tol && norm(DC - DC_lim, inf) <= tol)
        break;
    end

    % take inverses once.
    DR = 1 ./ DR;
    DC = 1 ./ DC;

    A = sparse(1:n1,1:n1,DR,n1,n1) * A * sparse(1:n2,1:n2,DC,n2,n2);

    d_left_opt = d_left_opt .* DR;
    d_right_opt = d_right_opt .* DC;

    iters = iters + 1;
end

A_scaled = A;

if(0 < n_opt_out)
    varargout(1) = {d_left_opt};
end

if(1 < n_opt_out)
    varargout(2) = {d_right_opt};
end

if(2 < n_opt_out)
    varargout(3) = {iters};
end

