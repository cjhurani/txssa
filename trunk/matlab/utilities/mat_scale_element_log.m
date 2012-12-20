function [A_scaled, varargout] = mat_scale_element_log(A, varargin)

% [A_scaled, d_left_opt, d_right_opt, iters_opt] =
%   mat_scale_element_log(A [, tol] [, max_iters])
%
% Compute a left and right diagonal scaling of a matrix A so that
% after scaling each entry is close to 1.  This is an iterative algorithm.
%
% varargin can be tol, max_iters.
% tol must be a 'small' non-negative real number.  Default is 1E-6.
% tol is used in an iterative preconditioned conjugate gradient algorithm.
% It is a relative tolerance, which means it is independent of the magnitude
% of matrix enrties.
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
% Algorithm inspired from
% On the automatic scaling of matrices for Gaussian elimination
% Curtis and Reid, J. Inst. Maths. Applics. (1972), 10, pp. 118-124.
% and
% http://www.hsl.rl.ac.uk/specs/mc29.pdf 
%
% If A is Hermitian or skew-Hermitian or complex-symmetric (symmetric or
% skew-symmetric when real), then the scaled matrix also satisfies the
% same properties and d_left_opt and d_right_opt are equal.  If A is
% multiplied by a non-zero number, then the scaled matrix for the
% multiple of A does not change.  This holds only for the exact solution
% and may not hold if the solution is computed iteratively as done below.

n_opt_out = max(nargout - 1, 0);
assert(n_opt_out <= 3, 'mat_scale_element_log: Too many optional output arguments (%d)', n_opt_out);

optargin = size(varargin, 2);
assert(optargin <= 2, 'mat_scale_element_log: Too many optional input arguments (%d)', optargin);

if(optargin == 0 || isempty(varargin{1}))
    tol = 1E-6;
else
    tol = varargin{1};
    assert(numel(tol) == 1, 'mat_scale_element_log: second input must be a number');
    assert(isreal(tol),     'mat_scale_element_log: second input must be real');
    assert(0 <= tol,        'mat_scale_element_log: second input (tol = %f) must not be less than 0', tol);
end

if(optargin <= 1 || isempty(varargin{2}))
    max_iters = 10;
else
    max_iters = varargin{2};
    assert(numel(max_iters) == 1, 'mat_scale_element_log: third input must be a number');
    assert(isreal(max_iters),     'mat_scale_element_log: third input must be real');
    assert(0 <= max_iters,        'mat_scale_element_log: third input (max_iters = %d) must not be less than 0', max_iters);
    assert(~(max_iters <= int32(inf) && max_iters ~= int32(max_iters)), 'mat_scale_element_log: third input (max_iters = %d) must be an integral value or infinity', max_iters);
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

pat = sparse(A ~= 0);

H = [ ...
   spdiags(sum(pat,2),0,n1,n1) pat;
   pat'                        spdiags(sum(pat,1)',0,n2,n2)];

B = A;
B(A ~= 0) = log(abs(A(A ~= 0)));

rhs = - full([sum(B,2); sum(B,1)']);

[x,flag,relres,iters] = pcg(H, rhs, tol, max_iters);

d_left_opt = exp(x(1:n1));
d_right_opt = exp(x(n1+1:end));

A_scaled = sparse(1:n1,1:n1,d_left_opt,n1,n1) * A * sparse(1:n2,1:n2,d_right_opt,n2,n2);

if(0 < n_opt_out)
    varargout(1) = {d_left_opt};
end

if(1 < n_opt_out)
    varargout(2) = {d_right_opt};
end

if(2 < n_opt_out)
    varargout(3) = {iters};
end

