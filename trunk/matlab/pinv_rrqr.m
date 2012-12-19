function [pinv_A, varargout] = pinv_rrqr(A)

% [pinv_A, rnull_optional, lnull_optional] = pinv_rrqr(A)
%
% Computes Moore-Penrose pseudoinverse of A using a rank-revealing QR
% factorization and optionally returns an orthonormal basis for right
% and left null-space.  This can be significantly faster than using
% SVD for computing the pseudoinverse.  Moreover, A can be sparse or
% dense.

n_opt_out = max(nargout - 1, 0);

assert(n_opt_out <= 2, 'pinv_rrqr: Too many optional output arguments (%d)', n_opt_out);
assert(ndims(A) <= 2,   'pinv_rrqr: ndims(A) > 2');

[Q R P] = qr(full(A));
abs_diag_R = abs(diag(R));

tol = 100 * min(size(A)) * max(abs_diag_R) * eps(class(A));  % MAGIC CONSTANT
rnk = sum(tol < abs_diag_R);

[n1 n2] = size(R);

% linsolve not implemented in octave-3.6.2, so just default to
% '\' for octave.

if(exist('OCTAVE_VERSION','builtin') ~= 0)
    S = R(1:rnk,1:rnk) \ R(1:rnk,rnk+1:end);
    T = R(1:rnk,1:rnk) \ Q(:,1:rnk)';
    X = ((S'*S + speye(size(S,2))) \ S') * T;
else
    opts_R11.UT = true;
    S = linsolve(R(1:rnk,1:rnk), R(1:rnk,rnk+1:end), opts_R11);
    T = linsolve(R(1:rnk,1:rnk), Q(:,1:rnk)', opts_R11);
    opts_STS_I.POSDEF = true;
    opts_STS_I.SYM = true;
    if(isequal(class(A), 'single'))
        X = linsolve(S'*S + eye(size(S,2),'single'), S', opts_STS_I) * T;
    else
        X = linsolve(S'*S + speye(size(S,2)), S', opts_STS_I) * T;
    end
end

if(isequal(class(A), 'single'))
    sP = P;
else
    sP = sparse(P);
end

pinv_A = sP * [T - S*X; X];

right_null_size = n2 - rnk;
left_null_size  = n1 - rnk;

if(0 < n_opt_out)
    if(0 < right_null_size)
        [rnull, junk, junk] = qr(sP * [S; -eye(right_null_size)], 0);
        size(junk); % no unused variable warning
    else
        rnull = zeros(n2,0);
    end

    varargout(1) = {rnull};
end

if(1 < n_opt_out)
    if(0 < left_null_size)
        lnull = Q(:,end - (left_null_size - 1): end);
    else
        lnull = zeros(n1,0);
    end

    varargout(2) = {lnull};
end

