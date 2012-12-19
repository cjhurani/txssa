function [X, varargout] = ssa_impose_action(Y, A, R_mat, L_mat, varargin)

% [X, iters_optional, resvec_optional] = ...
%     ssa_impose_action(Y, A, R_mat, L_mat [, tol] [, max_iters])
%
% 'ssa' stands for sparse spectral approximation.
%
% Computes a matrix X of the same size, field, and pattern, as Y such that X
% approximates Y and has same right and left actions on L_mat and R_mat as A,
% respectively.  All arguments must be real or complex.  We use the Uzawa-CG
% iterative algorithm.  It helps to have orthonormal R_mat and L_mat for
% faster convergence.
% 
% The computed X is such that (X - A) * R_mat = 0, and (X' - A') * L_mat = 0 (up
% to tol). Thus, R_mat and L_mat must be of appropriate sizes.  They could be
% empty too.  Typically, R_mat and L_mat form a basis for the right and left
% null-spaces of A, respectively.  The optimization process minimizes
%
%   (1/2) * || X - Y ||_F^2
%
% such that X and Y have the same sparsity pattern and
% (X - A) * R_mat = 0, and (X' - A') * L_mat = 0.
%
% varargin can be tol, max_iters.
% tol must be a 'small' positive real number.  Default is
% eps(class(Y)).  Iterations stop when
% || (X  - A ) * R_mat ||_F <= tol * size(X,2) * || (X  - A ) ||_F * || R_mat ||_F
% and
% || (X' - A') * L_mat ||_F <= tol * size(X,1) * || (X' - A') ||_F * || L_mat ||_F
%
% max_iters must be a positive integer.  Default is 1000.
% Optionally, the routine also returns the number of iterations actually
% taken.  resvec (residual vector) is returned if asked for.  It is of size
% (iters + 1) x 2 and contains Frobenius norms of the mismatch (X - A) * R_mat and
% (X' - A') * L_mat.  The first row has the initial residual norms and
% the last row has the final residual norms (computed using the returned X).
%
% Note that since we preserve the pattern, the zeros in a dense Y are maintained
% as zeros.

[n1 n2] = size(Y);

n_opt_out = max(nargout - 1, 0);

assert(n_opt_out <= 2, 'ssa_impose_action: Too many optional output arguments (%d)', n_opt_out);

optargin = size(varargin, 2);

assert(optargin <= 2, 'ssa_impose_action: Too many optional input arguments (%d)', optargin);

if(optargin == 0 || isempty(varargin{1}))
    tol = eps(class(Y));
else
    tol = varargin{1};
    assert(numel(tol) == 1 && isreal(tol) && 0 < tol, 'ssa_impose_action: tol should be a real number in (0,inf]');
end

if(optargin <= 1 || isempty(varargin{2}))
    max_iters = 1000;
else
    max_iters = varargin{2};
    assert(numel(max_iters) == 1 && isreal(max_iters) && 0 < max_iters ...
        && (max_iters == inf || max_iters == int32(max_iters)), 'ssa_impose_action: max_iters should be a positive integral number or inf');
end

is_real_Y = isreal(Y);

if(isempty(R_mat))
    if(is_real_Y)
        R_mat = zeros(n2,0);
    else
        R_mat = complex(zeros(n2,0));
    end
end

if(isempty(L_mat))
    if(is_real_Y)
        L_mat = zeros(n1,0);
    else
        L_mat = complex(zeros(n1,0));
    end
end

assert(all(size(Y) == size(A)), 'ssa_impose_action: Y not the same size as A');
assert(~xor(is_real_Y, isreal(A)), 'ssa_impose_action: A must be a of the same field as Y');

[nR1 nR2] = size(R_mat);
[nL1 nL2] = size(L_mat);

assert(n2 == nR1, 'ssa_impose_action: Y * R_mat not possible');
assert(n1 == nL1, 'ssa_impose_action: ctranspose(Y) * L_mat not possible');
assert(isempty(R_mat) || ~xor(is_real_Y, isreal(R_mat)), 'ssa_impose_action: R_mat must be a of the same field as Y');
assert(isempty(L_mat) || ~xor(is_real_Y, isreal(L_mat)), 'ssa_impose_action: L_mat must be a of the same field as Y');

% Quantities that don't change when iterating.
norm_R_mat = norm(R_mat, 'fro');
norm_L_mat = norm(L_mat, 'fro');

A_R_mat = A  * R_mat;
A_L_mat = A' * L_mat;

if(is_real_Y)
    nnz_Y = nnz(Y);
    [Y_i Y_j] = find(Y);
else
    real_R_mat = real(R_mat);
    imag_R_mat = imag(R_mat);

    real_L_mat = real(L_mat);
    imag_L_mat = imag(L_mat);

    real_Y = real(Y);
    imag_Y = imag(Y);

    nnz_real_Y = nnz(real_Y);
    nnz_imag_Y = nnz(imag_Y);

    [real_Y_i real_Y_j] = find(real_Y);
    [imag_Y_i imag_Y_j] = find(imag_Y);
end

% Quantities that do change when iterating.

% This will work even in complex case. No need to create complex matrices.
Lag_R = zeros(n1, nR2);
Lag_L = zeros(n2, nL2);

X = Y;

d_R = X  * R_mat - A_R_mat;
d_L = X' * L_mat - A_L_mat;

q_R = - d_R;
q_L = - d_L;

if(is_real_Y)
    d_R_proj_vec = zeros(nnz_Y, 1);
    d_L_proj_vec = zeros(nnz_Y, 1);
    d_proj_vec   = zeros(nnz_Y, 1);
else
    real_d_R_proj_vec = zeros(nnz_real_Y, 1);
    imag_d_R_proj_vec = zeros(nnz_imag_Y, 1);

    real_d_L_proj_vec = zeros(nnz_real_Y, 1);
    imag_d_L_proj_vec = zeros(nnz_imag_Y, 1);

    real_d_proj_vec = zeros(nnz_real_Y, 1);
    imag_d_proj_vec = zeros(nnz_imag_Y, 1);
end

if(1 < n_opt_out)
    resvec = zeros(max_iters + 1,2);
end

% start iterations
iters = 0;

for i = 1:max_iters
    norm_X_m_A = norm(X - A, 'fro');
    norm_XT_m_AT = norm_X_m_A;
    norm_d_R = norm(d_R, 'fro');
    norm_d_L = norm(d_L, 'fro');

    if(norm_d_R <= tol * n2 * norm_X_m_A * norm_R_mat && norm_d_L <= tol * n1 * norm_XT_m_AT * norm_L_mat)
        break;
    end

    if(1 < n_opt_out)
        resvec(i,1) = norm_d_R;
        resvec(i,2) = norm_d_L;
    end
    
    if(is_real_Y)
        for k = 1:nnz_Y
            d_R_proj_vec(k) = dot(R_mat(Y_j(k),:), d_R(Y_i(k),:));
            d_L_proj_vec(k) = dot(L_mat(Y_i(k),:), d_L(Y_j(k),:));
        end

        d_proj_vec = d_R_proj_vec + d_L_proj_vec;
    else
        real_d_R = real(d_R);
        imag_d_R = imag(d_R);

        real_d_L = real(d_L);
        imag_d_L = imag(d_L);

        for k = 1:nnz_real_Y
            real_Y_i_k = real_Y_i(k);
            real_Y_j_k = real_Y_j(k);

            real_d_R_proj_vec(k) = ...
                dot(real_R_mat(real_Y_j_k,:), real_d_R(real_Y_i_k,:)) ...
              + dot(imag_R_mat(real_Y_j_k,:), imag_d_R(real_Y_i_k,:));
            real_d_L_proj_vec(k) = ...
                dot(real_L_mat(real_Y_i_k,:), real_d_L(real_Y_j_k,:)) ...
              + dot(imag_L_mat(real_Y_i_k,:), imag_d_L(real_Y_j_k,:));
        end

        for k = 1:nnz_imag_Y
            imag_Y_i_k = imag_Y_i(k);
            imag_Y_j_k = imag_Y_j(k);

            imag_d_R_proj_vec(k) = ...
                dot(real_R_mat(imag_Y_j_k,:), imag_d_R(imag_Y_i_k,:)) ...
              - dot(imag_R_mat(imag_Y_j_k,:), real_d_R(imag_Y_i_k,:));
            imag_d_L_proj_vec(k) = - (...
                dot(real_L_mat(imag_Y_i_k,:), imag_d_L(imag_Y_j_k,:)) ...
              - dot(imag_L_mat(imag_Y_i_k,:), real_d_L(imag_Y_j_k,:)));
        end

        real_d_proj_vec = real_d_R_proj_vec + real_d_L_proj_vec;
        imag_d_proj_vec = imag_d_R_proj_vec + imag_d_L_proj_vec;
    end

    norm_q_R = norm(q_R, 'fro');
    norm_q_L = norm(q_L, 'fro');

    if(is_real_Y)
        alpha = (norm_q_R^2 + norm_q_L^2) / norm(d_proj_vec, 'fro')^2;
    else
        alpha = (norm_q_R^2 + norm_q_L^2) / ...
            (norm(real_d_proj_vec, 'fro')^2 + norm(imag_d_proj_vec, 'fro')^2);
    end

    Lag_R = Lag_R + alpha * d_R;
    Lag_L = Lag_L + alpha * d_L;

    %Do it this way because MATLAB may change sparsity pattern of X if an element becomes 0.

    if(is_real_Y)
        for k = 1:nnz_Y
            X(Y_i(k), Y_j(k)) = X(Y_i(k), Y_j(k)) - alpha * d_proj_vec(k);
        end
    else
        for k = 1:nnz_real_Y
            X(real_Y_i(k), real_Y_j(k)) = ...
            X(real_Y_i(k), real_Y_j(k)) - alpha * real_d_proj_vec(k);
        end

        for k = 1:nnz_imag_Y
            X(imag_Y_i(k), imag_Y_j(k)) = ...
            X(imag_Y_i(k), imag_Y_j(k)) - complex(0,alpha * imag_d_proj_vec(k));
        end
    end

    norm2_q_RL_old = norm(q_R, 'fro')^2 + norm(q_L, 'fro')^2;

    q_R = A_R_mat - X  * R_mat;
    q_L = A_L_mat - X' * L_mat;

    norm2_q_RL_new = norm(q_R, 'fro')^2 + norm(q_L, 'fro')^2;

    beta = norm2_q_RL_new / norm2_q_RL_old;
    d_R = beta * d_R - q_R;
    d_L = beta * d_L - q_L;

    iters = iters + 1;
end

if(0 < n_opt_out)
    varargout(1) = {iters};
end

if(1 < n_opt_out)
    resvec(iters + 1,1) = norm((X-A)*R_mat,'fro');
    resvec(iters + 1,2) = norm((X'-A')*L_mat,'fro');
    resvec = resvec(1:iters + 1,:);
    varargout(2) = {resvec};
end


