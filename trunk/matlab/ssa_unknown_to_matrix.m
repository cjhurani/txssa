function X = ssa_unknown_to_matrix(x, A_id)

% X = ssa_unknown_to_matrix(x, A_id)
%
% This function takes the real vector of unknowns 'x' and a sparse matrix
% A_id with locations where values from 'x' should be kept.  It returns
% a matrix of size and field same as A_id.  A_id can be real or complex.

class_A_id = class(A_id);
is_real_A_id = isreal(A_id);
[n1 n2] = size(A_id);

if(is_real_A_id)
    [A_id_i A_id_j A_id_v] = find(A_id);
    A_id_nnz = nnz(A_id);

    X_val = zeros(A_id_nnz,1,class_A_id);

    % copy values to output
    for k = 1:A_id_nnz
        X_val(k) = x(A_id_v(k));
    end

    X = sparse(A_id_i, A_id_j, X_val, n1, n2);
else
    
    real_A_id = real(A_id);
    imag_A_id = imag(A_id);

    [real_A_id_i real_A_id_j real_A_id_v] = find(real_A_id);
    [imag_A_id_i imag_A_id_j imag_A_id_v] = find(imag_A_id);
    
    real_A_id_nnz = nnz(real_A_id);
    imag_A_id_nnz = nnz(imag_A_id);
    
    real_X_val = zeros(real_A_id_nnz,1,class_A_id);
    imag_X_val = zeros(imag_A_id_nnz,1,class_A_id);

    % copy values to output
    for k = 1:real_A_id_nnz
        real_X_val(k) = x(real_A_id_v(k));
    end

    for k = 1:imag_A_id_nnz
        imag_X_val(k) = x(imag_A_id_v(k));
    end

    real_X = sparse(real_A_id_i, real_A_id_j, real_X_val, n1, n2);
    imag_X = sparse(imag_A_id_i, imag_A_id_j, imag_X_val, n1, n2);

    X = real_X + 1i * imag_X;
end
