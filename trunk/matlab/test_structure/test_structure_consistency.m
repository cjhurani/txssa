function test_structure_consistency(sizes)

H = ...
{
    'centrosymmetric',              @rand_centrosymmetric,              @diff_centrosymmetric,        @default_nearest_even_ranked,           [];
    'centrosymmetric_complex',      @rand_centrosymmetric_complex,      @diff_centrosymmetric,        @default_nearest_even_ranked,           [];
    'circulant',                    @rand_circulant,                    @diff_circulant,              @default_make_circulant_rank_deficient, [];
    'circulant_complex',            @rand_circulant_complex,            @diff_circulant,              @default_nearest_even_ranked,           [];
    'complex_symmetric',            @rand_complex_symmetric,            @diff_complex_symmetric,      @default_nearest_even_ranked,           [];
    'hamiltonian',                  @rand_hamiltonian,                  @diff_hamiltonian,            @default_nearest_even_ranked,           @is_even;
    'hamiltonian_complex',          @rand_hamiltonian_complex,          @diff_hamiltonian,            @default_nearest_even_ranked,           @is_even;
    'hermitian',                    @rand_hermitian,                    @diff_hermitian,              @default_nearest_even_ranked,           [];
    'persymmetric',                 @rand_persymmetric,                 @diff_persymmetric,           @default_nearest_even_ranked,           [];
    'persymmetric_complex',         @rand_persymmetric_complex,         @diff_persymmetric,           @default_nearest_even_ranked,           [];
    'skew_centrosymmetric',         @rand_skew_centrosymmetric,         @diff_skew_centrosymmetric,   @default_nearest_even_ranked,           [];
    'skew_centrosymmetric_complex', @rand_skew_centrosymmetric_complex, @diff_skew_centrosymmetric,   @default_nearest_even_ranked,           [];
    'skew_circulant',               @rand_skew_circulant,               @diff_skew_circulant,         @default_make_circulant_rank_deficient, [];
    'skew_circulant_complex',       @rand_skew_circulant_complex,       @diff_skew_circulant,         @default_nearest_even_ranked,           [];
    'skew_complex_symmetric',       @rand_skew_complex_symmetric,       @diff_skew_complex_symmetric, @default_nearest_even_ranked,           [];
    'skew_hamiltonian',             @rand_skew_hamiltonian,             @diff_skew_hamiltonian,       @default_nearest_even_ranked,           @is_even;
    'skew_hamiltonian_complex',     @rand_skew_hamiltonian_complex,     @diff_skew_hamiltonian,       @default_nearest_even_ranked,           @is_even;
    'skew_hermitian',               @rand_skew_hermitian,               @diff_skew_hermitian,         @default_nearest_even_ranked,           [];
    'skew_persymmetric',            @rand_skew_persymmetric,            @diff_skew_persymmetric,      @default_nearest_even_ranked,           [];
    'skew_persymmetric_complex',    @rand_skew_persymmetric_complex,    @diff_skew_persymmetric,      @default_nearest_even_ranked,           [];
    'skew_symmetric',               @rand_skew_symmetric,               @diff_skew_symmetric,         @default_nearest_even_ranked,           [];
    'symmetric',                    @rand_symmetric,                    @diff_symmetric,              @default_nearest_even_ranked,           [];
};

transformations = ...
{
    'identity',                           @identity;
    'transpose',                          @transpose;
    'Hermitian-transpose',                @ctranspose;
    'default_mat_scale_row_col_norm',     @mat_scale_row_col_norm;
    'default_mat_scale_element_log',      @mat_scale_element_log;
    'pinv_rrqr',                          @pinv_rrqr;
    'default_ssa_compute_bin_no_null',    @default_ssa_compute_bin_no_null;
    'default_ssa_compute_no_bin_no_null', @default_ssa_compute_no_bin_no_null;
    'default_ssa_compute_exact',          @default_ssa_compute_exact;
};

num_types = size(H,1);
num_transformations = size(transformations,1);

max_rel_errors = zeros(num_types, 2*num_transformations + 1);

for i = 1:num_types
    for sz = sizes
        if(isempty(H{i,5}) || H{i,5}(sz))
            A = H{i,2}(sz);
            rank_def_A = H{i,4}(A);
            norm_trans_rank_def_A = norm(rank_def_A,'fro');

            for i_trans = 1:num_transformations
                trans = transformations{i_trans,2};
                trans_A = trans(A);
                trans_rank_def_A = trans(A);

                norm_trans_A = norm(trans_A,'fro');
                norm_trans_A_diff = H{i,3}(trans_A);
                norm_trans_rank_def_A_diff = H{i,3}(trans_rank_def_A);

                rel_error = norm_trans_A_diff/norm_trans_A;
                rel_error_rank_def = norm_trans_rank_def_A_diff/norm_trans_rank_def_A;

                if(max_rel_errors(i, i_trans) < rel_error)
                    max_rel_errors(i, i_trans) = rel_error;
                end

                if(max_rel_errors(i, i_trans + num_transformations) < rel_error_rank_def)
                    max_rel_errors(i, i_trans + num_transformations) = rel_error_rank_def;
                end
            end

            tmp = H{i,2}(sz);
            imposed_action_A = ssa_impose_action(tmp, rank_def_A, null(rank_def_A), null(rank_def_A'));
            imposed_action_A_diff = H{i,3}(imposed_action_A);

            rel_error_after_imposed = imposed_action_A_diff/norm(tmp,'fro');

            if(max_rel_errors(i, 1 + 2*num_transformations) < rel_error_after_imposed)
                max_rel_errors(i, 1 + 2*num_transformations) = rel_error_after_imposed;
            end
        end
    end
end

log_errors = round(log10(max_rel_errors + eps('double')));
log_errors(log_errors > -10) = 1;
log_errors(log_errors <= -10) = 0;

disp('________________________________________________________');
disp('A: identity');
disp('B: transpose');
disp('C: Hermitian-transpose');
disp('D: mat_scale_row_col_norm');
disp('E: mat_scale_element_log');
disp('F: pinv_rrqr');
disp('G: ssa_compute_bin_no_null');
disp('H: ssa_compute_no_bin_no_null');
disp('I: ssa_compute_exact');
disp('J: identity (for rank-deficient)');
disp('K: transpose (for rank-deficient)');
disp('L: Hermitian-transpose (for rank-deficient)');
disp('M: mat_scale_row_col_norm (for rank-deficient)');
disp('N: mat_scale_element_log (for rank-deficient)');
disp('O: pinv_rrqr (for rank-deficient)');
disp('P: ssa_compute_bin_no_null (for rank-deficient)');
disp('Q: ssa_compute_no_bin_no_null (for rank-deficient)');
disp('R: ssa_compute_exact (for rank-deficient)');
disp('S: ssa_impose_action (for rank-deficient)');
disp('____________________________________________________________');
display('0 means no error.  1 means error.');
disp('____________________________________________________________');
disp('                               ABCDEFGHIJKLMNOPQRS');

for i = 1:size(log_errors,1)
    str = sprintf('%-30s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d', H{i,1}, log_errors(i,:));
    disp(str)
end
