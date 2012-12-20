function X = default_ssa_compute_no_bin_no_null(A)
    X = ssa_compute(A, 0.4, inf, 0, false);
    return
