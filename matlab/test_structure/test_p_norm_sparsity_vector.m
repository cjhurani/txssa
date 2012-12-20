function test_p_norm_sparsity_vector

p_vals = [ 0 exp(-5:0.2:0) exp(0.0:0.2:2) inf ];
ratio_vals = 0.0:0.05:1;
n_vals = 1:30;

num_tests = 0;
num_errors = 0;

for ip = 1:length(p_vals)
    p = p_vals(ip);

    for ir = 1:length(ratio_vals)
        r = ratio_vals(ir);

        for in = 1:length(n_vals)

            a = rand(n_vals(in),1) - 0.5;
            a = [a; a; a];  % triplicate values

            for m = 0:nnz(a)
                pat = p_norm_sparsity_vector(a, r, p, m);
                num_tests = num_tests + 1;

                vec_norm = vec_p_norm(a, p);
                discarded_part_norm = vec_p_norm(a - a .* pat, p);

                have_err = 0;

                if(r == 1 && discarded_part_norm > 0)
                    display('error 1');
                    have_err = 1;
                end

                if(r < 1 && discarded_part_norm > (1 - r)*vec_norm*(1 + 2*eps))
                    display('error 2');
                    have_err = 1;
                end

                if(nnz(pat) < m)
                    display('error 3');
                    have_err = 1;
                end

                % make sure that discarding a little more violates
                % our conditions on sparsity.  This works for m = 0.
                if(m == 0 && nnz(pat) > 0)
                    [pat_ids, junk] = find(pat);
                    to_remove = pat_ids(int32(max(1,length(pat_ids)*rand())));

                    if(pat(to_remove) == 1)
                        new_pat = pat;
                        new_pat(to_remove) = 0;
                        new_discarded_part_norm = vec_p_norm(a - a .* new_pat, p);

                        if(r < 1 && new_discarded_part_norm <= (1 - r)*vec_norm)
                            display('error 4');
                            have_err = 1;
                        end
                    end
                end

                if(have_err)
                    num_errors = num_errors + 1;
                end
            end
        end
    end
end

num_tests
num_errors

