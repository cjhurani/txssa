function [num_near_zero_rows num_near_zero_cols] = near_zero_row_col(A)

assert(ndims(A) == 2, 'near_zero_row_col: A must be a 2-D matrix');

max_in_col = max(abs(A));
max_in_row = max(abs(A'));
A_max = max(max_in_col);

tol = A_max * eps(class(A)) * 100;  % MAGIC CONSTANT

num_near_zero_cols = length(find(max_in_col <= tol));
num_near_zero_rows = length(find(max_in_row <= tol));

