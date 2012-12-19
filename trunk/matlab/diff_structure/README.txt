The functions in this directory take a matrix and return a number
that tells how far is the matrix far from satisfying a given property.
The numbers are absolute and not relative.  They are based on Frobenius
norm.  For example, diff_symmetric(A) := norm(A - A','fro');
Some functions depend on functions kept in the directory 'gen_structure'
and thus it must be in path.
