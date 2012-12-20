function new_bin_ids = bin_mapping(bin_ids)

% new_bin_ids = bin_mapping(bin_ids)
%
% This method reorders bins so that there are no 'holes' and all new_bin_ids
% are contiguous.
%
% bin_ids must be an int32 based array.  We don't assume that it is 0-based.
% But we assume that the minimum value in it should not be the minimum int32.
% new_bin_ids is also int32 based and 1-based.
%
% For example:
%  a =              5     1    -1     1    -3    -2     7     7     3    -3
%  bin_mapping(a)   1     2     3     2     4     5     6     6     7     4

assert(isvector(bin_ids), 'bin_mapping: input array must be a 1-D vector');
assert(isequal(class(bin_ids), 'int32'), 'bin_mapping: input array must be int32');

impossible_id = intmin('int32');
min_id = min(bin_ids);

assert(impossible_id ~= min_id, 'bin_mapping: input array must not contain the minimum int');

max_id = max(bin_ids);

num_bins = max_id - min_id + 1;
work_array = impossible_id * ones(num_bins, 1, 'int32');
n_unique_vals = 0;

n = length(bin_ids);

for i = 1:n
    id = work_array(bin_ids(i) - min_id + 1);

    if(id == impossible_id)
        n_unique_vals = n_unique_vals + 1;
        work_array(bin_ids(i) - min_id + 1) = n_unique_vals;
    end
end

new_bin_ids = bin_ids;

for i = 1:n
    new_bin_ids(i) = work_array(bin_ids(i) - min_id + 1);
end

