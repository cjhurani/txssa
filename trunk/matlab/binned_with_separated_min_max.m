function bin_ids = binned_with_separated_min_max(v, max_num_bins, ...
    min_left, max_left, min_right, max_right, varargin)

% bin_ids = binned_with_separated_min_max(v, max_num_bins, ...
%    min_left, max_left, min_right, max_right, varargin)
%
% Computes bins for values in vector v and using min/max left/right.
% Maximum number of allowable bins is given by max_num_bins.
% If max_num_bins == 0, all bin ids are sequentially allocated and
% are unique.
% varargin can be "separated_at", which is how the separated
% min/max were computed.  If not provided, default separated_at is 0.
%
% Note: bin_ids is 0-based and is an int32 array.

assert(isvector(v), 'binned_with_separated_min_max: first input (v) must be a 1-D vector');
assert(isreal(v),   'binned_with_separated_min_max: first input (v) must be real');

assert(~isempty(max_num_bins), 'binned_with_separated_min_max: second input (max_num_bins) must have some elements');
assert(isreal(max_num_bins),   'binned_with_separated_min_max: second input (max_num_bins) must be real');
assert(0 <= max_num_bins,      'binned_with_separated_min_max: second input (max_num_bins) must be non-negative');
assert(max_num_bins == int32(max_num_bins), 'binned_with_separated_min_max: second input (max_num_bins) must be integral');

optargin = size(varargin, 2);
assert(optargin <= 1, 'binned_with_separated_min_max: Too many optional input arguments (%d)', optargin);

if(optargin == 0 || isempty(varargin{1}))
    separated_at = 0;
else
    separated_at = varargin{1};
    assert(numel(separated_at) == 1, 'separated_min_max: separated_at must be a number');
    assert(isreal(separated_at),     'separated_min_max: separated_at must be real');
    assert(~isinf(separated_at),     'separated_min_max: separated_at must be finite');
end

if(min_left <= max_left)
    assert(all(min_left <= v),       'binned_with_separated_min_max: some element of v is less than min_left = %f', min_left);
    assert(min_left <= separated_at, 'binned_with_separated_min_max: incompatible min_left = %f, separated_at = %f', min_left, separated_at);
end

if(min_right <= max_right)
    assert(all(v <= max_right),      'binned_with_separated_min_max: some element of v is greater than max_right = %f', max_right);
    assert(separated_at <= max_right, 'binned_with_separated_min_max: incompatible separated_at = %f, max_right = %f', separated_at, max_right);
end

if(max_num_bins == 0)
    bin_ids = int32(0:(length(v)-1));
    if(size(v,2) < size(v,1)) % if a column instead of a row.
        bin_ids = bin_ids';
    end
    return;
end

% separate bin for values at separated_at
if(max_left == min_right)
    loc_max_num_left_right_bins = max_num_bins - 1;
else
    loc_max_num_left_right_bins = max_num_bins;
end

% separated_at may be <= min_left in case there are no values in v on the left
if(separated_at > min_left)
    left_dist = separated_at - min_left;
else
    left_dist = 0;
end

% max_right may be <= separated_at in case there are no values in v on the right
if(max_right > separated_at)
    right_dist = max_right - separated_at;
else
    right_dist = 0;
end

total_dist = left_dist + right_dist;

max_n_left_bins  = floor((loc_max_num_left_right_bins * left_dist ) / total_dist);
max_n_right_bins = floor((loc_max_num_left_right_bins * right_dist) / total_dist);

assert(max_n_left_bins + max_n_right_bins <= loc_max_num_left_right_bins);

inv_h_l = max_n_left_bins/left_dist;
inv_h_r = max_n_right_bins/right_dist;

fuzz = 1E2; % MAGIC CONSTANT

left_tol =  fuzz * eps(class(v)) * left_dist;
right_tol = fuzz * eps(class(v)) * right_dist;

n = length(v);

bin_ids = -ones(size(v), 'int32');

middle_id = max_n_left_bins + max_n_right_bins;

for i = 1:n
    vi = v(i);

    if(vi == separated_at)
        b = middle_id;
    elseif(vi < separated_at)
        if(separated_at - vi <= left_tol)
            b = middle_id;
        else
            b = floor((vi - min_left) * inv_h_l);
        end
    else % if(separated_at < vi)
        if(vi - separated_at <= right_tol)
            b = middle_id;
        else
            b = floor((max_right - vi) * inv_h_r) + max_n_left_bins;
        end
    end

    bin_ids(i) = b;
end

assert(all(0 <= bin_ids),         'binned_with_separated_min_max: Internal error 1');
assert(all(bin_ids <= middle_id), 'binned_with_separated_min_max: Internal error 2');
assert(length(unique(bin_ids)) <= max_num_bins, 'binned_with_separated_min_max: Internal error 3');

