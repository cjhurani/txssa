function [min_left, max_left, min_right, max_right, did_perturb] = ...
    separated_min_max(v, varargin)

% [min_left, max_left, min_right, max_right, did_perturb] = ...
%    separated_min_max(v [, separated_at] [, perturb_fuzz])
%
% Computes minimum and maximum values of a vector when it is divided
% into two by the value "separated_at".  Optionally (when 0 < perturb_fuzz),
% see if the min/max values can be perturbed so that they become symmetric
% around "separated_at".
%
% varargin can be separated_at, perturb_fuzz.
% separated_at default is 0.
% perturb_fuzz default is 0 (so no perturbation).  Use something like 1E2 or 1E3
% to avoid asymmetry arising due to floating point issues.
%
% If there are no values in v on one/both sides of separated_at, the
% corresponding min will be greater than corresponding max.

assert(isvector(v), 'separated_min_max: first input (v) must be a 1-D vector');
assert(isreal(v),   'separated_min_max: first input (v) must be real');

optargin = size(varargin, 2);
assert(optargin <= 2, 'separated_min_max: Too many optional input arguments (%d)', optargin);

if(optargin == 0 || isempty(varargin{1}))
    separated_at = 0;
else
    separated_at = varargin{1};
    assert(numel(separated_at) == 1, 'separated_min_max: second input must be a number');
    assert(isreal(separated_at),     'separated_min_max: second input must be real');
    assert(~isinf(separated_at),     'separated_min_max: second input must be finite');
end

if(optargin <= 1 || isempty(varargin{2}))
    perturb_fuzz = 0;
else
    perturb_fuzz = varargin{2};
    assert(numel(perturb_fuzz) == 1, 'separated_min_max: third input must be a number');
    assert(isreal(perturb_fuzz),     'separated_min_max: third input must be real');
    assert(0 <= perturb_fuzz,        'separated_min_max: third input (perturb_fuzz = %f) must not be less than 0', perturb_fuzz);
end

min_left  =  inf;
max_left  = -inf;
min_right =  inf;
max_right = -inf;

n = length(v);

for i = 1:n
    vi = v(i);
    if(vi == separated_at)
        max_left  = separated_at;
        min_right = separated_at;
    elseif(vi < separated_at)
        if(max_left < vi)
            max_left = vi;
        end
        if(vi < min_left)
            min_left = vi;
        end
    else % if(separated_at < vi)
        if(max_right < vi)
            max_right = vi;
        end
        if(vi < min_right)
            min_right = vi;
        end
    end
end

did_perturb = false;

if(0 < perturb_fuzz && min_left < separated_at && separated_at < max_right)
    max_left_dist  = separated_at - max_left;
    min_left_dist  = separated_at - min_left;
    max_right_dist = max_right - separated_at;
    min_right_dist = min_right - separated_at;

    large_diff = abs(max_right_dist - min_left_dist);
    small_diff = abs(min_right_dist - max_left_dist);

    rel_tol = perturb_fuzz * eps(class(v));

    large_tol = rel_tol * (max_right_dist + min_left_dist);
    small_tol = rel_tol * (min_right_dist + max_left_dist);

    if(large_diff <= large_tol && small_diff <= small_tol)
        % We have a sequence such that when viewed on both sides of
        % separated_at, the values are roughly within same min/max
        % range.  This does not mean that for every left value there
        % is one and only one right value.  But it makes sense that we
        % should produce a binning that tries to adjust somewhat for
        % this scenario even in presence of floating point errors in
        % incoming data.

        % Choose value that is closer to separated_at between max_left, min_right
        if(max_left_dist < min_right_dist)
            min_right = separated_at + max_left_dist;
        else
            max_left  = separated_at - min_right_dist;
        end

        % Choose value that is farther from separated_at between max_right, min_left
        if(min_left_dist < max_right_dist)
            min_left  = separated_at - max_right_dist;
        else
            max_right = separated_at + min_left_dist;
        end

        % If the code above works, these asserts should not fail even in
        % presence of any floating point inaccuracies.
        assert((max_right - separated_at) == (separated_at - min_left), 'separated_min_max: Internal error 1');
        assert((min_right - separated_at) == (separated_at - max_left), 'separated_min_max: Internal error 2');

        did_perturb = true;
    end
end
