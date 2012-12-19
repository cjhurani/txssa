function pos = lower_bound(v, value)

% pos = lower_bound(v, value)
%
% For any real 'value' and v: each of v(1:lower_bound(v, value)-1) < value is
% true.  Additionally, lower_bound(v, value) is the largest possible integer
% that satisfies this.
% Preconditions: v must be sorted.  pos is 1-based.
% Note: all v(lower_bound(v,value):upper_bound(v,value)-1) ~= value.

assert(isvector(v),       'lower_bound: first input must be a vector');
assert(issorted(v),       'lower_bound: first input must be sorted');
assert(isreal(v),         'lower_bound: first input must be real');
assert(numel(value) == 1, 'lower_bound: second input must be a number');
assert(isreal(value),     'lower_bound: second input must be real');

n = length(v);

if(n == 0)
    pos = 0;
else
    low = 1;
    high = n;

    while(low <= high)
        mid = floor((low + high) / 2);
        if(~(v(mid) < value))
            high = mid - 1;
        else
            low = mid + 1;
        end
    end

    pos = low;

    assert(all(v(1:pos-1) < value), 'lower_bound: internal error 1');
    if(1 < pos)
        assert(v(pos-1) < value, 'lower_bound: internal error 2');
    end

    v = v(:); % so concatenation works below
    assert(issorted([v(1:pos-1); value; v(pos:end)]), ...
        'lower_bound: internal error 3');
    if(1 < pos)
        assert(~issorted([v(1:pos-2); value; v(pos-1:end)]), ...
           'lower_bound: internal error 4');
    end
end
