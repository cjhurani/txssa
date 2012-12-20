function pos = upper_bound(v, value)

% pos = upper_bound(v, value)
%
% For any real 'value' and v: each of value < v(1:upper_bound(v, value)-1) is
% false.  Additionally, upper_bound(v, value) is the largest possible integer
% that satisfies this.
% Preconditions: v must be sorted.  pos is 1-based.
% Note: all v(lower_bound(v,value):upper_bound(v,value)-1) ~= value.

assert(isvector(v),       'upper_bound: first input must be a vector');
assert(issorted(v),       'upper_bound: first input must be sorted');
assert(isreal(v),         'upper_bound: first input must be real');
assert(numel(value) == 1, 'upper_bound: second input must be a number');
assert(isreal(value),     'upper_bound: second input must be real');

n = length(v);

if(n == 0)
    pos = 0;
else
    low = 1;
    high = n;

    while(low <= high)
        mid = floor((low + high) / 2);
        if(value < v(mid))
            high = mid - 1;
        else
            low = mid + 1;
        end
    end

    pos = low;

    assert(~any(value < v(1:pos-1)), 'upper_bound: internal error 1');
    if(pos <= n)
        assert(value < v(pos), 'upper_bound: internal error 2');
    end

    v = v(:); % so concatenation works below
    assert(issorted([v(1:pos-1); value; v(pos:end)]), ...
        'upper_bound: internal error 3');
    if(pos <= n)
        assert(~issorted([v(1:pos); value; v(pos+1:end)]), ...
            'upper_bound: internal error 4');
    end
end
