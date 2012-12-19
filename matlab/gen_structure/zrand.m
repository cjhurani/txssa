function R = zrand(varargin)

% R = zrand
% R = zrand(m)
% R = zrand(m,n)
%
% Creates a random matrix R which has a reasonable probability of containing
% elements exactly equal to 0.5.  This is so that if 0.5 is subtracted, the
% new matrix will have some exactly zero elements.

optargin = size(varargin, 2);

assert(optargin <= 2, 'zrand: Too many optional input arguments (%d)', optargin);

if(optargin == 0 || isempty(varargin{1}))
    m = 1;
    n = 1;
elseif(optargin == 1 && numel(varargin{1}) == 1)
    m = varargin{1};
    n = varargin{1};
elseif(optargin == 2 && numel(varargin{1}) == 1 && numel(varargin{2}) == 1)
    m = varargin{1};
    n = varargin{2};
else
    assert(false, 'zrand: Bad arguments');
end

R = rand(m,n);
R(abs(R - 0.5) < 0.05) = 0.5;
