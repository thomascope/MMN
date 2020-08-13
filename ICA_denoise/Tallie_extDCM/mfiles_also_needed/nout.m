function value = nout(N,fcn,varargin)


if isnumeric(N)
    [value{1:N}] = fcn(varargin{:});
    value = value{N};
elseif iscell(N)
    [value{1:N{1}}] = fcn(varargin{:});
    value = value.(N{2});
elseif ischar(N)
    value = fcn(varargin{:});
    value = value.(N);
else err('unknown first input argument type for nout')
end


end

