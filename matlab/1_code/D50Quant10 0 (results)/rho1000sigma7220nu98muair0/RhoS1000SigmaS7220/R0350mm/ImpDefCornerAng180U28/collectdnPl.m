function y = collectdnPl(nmax, x, varargin) % TODO: Change to cell of function handles created on compile time. It is more efficient

    if nargin == 2 || varargin{1} == 1
        % THe following reucrrence relation may be derived from the
        % original legendre polynomial's recurrence relation:
         % n P'_{n+1}(x) = (2n+1) * x *P'_n(x) - (n+1)* P'_{n-1}(x)
        y = ones(nmax, length(x));
        if nmax >= 2
            y(2, :) = legendre_dx(2, x);
        end
        
        for idx = 2:(nmax-1)
            y(idx + 1, :) = (2*idx + 1) / idx * x .* y(idx, :) - (idx+1)/idx * y(idx-1, :);
        end
    elseif varargin{1} == 2
        
        % Using (d/dx on the expression on the if block above, we get:
        % (n-1) P'_{n+1}(x) = (2n+1) * x *P'_n(x) - (n+2)* P'_{n-1}(x)
        % DEP : v = collectdnPl(nmax, x);
        y = zeros(nmax, length(x));
        if nmax >=2; y(2, :) = 3; end
        if nmax >=3; y(3, :) = legendre_ddx(3, x); end
        for idx = 3:(nmax-1)
           y(idx + 1, :) = (2*idx+1)/(idx-1) * x .* y(idx, :) - (idx + 2)/(idx - 1) * y(idx-1, :);
        end
    end
end