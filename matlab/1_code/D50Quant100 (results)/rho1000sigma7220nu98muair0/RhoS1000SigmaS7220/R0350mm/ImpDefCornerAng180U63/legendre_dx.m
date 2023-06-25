function y = legendre_dx(n, x)
    if n == 0
        y = zeros(size(x));
    else
        % The following is true
        % P^{\prime}_{n}(x) = \frac{d}{dx} \left(x P_n(x) - P_{n-1}(x) \right)
        % and, for the endpoints
        % \frac{d}{dx} P_{n+1}(x) = (n+1)P_n(x) + x \frac{d}{dx}P_{n}(x)
        y = n * (x .* my_legendre(n, x) - my_legendre(n-1, x)) ./ (x.^2 - 1);
        y(or(isnan(y), isinf(y))) = 0;
        y = y + n * (n+1) / 2 * ((x == 1) + (-1)^(n+1) * (x == -1));
        %y = y(~isnan(y)) + 
    end
end