function y = legendre_ddx(n, x)
    if n<=1
       y = zeros(size(x));
    else
        % Following recurrence relation derived from
        % https://en.wikipedia.org/wiki/Legendre_polynomials#Recurrence_relations
        % 

        y = (- n * (x.^2 + 1) .* my_legendre(n, x) ...
            + 2 * n * x .* my_legendre(n-1, x))./((x.^2 - 1).^2) ...
            + n * (x.* legendre_dx(n, x) - legendre_dx(n-1, x)) ./(x.^2 - 1);
        y(or(isnan(y), isinf(y))) = 0;
        y = y + (n-1) * n * (n+1) * (n+2) / 8 .* ((x==1) + (-1)^n .* (x == -1) );
    end

end