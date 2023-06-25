function y = collectPl(n, x)
    % Using the recurrence relation
    %(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)
    if size(x, 1) > 1; x = x'; end
    y = zeros(n, length(x));
    y(1, :) = x;
    if n >= 2
        y(2, :) = my_legendre(2, x);
    end
    for idx = 3:n
        y(idx, :) = (2 * idx - 1) / idx * x .* y(idx-1, :) - (idx - 1)/idx * y(idx-2, :);
    end
end