function y = collectPl(n, x)
    % Using the recurrence relation
    %(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)
    if size(x, 1) > 1; x = x'; end
    y = ones(n+1, length(x));
    y(2, :) = x;
%     if n >= 2
%         y(2, :) = my_legendre(2, x);
%     end
    for idx = 3:(n+1)
        y(idx, :) = (2 * idx - 3) / (idx-1) * x .* y(idx-1, :) - (idx - 2)/(idx-1) * y(idx-2, :);
    end
    y = y(2:end, :);
end