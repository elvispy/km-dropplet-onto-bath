function y = my_legendre(n, x)
    if n < 0; y = zeros(length(n), length(x)); return; end
    A = legendre(n, x);
    y = A(1, :);
end