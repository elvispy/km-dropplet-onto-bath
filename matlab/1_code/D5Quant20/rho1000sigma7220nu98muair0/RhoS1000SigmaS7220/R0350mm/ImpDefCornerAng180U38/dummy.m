

N = 1;

angles = [3*pi/4, pi];
values = [0, 1];
f = @(thetas) ff(thetas, angles(1), angles(2));
custom = custom_project_amplitudes(angles, values, N, NaN, NaN);
old    = project_amplitudes(f, N, angles(1:2), NaN, 1);
M = 1;
fprintf("%gth custom result %.5f \n", M, custom(M));
fprintf("%gth old    result %.5f \n", M, old(M));
fprintf("%gth manual result %.5f \n", M, manual(angles, values, M));

disp(norm(old - custom));

%syms y;
G = @(n, a) (legendreP(n+1, a) - (n>0) * legendreP(n-(n > 0), a))/(2*n+1);
intPab = @(n, a, b) G(n, b) - G(n, a);

function res = ff(thetas, a, b)
    res = (cos(thetas) - cos(a))/(cos(b) - cos(a));
    res = res .* (thetas >= a) .* (thetas <= b);
end

function coef = manual(angs, values, N)
    f = @(x) interp1(cos(angs), values, x, 'linear', 0);
    coef = (2*N+1)/2 * integral(@(x) f(x) .* legendreP(N, x), -1, 1, 'AbsTol', 1e-4);
end