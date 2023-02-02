function plot_amplitudes()
syms x
legendreP(0.5, [x, x])
clf;
syms ang;
x_coordinate = @(angle) double(subs(sin(ang) * (1 + sum(amplitudes .* legendreP(0:(length(amplitudes)-1), cos(ang)))), ang, angle));
y_coordinate = @(angle) double(subs(cos(ang) * (1 + sum(amplitudes .* legendreP(0:(length(amplitudes)-1), cos(ang)))), ang, angle));
angles = linspace(0, 2*pi, 100);
X_cord = x_coordinate(angles);
Y_cord = y_coordinate(angles);
max_idx = find(X_cord == max(X_cord), 1);
fprintf("The maximum radius possible with given amplitudes is r = %.6g\n", X_cord(max_idx));
fprintf("This maximum radius corresponds to the following angle: %.6g\n", angles(max_idx));
plot(X_cord, Y_cord);
hold on;
line([0, X_cord(max_idx)], [0, Y_cord(max_idx)]);
xline(X_cord(max_idx));
xlim([-2, 2]);
ylim([-2, 2]);
end