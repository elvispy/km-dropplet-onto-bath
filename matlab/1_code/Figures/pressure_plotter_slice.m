%% Panel (c) recreation: pressure slices at t1=0, t2=0.6, t3=2.1 (nondimensional)
clc; close all;

curr = pwd;
addpath(fullfile(curr, "..", "simulation_code"));

%% --- case folder ---
case_dir = fullfile(curr, "..", "D50Quant100", "rho1000sigma7220nu98muair0", ...
    "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U39", "N=20tol=5.00e-05");
cd(case_dir);

global errored
errored = ~isfile('z.mat');

%% --- load case-level conditions ---
try
    load('ProblemConditions.mat'); %#ok<NASGU>
catch
    load_vars('U0.mat');
    load_vars('Fr.mat');
    disp("Couldn't find ProblemConditions.mat; loaded U0/Fr instead.");
end
load_vars('vz.mat'); Vo = abs(vz(1)); %#ok<NASGU>

%% --- load fields needed for pressure reconstruction ---
load_vars('oscillation_amplitudes.mat');
load_vars('pressure_amplitudes.mat');
load_vars('numl.mat');
load_vars('tvec.mat');

%% --- load shared parameters using same folder-walk pattern ---
cd ..;  if ~isfile('Ro.mat'), cd ..; end;  load('Ro.mat','Ro'); %#ok<NASGU>
cd ..;  load('rhoS.mat','rhoS'); %#ok<NASGU>
cd ..;  load('rho.mat','rho');   %#ok<NASGU>
        load('nu.mat','nu');     %#ok<NASGU>
        load('muair.mat','muair'); %#ok<NASGU>
        load('g.mat','g');       %#ok<NASGU>
cd ..;  load('nr.mat','nr');     %#ok<NASGU>
        load('dr.mat','dr');
        load('r.mat','r');       %#ok<NASGU>
        load('zs.mat','zs');     %#ok<NASGU>
        load('xplot.mat','xplot'); %#ok<NASGU>

cd(case_dir);

%% --- targets (nondimensional times) ---
t_targets = [0.0, 0.6, 2.1];
labels_DNS = arrayfun(@(i,t) sprintf('t_{%d} = %.1f (DNS)', i, t), ...
                  1:numel(t_targets), t_targets, 'UniformOutput', false);
labels_KM = arrayfun(@(i,t) sprintf('t_{%d} = %.1f (KM)', i, t), ...
                  1:numel(t_targets), t_targets, 'UniformOutput', false);

% nearest indices in tvec
j_targets = zeros(size(t_targets));
for k = 1:numel(t_targets)
    [~, j_targets(k)] = min(abs(tvec - t_targets(k)));
end

% radial grid up to max contact extent among the 3 times
%max_numl = max(numl(j_targets));
rvec = 0:dr:1;
rc   = numl(j_targets) * dr; %#ok<NASGU>

%% --- compute pressure slices ---
P = zeros(numel(rvec), numel(t_targets));
for k = 1:numel(t_targets)
    jj = j_targets(k);

    f = zeta_generator(pressure_amplitudes(:, jj));
    thetaplots = theta_from_cylindrical(dr .* (0:numl(jj)), oscillation_amplitudes(:, jj));

    pfield_now = f(thetaplots) - sum(pressure_amplitudes(:, jj));
    pmean      = mean(f(linspace(0, pi/10, 20)) - sum(pressure_amplitudes(:, jj)));

    p_slice = (pfield_now' - pmean);
    P(1:(numl(jj)+1), k) = p_slice;
end

% OPTIONAL normalization:
% P = P ./ (rho * Vo^2);

%% --- plot: import paper fig and overlay ---
cd(curr);
figc = openfig(fullfile(curr, "..", "..", "0_data", "manual", "AirGap_Pressure.fig"));

ax = gca;
hold(ax,'on');




% Get the existing legend ONCE and keep its handle
lgd = legend(ax);

% --- LOCK AXES SIZE (pixels) and create right margin for outside legend ---
figc.Units = 'pixels';
ax.Units   = 'pixels';

axPos = ax.Position;          % fixed axes size in pixels
figPos = figc.Position;

rightMargin = 130;            % pixels reserved for legend (tune 260–380)
figPos(3) = figPos(3) + rightMargin;   % widen figure only
figc.Position = figPos;

ax.Position = axPos;          % restore axes to original pixel size
ax.PositionConstraint = 'innerposition';


% Legend children include ConstantLine (xline) + 3 DNS marker series (Line)
hAll = flipud(lgd.PlotChildren);

% Keep only DNS series entries (Line objects), exclude ConstantLine
isDNS = arrayfun(@(h) isa(h,'matlab.graphics.chart.primitive.Line'), hAll);
hDNS  = hAll(isDNS);
hDNS  = hDNS(3:-1:1);

% Relabel DNS legend entries and capture their colors
for k = 1:3
    hDNS(k).DisplayName = labels_DNS{k};
    hDNS(k).MarkerSize = k+2;   % t1 small, t2 medium, t3 large

end
dnsColor = vertcat(hDNS.Color);

% Overlay simulations using DNS colors (keep out of legend for now)
for k = 3:-1:1
    valid = ~isnan(P(:,k));
    plot(ax, rvec(valid), P(valid,k), '-', ...
        'Color', dnsColor(k,:), ...
        'LineWidth', 3*k-1, ...
        'DisplayName', labels_KM{k}, ...
        'HandleVisibility','on');
end

% --- lock axes size, then enlarge figure and place legend outside ---
figc.Units = 'pixels';
figPos = figc.Position;
figPos(3) = figPos(3) + 90;   % widen for legend
figc.Position = figPos;

% Get Contact radius handle (ConstantLine) from axes (more robust than legend children)
hContact = findobj(ax, 'Type', 'ConstantLine');
hContact = hContact(1);
hContact.DisplayName = 'Contact radius (DNS)';

% One KM legend entry (dummy handle only)
%hKM = plot(ax, NaN, NaN, '-', 'Color', [0 0 0], 'LineWidth', 3,...
%    'DisplayName', 'Kinematic Match');

% Build legend explicitly (DNS + contact + KM), outside
%lgd = legend(ax, [hDNS(:); hContact; hKM], 'Location','northeastoutside');
lgd = legend; legend(lgd.PlotChildren([1:end-3, end:-1:end-2]), lgd.String([1:end-3, end:-1:end-2]));
set(lgd, 'Interpreter','tex', 'FontSize',16, 'FontName','Cambria Math');

% Force legend into the new right margin (pixels)
lgd.Units = 'pixels';
lgd.Position(1) = axPos(1) + axPos(3) + 30;   % 20 px gap to the right of axes
lgd.Position(2) = 185;
axPos(2) = axPos(2) + 15;
drawnow;
ax.Position = axPos;  % restore axes size (prevents axes from expanding/shrinking)

% Place legend in the reserved right margin (outside the axes but inside figure)
%set(lgd, 'Units','normalized');
%lgd.Position = [0.79 0.62 0.20 0.30];   % [left bottom width height] (tune)
lgd.Box = 'on';
lgd.LineWidth = 1.2;

% --- Paper style: fonts, ticks, box, linewidths ---
ax.LineWidth = 1.6;
ax.FontName  = 'Times';
ax.FontSize  = 20;                 % tick label size (tune)

ax.TickDir   = 'in';               % matches your second figure
ax.TickLength = [0.012 0.012];     % shorter inward ticks

ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

% If you want no grid (like your second figure):
grid(ax,'off');

% Labels (roman-style like your paper)
xlabel(ax, '$r/R$', 'Interpreter','latex', ...
       'FontName','Latin Modern Roman', 'FontSize',26);
ylabel(ax, '$p/(\rho V^2)$', 'Interpreter','latex', ...
       'FontName','Latin Modern Roman', 'FontSize',26);

% Legend font to match
lgd.FontName = 'Cambria Math';
lgd.FontSize = 20;

drawnow;


box(ax,'on');
xlim(ax,[0, max(numl(j_targets)*dr)*1.25]);
ylim(ax, [-0.05, 2])

%% --- export ---
outbase = fullfile(curr, "..", "..", "0_data", "manual", "panel_c_pressure_slices");
saveas(figc, outbase, 'fig');
exportgraphics(figc, outbase + ".eps", 'ContentType','vector');

%% --- helper ---
function load_vars(str)
    global errored
    if errored, str = "errored_" + str; end
    vars = load(str);
    fn = fieldnames(vars);
    for ii = 1:numel(fn)
        assignin('caller', fn{ii}, vars.(fn{ii}));
    end
end
