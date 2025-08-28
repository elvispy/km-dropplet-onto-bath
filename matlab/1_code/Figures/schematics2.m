
close all; clear all;
addpath(fullfile(pwd, "..", "simulation_code"));
curr = pwd;
p = fullfile(pwd, "..", "D50Quant100", "rho1000sigma7220nu98muair0", "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U39", "N=20tol=5.00e-05"); % uigetdir();
cd(p);

saving = true;
% global errored
% errored = ~isfile('z.mat');
% try
% 
%     load('ProblemConditions.mat'); NN = N;
%     disp("Starting impact panel for the following parameters:");
%     fprintf("Re = %g\n", Re);
%     fprintf("We = %g\n", We);
%     fprintf("Fr = %g\n", Fr);
%     fprintf("U0 = %g cm/s\n", U0);
% catch
%     load('U0.mat');
%     load('Fr.mat');
%     disp("Couldn't find Problem Conditions");
% end
% load_vars('vz.mat'); Vo = abs(vz(1));
% 
% 
% try
%     load_vars('etas.mat');
%     etaMatPer = etas;
% catch
%     files = dir(fullfile(pwd, "etaMatPer*.mat"));
%     N = length(files);
%     etaAux = [];
%     for i = 1:N
%         load_vars(files(i).name);
%         etaAux = [etaAux, etaMatPer];
%     end
%     etaMatPer = etaAux;
% end
% load_vars('z.mat')
% load_vars('etaOri.mat')
% load_vars('tvec.mat')
% 
% load_vars('oscillation_amplitudes.mat');
% load_vars('pressure_amplitudes.mat');
% load_vars('numl.mat');
% Rv = zeros(1, size(oscillation_amplitudes, 2));
% for ii = 1:size(oscillation_amplitudes, 2)
%     Rv(ii) = zs_from_spherical(pi, oscillation_amplitudes(:, ii));
% end
% 
% %%
% 
% cd ..
% try
%     load('Ro.mat','Ro')%Sphere's radius in CGS
% catch
%     cd ..
%     load('Ro.mat','Ro')
% end
% cd ..
% load('rhoS.mat','rhoS')%Sphere density
% %load('sigmaS.mat')%Sphere's surface tension
% 
% cd ..
% load('rho.mat','rho')
% %load('sigma.mat','sigma')
% load('nu.mat','nu')
% load('muair.mat')
% load('g.mat','g') %gravitational constant
% 
% cd ..
% %load('D.mat')%Domain diameter in units of droplet radii
% %load('quant.mat')%number of dr's contained in an undeformed dropelt radius
% load('nr.mat','nr')
% load('dr.mat','dr')
% load('r.mat')
% load('zs.mat','zs')
% 
% %xplot = dr*(0:nr-1); save('xplot.mat','xplot')%I might remove or relocate this
% load('xplot.mat')
% cd(curr);
% %load('var.mat');
cd(p);


myFont = "Arial";
n=8000;
t=0:2*pi/n:2*pi;
r = ones(size(t)); 
x = r .* sin(t);
CENTER = 0.2;
y = r .* cos(t) + 1;
saving_figure = figure;
plot(x,y + CENTER,'Color', [.1 .08 .08] ,'LineWidth',4);
myax = gca; mypos = myax.Position;
axis equal
hold on
mylims = 1.95;
set(myax,'Xlim',[-mylims mylims],'Xtick',[],'Ytick',[]);
%set(myax, 'Ylim', [-1.4, -1.4+2*mylims]);

%S(t)
epsi = 1;
idx1 = floor((4-epsi)/8*n):floor((4+epsi)/8*n);
x1 = x(idx1);
y1 = y(idx1);
%plot(x1,y1-.02,'color',[.5, .5, .5], 'LineWidth', 2);%[0/256, 191/255, 255/255],'LineWidth',4)

%[xSt, ySt] = gca_to_Normalized(myax, [x1(end) + 0.02, x1(end) + 0.29], ...
%    [y1(1) - 0.28, y1(1)-0.30]);
%annotation('textarrow', xSt, ySt, 'HeadLength', 10, ...
%'HeadStyle', 'vback3', 'String', '$S(t)$', 'Interpreter', 'latex', ...
%'LineWidth', 1.5, 'FontSize', 14);
hold on


% Free surfaces (eta)
t2 = (-0.01):.01:mylims;
x2 = t2; %t2+x1(length(x1));
%y2 = y1(length(y1))+(1-exp(-t2/2.5))-.02;
ZERO = 0;

A = 0.5; k = 3; alpha = 1.5; omega = 0.5;
r = linspace(0, pi, length(x2));  % Radial distance
z =  -A * sin(k * r) .* (1-exp(-alpha * (r-0.5).^2)) .* (1-exp(-0.1 * r.^2));

cc = plot(x2,z,'color',[.5, .5, .5],'LineWidth',2);%free surfaces

%x3 = -x2;
plot(-x2,z,'color',[.5, .5, .5],'LineWidth',2)%free surfaces
hold on

% z = 0
plot([-x2(end) 0], [ZERO, ZERO], '--','color',[.8 .8 .8],'LineWidth',1);%free surfaces
aa = plot([x2(end) -0], [ZERO, ZERO], '--','color',[.8 .8 .8],'LineWidth',1);%free surfaces

% h(t)
%Centre mark

%plot([-val, val], [val+disloc, -val+disloc], 'k');
%plot([-val, val], [-val+disloc, val+disloc], 'k');

%plot([-0.02 -.15],[disloc disloc], 'Color', [.7 .7 .7], 'LineWidth', 1.1);%bottom measure transporter for h
[xht, yht] = gca_to_Normalized(myax, [sin(3*pi/4)*.99, ZERO], [(cos(3*pi/4)+1)*.99 + CENTER, ZERO+1+CENTER]);
annotation('doublearrow',xht, yht, ...
    'Head1Length', 9, 'Head2Length', 0, ...
    'Head1Style', 'vback3', 'Head2Style', 'vback3' ,'Linewidth',1 )
text(sin(3*pi/4)*0.27, .5+CENTER, "$R_d$", ...
    'Fontsize', 20, 'FontName', myFont, 'interpreter', 'latex')


% Define coordinates for the two semi-axes (using annotation)
[xht1, yht1] = gca_to_Normalized(myax, [sin(0)*.99, ZERO], [(cos(0)+1)*.99 + CENTER, ZERO+1+CENTER]);
annotation('doublearrow', xht1, yht1, 'Head1Length', 0, 'Head2Length', 0, ...
    'Head1Style', 'vback3', 'Linewidth', .75);

[xht2, yht2] = gca_to_Normalized(myax, [sin(pi/4)*.99, ZERO], [(cos(pi/4)+1)*.99 + CENTER, ZERO+1+CENTER]);
annotation('doublearrow', xht2, yht2, 'Head1Length', 7, 'Head2Length', 0, ...
    'Head1Style', 'vback3', 'Linewidth', 1);

% Plot angle arc using annotation for the angle (phi)
theta = linspace(0, pi/4, 20);
x_arc = 0.15*sin(theta);
y_arc = 0.15*cos(theta) + 1 + CENTER;
plot(x_arc, y_arc, 'k', 'LineWidth', 1.5);

% Text annotation for the angle \varphi
text(ZERO +.05, ZERO + 1 + CENTER + 0.25, '$\theta$', 'FontSize', 20, 'Interpreter', 'latex');
text(sin(pi/4)*.50, (cos(pi/4)+1)*.85, "$\xi'(t, \theta)$", 'FontSize', 16, 'Interpreter', 'latex');
%annotation('doublearrow', x_arc(end-1:end), y_arc(end-1:end), 'Head1Length', 9, 'Head2Length', 0, ...
%   'Head1Style', 'vback3', 'Linewidth', .75);
%r_c(t)
%plot([0 0],[0.07 0.16],'Color', [.7 .7 .7], 'LineWidth', 1.1);%bottom measure transporter for r_c(t)
%[xrct, yrct] = gca_to_Normalized(myax, [0.01 -x1(end)-0.02], [.1,  .1]);
%annotation('doublearrow',xrct, yrct, ...
%    'Head1Length', 0, 'Head2Length', 7, ...
%    'Head1Style', 'vback3', 'Head2Style', 'vback3' ,'Linewidth',1 );
%[Xproyrct, Yproyrct] = gca_to_Normalized(myax, [-x1(end), -x1(end)] ...
%    , [.16, ZERO+0.07]);
%annotation('arrow', Xproyrct, Yproyrct, 'Color', [.7 .7 .7], ...
%    'HeadLength', 15, 'HeadStyle', 'none', 'LineWidth', 1.5);

%text(-x1(end)/3, .2, "$r_c(t)$", ...
%    'Fontsize', 14, 'FontName', myFont, 'interpreter', 'latex');

% C(t)
%Ct = plot([x1(1)+.01, -x1(1)-0.01], [y1(1)-.05, y1(1)-.05], 'LineStyle', 'none', ...
%    'Marker', 's', 'LineWidth', 1.5, 'Color', [.2 .2 .2], ...
%    'MarkerSize', 13);

% L(t)
%Lt = plot([x1(1)+.01, -x1(1)-0.01], [ZERO, ZERO], 'LineStyle', 'none', ...
%    'Marker', '^', 'LineWidth', 1.1, 'Color', [.2 .2 .2], ...
%    'MarkerSize', 8);

% x, z axis
[xht, yht] = gca_to_Normalized(myax, [ZERO ZERO], [ZERO, 0.15]);
annotation('arrow',xht, yht, ...
    'HeadLength', 7, 'HeadStyle', 'vback3','Linewidth', 2.5)
text(.04, 0.092, "$z'$", ...
    'Fontsize', 16, 'FontName', myFont, 'interpreter', 'latex')

[xht, yht] = gca_to_Normalized(myax, [ZERO ZERO+0.15], [ZERO, ZERO]);
annotation('arrow',xht, yht, ...
    'HeadLength', 7, 'HeadStyle', 'vback3','Linewidth', 2)
text(0.04, -0.1, "$r'$", ...
    'Fontsize', 16, 'FontName', myFont, 'interpreter', 'latex')


% FLuid properties
text(-1.34, -0.2, "$(\sigma, \rho, \nu)$", ...
    'Fontsize', 20, 'FontName', myFont, 'interpreter', 'latex')
text(-0.84, 0.9, "$(\sigma_d, \rho_d, \nu_d)$", ...
    'Fontsize', 20, 'FontName', myFont, 'interpreter', 'latex')



legend([cc, aa], ["Free surface", "$z = 0$"], ...
    'FontSize', 16,  'interpreter', 'latex', 'AutoUpdate', 'off');
%legend([Ct, Lt, aa], ["$C(t)$", "$L(t)$", "$z = 0$"], ...
%    'FontSize', 11,  'interpreter', 'latex', 'AutoUpdate', 'off');
cd(curr);

if saving == true
    saveas(saving_figure, "../../0_data/manual/schematics2", 'fig');
    print(saving_figure, '-depsc', '-r300', "../../0_data/manual/schematics2.eps");
end

function [X, Y] = gca_to_Normalized(ca, xs, ys)
    % (xs(i), ys(i)) represents a point. For annotation function, xs(1),
    % ys(1) represent the initial point of the annotation and xs(2), ys(2)
    % the final point. 
    pos = ca.Position;
    xLims = ca.XLim;
    yLims = ca.YLim;
    X = zeros(size(xs));
    Y = zeros(size(ys));
    for i = 1:length(xs)
       pctX = (xs(i) - xLims(1))/(xLims(2) - xLims(1)); 
       pctY = (ys(i) - yLims(1))/(yLims(2) - yLims(1));
       
       X(i) = pos(1) + pctX * pos(3);
       Y(i) = pos(2) + pctY * pos(4);
    end
end

function load_vars(str)
    global errored
    
    if errored == true
        str = "errored_" + str; 
    end
    
    vars = load(str);
    fn = fieldnames(vars);
    for ii = 1:length(fn)
        assignin('caller', fn{ii}, vars.(fn{ii}));
    end
    
end
