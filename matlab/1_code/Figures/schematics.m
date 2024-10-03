
close all;
addpath(fullfile(pwd, "..", "simulation_code"));
curr = pwd;
p = fullfile(pwd, "..", "D50Quant100", "rho1000sigma7220nu98muair0", "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U39", "N=20tol=5.00e-05"); % uigetdir();
cd(p);

global errored
errored = ~isfile('z.mat');
try

    load('ProblemConditions.mat'); NN = N;
    disp("Starting impact panel for the following parameters:");
    fprintf("Re = %g\n", Re);
    fprintf("We = %g\n", We);
    fprintf("Fr = %g\n", Fr);
    fprintf("U0 = %g cm/s\n", U0);
catch
    load('U0.mat');
    load('Fr.mat');
    disp("Couldn't find Problem Conditions");
end
load_vars('vz.mat'); Vo = abs(vz(1));


try
    load_vars('etas.mat');
    etaMatPer = etas;
catch
    files = dir(fullfile(pwd, "etaMatPer*.mat"));
    N = length(files);
    etaAux = [];
    for i = 1:N
        load_vars(files(i).name);
        etaAux = [etaAux, etaMatPer];
    end
    etaMatPer = etaAux;
end
load_vars('z.mat')
load_vars('etaOri.mat')
load_vars('tvec.mat')

load_vars('oscillation_amplitudes.mat');
load_vars('pressure_amplitudes.mat');
load_vars('numl.mat');
Rv = zeros(1, size(oscillation_amplitudes, 2));
for ii = 1:size(oscillation_amplitudes, 2)
    Rv(ii) = zs_from_spherical(pi, oscillation_amplitudes(:, ii));
end

%%

cd ..
try
    load('Ro.mat','Ro')%Sphere's radius in CGS
catch
    cd ..
    load('Ro.mat','Ro')
end
cd ..
load('rhoS.mat','rhoS')%Sphere density
%load('sigmaS.mat')%Sphere's surface tension

cd ..
load('rho.mat','rho')
%load('sigma.mat','sigma')
load('nu.mat','nu')
load('muair.mat')
load('g.mat','g') %gravitational constant

cd ..
%load('D.mat')%Domain diameter in units of droplet radii
%load('quant.mat')%number of dr's contained in an undeformed dropelt radius
load('nr.mat','nr')
load('dr.mat','dr')
load('r.mat')
load('zs.mat','zs')
%load('Deload('zs.mat','zs')lta.mat','Delta')
%load('IntMat.mat','IntMat')
%load(sprintf('DTNnew345nr%dD%drefp10.mat', nr, D),'DTNnew345')
%DTN = DTNnew345;
%clear DTNnew345

%xplot = dr*(0:nr-1); save('xplot.mat','xplot')%I might remove or relocate this
load('xplot.mat')
cd(curr);
%load('var.mat');
cd(p);




myFont = "Arial";
n=8000;
t=0:2*pi/n:2*pi;
r = (1 + sum(aa .* collectPl(length(aa), cos(t)), 1));
x = r .* sin(t);
disloc = -0.25;
y = r .* cos(t) + disloc;
saving_figure = figure;
plot(x,y,'k','LineWidth',4);
myax = gca; mypos = myax.Position;
axis equal
hold on
mylims = 1.6;
set(myax,'Xlim',[-mylims mylims],'Xtick',[],'Ytick',[]);
%set(myax, 'Ylim', [-1.4, -1.4+2*mylims]);

%S(t)
%t1 = 5*pi/4:2*pi/n:7*pi/4;
epsi = 1.4;
idx1 = floor((4-epsi)/8*n):floor((4+epsi)/8*n);
x1 = x(idx1);
y1 = y(idx1);
plot(x1,y1-.08,'color',[.5, .5, .5], 'LineWidth', 5);%[0/256, 191/255, 255/255],'LineWidth',4)

%[xSt, ySt] = gca_to_Normalized(myax, [x1(1) - 0.25, x1(1) + 0.02], ...
%    [y1(1) + 0.02, y1(1)-0.12]);
[xSt, ySt] = gca_to_Normalized(myax, [x1(end) + 0.02, x1(end) + 0.29], ...
    [y1(1) - 0.28, y1(1)-0.30]);
annotation('textarrow', xSt, ySt, 'HeadLength', 10, ...
'HeadStyle', 'vback3', 'String', '$S(t)$', 'Interpreter', 'latex', ...
'LineWidth', 1.5, 'FontSize', 14);
hold on


% Free surfaces (eta)
t2 = (-0.01):.01:20;
x2 = t2+x1(1); %t2+x1(length(x1));
y2 = y1(length(y1))+(1-exp(-t2))-.075;
ZERO = 0;
CENTRE = 0;
plot(x2,y2,'color',[.5, .5, .5],'LineWidth',2)%free surfaces

x3 = -x2;
%y3 = y1(length(y1))+(1-exp(-t2))-.05;
plot(x3,y2,'color',[.5, .5, .5],'LineWidth',2)%free surfaces
hold on


% z = 0
plot([-x2(end) x1(1)], [ZERO, ZERO], '--','color',[.8 .8 .8],'LineWidth',1);%free surfaces
aa = plot([x2(end) -x1(1)], [ZERO, ZERO], '--','color',[.8 .8 .8],'LineWidth',1);%free surfaces


% A(t)
plot([x1(1),-x1(1)],[ZERO ZERO],'-.','color',[.3 .3 .3],'LineWidth',4);%'color',[0 0 0],'LineWidth',4)
[xAt, yAt] = gca_to_Normalized(myax, [x1(end)+.3, x1(end) + .5] ...
    , [ZERO - 0.15, ZERO-0.05]);
annotation('textarrow', xAt, yAt, ...
    'HeadLength', 10, 'HeadStyle', 'vback3', 'String', '$A(t)$', ...
    'Interpreter', 'latex', 'LineWidth', 1.5, 'FontSize', 14);

[leftxAt, leftyAt] = gca_to_Normalized(myax, [x1(1)+.01, x1(1)+.01] ...
    , [y1(1) + .05, ZERO-0.04]);
annotation('arrow', leftxAt, leftyAt, 'Color', [.7 .7 .7], ...
    'HeadLength', 15, 'HeadStyle', 'none', 'LineWidth', 1.5);
[rightxAt, rightyAt] = gca_to_Normalized(myax, [-x1(1)-.01, -x1(1)-.01] ...
    , [y1(1) + .05, ZERO-0.04]);
annotation('arrow', rightxAt, rightyAt, 'Color', [.7 .7 .7], ...
    'HeadLength', 15, 'HeadStyle', 'none', 'LineWidth', 1.5);

% h(t)
%Centre mark
val = 1/30;
bb = scatter(0, disloc, 20, 'k', "filled"); 
%plot([-val, val], [val+disloc, -val+disloc], 'k');
%plot([-val, val], [-val+disloc, val+disloc], 'k');

plot([-0.02 -.15],[disloc disloc], 'Color', [.7 .7 .7], 'LineWidth', 1.1);%bottom measure transporter for h
[xht, yht] = gca_to_Normalized(myax, [-.13 -.13], [disloc, ZERO]);
annotation('doublearrow',xht, yht, ...
    'Head1Length', 5, 'Head2Length', 0, ...
    'Head1Style', 'vback3', 'Head2Style', 'vback3' ,'Linewidth',1 )
text(-.4, disloc*.6, "$h(t)$", ...
    'Fontsize', 14, 'FontName', myFont, 'interpreter', 'latex')

%r_c(t)
plot([0 0],[0.07 0.16],'Color', [.7 .7 .7], 'LineWidth', 1.1);%bottom measure transporter for r_c(t)
[xrct, yrct] = gca_to_Normalized(myax, [0.01 -x1(end)-0.02], [.1,  .1]);
annotation('doublearrow',xrct, yrct, ...
    'Head1Length', 0, 'Head2Length', 7, ...
    'Head1Style', 'vback3', 'Head2Style', 'vback3' ,'Linewidth',1 );
[Xproyrct, Yproyrct] = gca_to_Normalized(myax, [-x1(end), -x1(end)] ...
    , [.16, ZERO+0.07]);
annotation('arrow', Xproyrct, Yproyrct, 'Color', [.7 .7 .7], ...
    'HeadLength', 15, 'HeadStyle', 'none', 'LineWidth', 1.5);

text(-x1(end)/3, .2, "$r_c(t)$", ...
    'Fontsize', 14, 'FontName', myFont, 'interpreter', 'latex');

% C(t)
%Ct = plot([x1(1)+.01, -x1(1)-0.01], [y1(1)-.05, y1(1)-.05], 'LineStyle', 'none', ...
%    'Marker', 's', 'LineWidth', 1.5, 'Color', [.2 .2 .2], ...
%    'MarkerSize', 13);

% L(t)
%Lt = plot([x1(1)+.01, -x1(1)-0.01], [ZERO, ZERO], 'LineStyle', 'none', ...
%    'Marker', '^', 'LineWidth', 1.1, 'Color', [.2 .2 .2], ...
%    'MarkerSize', 8);

% x', z' axis
[xht, yht] = gca_to_Normalized(myax, [0 0], [disloc, disloc/2]);
annotation('arrow',xht, yht, ...
    'HeadLength', 3, 'HeadStyle', 'vback3','Linewidth', 2)
text(.04, disloc*0.4, "$z'$", ...
    'Fontsize', 14, 'FontName', myFont, 'interpreter', 'latex')

[xht, yht] = gca_to_Normalized(myax, [ZERO ZERO-disloc/2], [disloc, disloc]);
annotation('arrow',xht, yht, ...
    'HeadLength', 3, 'HeadStyle', 'vback3','Linewidth', 2)
text(-disloc/2, -0.34, "$x'$", ...
    'Fontsize', 14, 'FontName', myFont, 'interpreter', 'latex')




legend([aa, bb], ["$z = 0$", "Centre of mass"], ...
    'FontSize', 14,  'interpreter', 'latex', 'AutoUpdate', 'off');
%legend([Ct, Lt, aa], ["$C(t)$", "$L(t)$", "$z = 0$"], ...
%    'FontSize', 11,  'interpreter', 'latex', 'AutoUpdate', 'off');
cd(curr);

saveas(saving_figure, "../../0_data/manual/schematics", 'fig');
print(saving_figure, '-depsc', '-r300', "../../0_data/manual/schematics.eps");


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
