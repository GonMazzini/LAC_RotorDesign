clear all
close all
clc

set(0,'defaulttextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',12);
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultFigureRenderer', 'painters');
set(0,'DefaultFigureWindowStyle','docked')

%%

data = readtable('st_new_flexible.dat', 'Filetype', 'text');
data = table2array(data);
data(1,:) = [];

data_original = readtable('st_original_flexible.dat', 'Filetype', 'text');
data_original = table2array(data_original);

[row, col] = size(data);
labels = {'$r$', '$m$', '$x_{cg}$', '$y_{cg}$', '$ri_x$', '$ri_y$', '$x_{sh}$', '$y_{sh}$', '$E$', '$G$', '$I_x$', '$I_y$', '$I_p$', '$k_x$', '$k_y$', '$A$', '$\theta$', '$x_e$', '$y_e$'};

for i=2:col
   figure;
   plot(data(:,1)/97.77, data(:,i),  'DisplayName', 'DTU 10MW Redesign')
   hold on
   grid on
   plot(data_original(:,1)/(178.3/2), data_original(:,i), 'DisplayName', 'DTU 10MW')
   xlabel('r/R [-]')
   ylabel(labels{i})
   legend
   
end

%%

[EIxx, EIyy, GJ, EA] = deal(zeros(2,length(data(:,1))));
for i=1:2
    if i==1
        df = data;
    else
        df = data_original;
    end
    EIxx(i,:) = df(:,9).*df(:,11);
    EIyy(i,:) = df(:,9).*df(:,12);
    GJ(i,:) = df(:,10).*df(:,13);
    EA(i,:) = df(:,9).*df(:,16);
end

variables = [EIxx; EIyy; GJ; EA];
label = {'$EI_{xx}$ [$Nm^2$]', '$EI_{yy}$ [$Nm^2$]', '$GJ$ [$Nm^2rad^{-1}$]', '$EA$ [N]'};
for i=1:4
    figure
    idx = (i-1)*2+1;
    var = variables(idx:idx+1,:);
    plot(data(:,1), var(1,:), 'DisplayName', 'New design');
    hold on
    plot(data_original(:,1), var(2,:), 'DisplayName', 'Original design');
    grid on
    xlabel('r [m]')
    ylabel(label(i))
    legend
end

