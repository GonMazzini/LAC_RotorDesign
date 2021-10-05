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

data = readtable('st_new.dat', 'Filetype', 'text');
data = table2array(data);
data(1,:) = [];

data_original = readtable('st_original.dat', 'Filetype', 'text');
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
