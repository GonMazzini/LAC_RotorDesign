clear all
close all
clc

set(0,'defaulttextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',12);
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultFigureRenderer', 'painters');
set(0,'DefaultFigureWindowStyle','normal')

%% Part 1

% data = readtable('Part1/DTU_10MW_rigid_hawc2s_u8000.ind','Filetype', 'text');
% data = table2array(data);

% name = {'a [-]', 'a''[-]', '$C_l$ [-]', '$C_d$ [-]', '$C_P$ [-]', '$C_T$ [-]', '$\phi [deg]$'};
% pos = [2, 3, 17, 18, 34, 33, 5];
% r = data(:,1);
% for i=1:length(pos)
%     figure;
%     if pos(i)==5
%         plot(r, rad2deg(data(:,pos(i))));
%     else
%     plot(r, data(:,pos(i)));
%     end
%     grid on
%     xlabel('r [m]');
%     ylabel(name(i));
% end

%% Part 2a

TSR = 6:0.5:9;
R = 99.77;
U_inf = 8;
omega = TSR*U_inf/R*30/pi;


% use this for the original values:  'HAWC_results/DTU_10MW_rigid_hawc2s.pwr'

data = readtable('HAWC_results/DTU_10MW_rigid_hawc2s.pwr', 'Filetype', 'text');
data = table2array(data);

P = data(:,2);
CP = data(:,4);
T = data(:,3);
CT = data(:,5);

subplot(2,1,1);
yyaxis left
plot(TSR,P);
ylabel('P [kW]')
yyaxis right
plot(TSR,CP);
ylabel('$C_T$ [-]')
xlabel('TSR [-]')
grid on

subplot(2,1,2);
yyaxis left
plot(TSR,T);
ylabel('T [kN]')
yyaxis right
plot(TSR,CT);
ylabel('$C_T$ [-]')
xlabel('TSR [-]')
grid on

name = {'a [-]', 'a''[-]', '$C_l$ [-]', '$C_d$ [-]', '$C_P$ [-]', '$C_T$ [-]'};
pos = [2, 3, 17, 18, 34, 33];
for i=1:length(pos)
    figure;
    for j=0:6
        
    % use this for the original values: 'HAWC_results/DTU_10MW_rigid_hawc2s_u800%d.ind'
    filename = sprintf('HAWC_results/DTU_10MW_rigid_hawc2s_u800%d.ind',j);
    data = readtable(filename,'Filetype', 'text');
    data = table2array(data);
    plot(data(:,1), data(:,pos(i)), 'DisplayName', num2str(TSR(j+1)));
    hold on
    end
    grid on
    xlabel('r [m]');
    ylabel(name(i));
    leg = legend;
    title(leg,'TSR')
end
%% for design TSR (8.0 RPM)
pos = 17;
j = 4;
filename = sprintf('HAWC_results/DTU_10MW_rigid_hawc2s_u800%d.ind',j);
data = readtable(filename,'Filetype', 'text');
data = table2array(data);
plot(data(:,1), data(:,pos), 'DisplayName', num2str(TSR(j+1)));
%load('rotor.mat')
% Side-by-side plots of the actual lift coefficient and the design lift coefficient versus relative thickness (left 
% plot) and versus radius (right plot) for design pitch and TSR


%% Part 3
data = readtable('Part1/data/DTU_10MW_RWT_ae', 'Filetype', 'text');
data = table2array(data(:,1:4));

figure;
subplot(2,2,1);
plot(data(:,1),data(:,2));
xlabel('r [m]');
ylabel('c [m]');
xlim([0,max(data(:,1))]);
grid on;

subplot(2,2,3);
plot(data(:,1),data(:,3));
xlabel('r [m]');
ylabel('t/c [\%]');
xlim([0,max(data(:,1))]);
grid on;

data2 = readtable('Part1/DTU_10MW_rigid_hawc2s_u8000.ind','Filetype', 'text');
data2 = table2array(data2);

subplot(2,2,2);
plot(data2(:,1),-rad2deg(data2(:,13)));
xlabel('r [m]');
ylabel('$\beta$ [deg]');
xlim([0,max(data2(:,1))]);
grid on;

subplot(2,2,4);
plot(data(:,1),data(:,3)/100.*data(:,2));
xlabel('r [m]');
ylabel('t [m]');
xlim([0,max(data(:,1))]);
grid on;

saveas(gcf,'image.svg')

%% Bonus

data1 = readtable('Part1/DTU_10MW_rigid_hawc2s_u8000.ind','Filetype', 'text');
data1 = table2array(data1);

data2 = readtable('Bonus/DTU_10MW_rigid_hawc2s_u8000.ind','Filetype', 'text');
data2 = table2array(data2);

name = {'a [-]', 'a''[-]', '$C_l$ [-]', '$C_d$ [-]', '$C_P$ [-]', '$C_T$ [-]'};
pos = [2, 3, 17, 18, 34, 33];
r = data1(:,1);
for i=1:length(pos)
    figure;
    plot(r, data1(:,pos(i)), 'DisplayName', 'No blade deflection');
    hold on
    plot(r, data2(:,pos(i)), 'DisplayName', 'Blade deflection');
    grid on
    xlabel('r [m]');
    ylabel(name(i));
    legend;
end

