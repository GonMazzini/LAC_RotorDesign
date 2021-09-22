% LAC Assignment 1
% Group 0
% Aerodynamic redesign of DTU 10MW
%% Content of this script:  
% - tip radius scaling (line 16 until 37)
% - read airfoil data and plot polars (line 39 until 88)
clear all; close all; clc

set(0,'defaulttextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',12);
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultFigureRenderer', 'painters');
set(0,'DefaultFigureWindowStyle','normal') % docked
%% Tip radius scaling

R1 = 178.3/2; % original radius
V1 = 11.4; % original rated speed
I1 = 0.16; % original turbulence intensity (class IA)
V1_max = V1*(1+2*I1); % max wind speed wit 2 stdev of TI.
I2 = 0.14; % target turbulence intensity (class IIIB)

R2_guess = 90; % new radius guess

dif = 1e12;
it = 0;
while dif>1e-6
    
    V2 = (V1^3*R1^2/R2_guess^2)^(1/3); % new rated wind speeed
    V2_max = V2*(1+2*I2); % new maximum wind speed

    R2 = R1*V1_max/V2_max; % corrected radius
    dif = abs(R2-R2_guess); % residual
    R2_guess = R2; % update radius
    it = it + 1;
end

%% Airfoil selection

% loading airfoil data
path = 'airfoil_data/';
fileinfo = dir(path);
filenames = {fileinfo.name};
filenames = filenames([fileinfo.bytes]>0);

%Airfoil models --> 24%, 30%, 36%, 48%, 60%, cylinder
polars = true; % set true to see polars
cl_des = zeros(1,length(filenames)-2); % exclude last 2 airfoils

shift = [0.38, 0.3, 0.42, 0.2]; % shift cl

for i=1:length(filenames)-2
    % open data
    data = readtable(fullfile(path,filenames{i}));
    data = table2array(data);
    % restrict range of AoA
    [~, idx_min] = min(abs(data(:,1)+15));
    [~, idx_max] = min(abs(data(:,1)-30));
    data = data(idx_min:idx_max,:);
    % find cl max in expected range of AoA
    [~, idx_min] = min(abs(data(:,1)+5));
    [~, idx_max] = min(abs(data(:,1)-15));
    [cl_max, idx] = max(data(idx_min:idx_max,2));
    % apply shift to find design  cl
    cl_des(1,i) = cl_max - shift(i);
    
    % Find interval
%     Cl_curve = data(idx_min:idx_max,2);
%     alpha_curve = data(idx_min:idx_max,1);
%     Cd_curve = data(idx_min:idx_max,3);   
    Cl_curve = data(idx_min:(idx_min+idx-1),2);
    alpha_curve = data(idx_min:(idx_min+idx-1),1);
    Cd_curve = data(idx_min:(idx_min+idx-1),3);  
    
    alpha_des(i) = interp1(Cl_curve,alpha_curve,cl_des(1,i)); 
    cd_des(i) = interp1(Cl_curve,Cd_curve,cl_des(1,i));
    
    % plotting
    if polars
        figure;
        subplot(1,2,2)
        plot(data(:,1), data(:,2), '-O');hold on; 
        yline(cl_des(1,i), '--r'); hold on;
        xline(alpha_des(i),'--r'); hold on;
        scatter(alpha_des(i),cl_des(i),'xr','linewidth',1)
        xlabel('$\alpha$ [deg]')
        ylabel('$C_l$ [-]')
        title(filenames{i}(1:end-4));
        grid on
        xlim([min(data(:,1)), max(data(:,1))])
        
%         figure;
        subplot(1,2,1)
        plot(data(:,3), data(:,2), '-O');hold on;
        yline(cl_des(1,i), '--r');hold on;
        xline(cd_des(i),'--r');hold on;
        scatter(cd_des(i),cl_des(i),'xr','linewidth',1)
        xlabel('$C_d$ [-]')
        ylabel('$C_l$ [-]')
        title(filenames{i}(1:end-4));
        grid on
    end
end

tcratio = [24, 30, 36, 48, 100];
cl_des(1,end+1) = 0;

% Fitting first 4 points with polynomial
p1_cl = polyfit(tcratio(1:end-1), cl_des(1:end-1), 3);
x1 = linspace(tcratio(1), tcratio(end-1), 40);
y1 = polyval(p1_cl,x1);

% Fitting last 2 points with straight line
p2_cl = polyfit(tcratio(end-1:end), cl_des(end-1:end), 1);
x2 = linspace(tcratio(end-1), tcratio(end), 40);
y2 = polyval(p2_cl,x2);

figure
scatter(tcratio, cl_des)
hold on
plot(x2, y2,'--k')
plot(x1, y1,'--k')
grid on
xlabel('$t/c$ [\%]')
ylabel('Design $C_l$')

% Read corresponding design alpha and design cd 
% Alpha design:
alpha_des(end+1)=0;
cd_des(end+1)=0.6;
% alpha_des = [9.3241, 9.2707, 6.0528, 3.7795, 0];
% cd_des = [0.01372, 0.01697, 0.02137, 0.03397, 0.6];
clcd_des = cl_des./cd_des;

% Fit curve to alpha design
% Fitting first 4 points with polynomial
p1_alpha = polyfit(tcratio(1:end-1), alpha_des(1:end-1), 3);
x1 = linspace(tcratio(1), tcratio(end-1), 40);
y1 = polyval(p1_alpha,x1);

% Fitting last 2 points with straight line
p2_alpha = polyfit(tcratio(end-1:end), alpha_des(end-1:end), 1);
x2 = linspace(tcratio(end-1), tcratio(end), 40);
y2 = polyval(p2_alpha,x2);

figure
scatter(tcratio, alpha_des)
hold on
plot(x2, y2,'--k')
plot(x1, y1,'--k')
grid on
xlabel('$t/c$ [\%]')
ylabel('Design $\alpha$')

% Fitting first 4 points with polynomial
p1_clcd = polyfit(tcratio(1:end-1), clcd_des(1:end-1), 3);
x1 = linspace(tcratio(1), tcratio(end-1), 40);
y1 = polyval(p1_clcd,x1);

% Fitting last 2 points with straight line
p2_clcd = polyfit(tcratio(end-1:end), clcd_des(end-1:end), 1);
x2 = linspace(tcratio(end-1), tcratio(end), 40);
y2 = polyval(p2_clcd,x2);

figure
scatter(tcratio, clcd_des)
hold on
plot(x2, y2,'--k')
plot(x1, y1,'--k')
grid on
xlabel('$t/c$ [\%]')
ylabel('Design $C_l/C_d$')

%% Thickness scaling

data = readtable('DTU_10MW_RWT_ae', 'Filetype', 'text');
data = table2array(data(:,1:4));

t_dtu = data(:,3).*data(:,2)/100;
t = t_dtu*R2/R1;
r = linspace(0, R2 ,40);

p_t = polyfit(r, t, 6);
y = polyval(p_t, linspace(0, R2, 100));

figure;
plot(data(:,1), t_dtu, 'x', 'DisplayName', 'DTU 10MW')
hold on
plot(r, t, 'x','DisplayName', 'New design')
plot(linspace(0, R2, 100), y, '--k',  'DisplayName', 'Fit')
xlabel('r [m]')
ylabel('t [m]')
grid on
legend

%% Save polynomials
curves.polynomials.t = p_t;
curves.polynomials.cl1 = p1_cl;
curves.polynomials.cl2 = p2_cl;
curves.polynomials.clcd1 = p1_clcd;
curves.polynomials.clcd2 = p2_clcd;
curves.polynomials.alpha1 = p1_alpha;
curves.polynomials.alpha2 = p2_alpha;

save('polynomials.mat','-struct','curves')