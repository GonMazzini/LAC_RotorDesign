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


%% Folder and file settings

%File names
path1 = 'HAWC_inputs/DTU_10MW_redesign_flexible_hawc2s'; %Main file 
path2 = 'HAWC_inputs/DTU_10MW_redesign_rigid_hawc2s'; % rigid to correct deflections 

%Operational file path:
operational_file = 'HAWC_inputs/DTU_10MW_redesign_flexible_hawc2s.opt';
%operational_file = 'HAWC_inputs/data/operation_7pt.dat';
R = 97.77; %Rotor radius

%Making operational vectors
data = readtable(operational_file, 'Filetype', 'text');
data = table2array(data);
WSP = data(:,1);
omega = data(:,3)*pi/30;
pitch = data(:,2);
TSR = round(omega*R./round(WSP,1),2);
    
%% P, T, CP, CT vs TSR

data = readtable(append(path1,'.pwr'), 'Filetype', 'text');
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
ylabel('$C_P$ [-]')
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

%% Cl, alpha Cl/Cd vs t/c (design and actual)

load('Aerodynamics\polynomials.mat');

name = {'$\alpha$ [deg]' '$C_l$ [-]', '$C_l/C_d$ [-]'};
pos = [5, 17, 18];
desired_TSR = 6.75; %Choose desired TSR value

[tc, alpha, cl, clcd] = deal(zeros(length(data(:,1)), 1));

[~, idx] = min(abs(TSR - TSRs(j)));
filename = append(path,'_u',erase(num2str(WSP(idx)),'.'),'.ind');
data = readtable(filename,'Filetype', 'text');
data = table2array(data);
r = data(:,1);
for k=1:length(data(:,1))
    t = thickness(data(k,1), polynomials);
    if t>5.38
        t = 5.38;
    end
    
    tc(k,1) = t/data(k,32)*100;

    if tc(k,1) < 24.1
        tc(k,1) = 24.1;
    end
      alpha(k,1) = alpha_des(tc(k,1), polynomials);
      cl(k,1) = cl_des(tc(k,1), polynomials);
      clcd(k,1) = clcd_des(tc(k,1), polynomials);
end



for i=1:length(pos)
    figure;
    subplot(2,1,1)
    
%     tc = zeros(length(data(:,1)),1);
    

   
    if pos(i) == 5 % alpha
        plot(tc, rad2deg(data(:,pos(i))), 'DisplayName', 'actual');
        ylabel(name(i));
        xlabel('t/c [\%]');
        hold on
        grid on
        plot(tc, alpha, 'DisplayName', 'design');
        subplot(2,1,2)
        plot(r, rad2deg(data(:,pos(i))), 'DisplayName', 'actual');
        hold on
        plot(r, alpha, 'DisplayName', 'design');
        
    elseif pos(i) == 18 % cl/cd
        plot(tc, data(:,pos(i)-1)./data(:,pos(i)), 'DisplayName', 'actual');
        ylabel(name(i));
        xlabel('t/c [\%]');
        hold on
        grid on
        plot(tc, clcd, 'DisplayName', 'design');
        subplot(2,1,2)
        plot(r, data(:,pos(i)-1)./data(:,pos(i)), 'DisplayName', 'actual');
        hold on
        plot(r, clcd, 'DisplayName', 'design');
    elseif pos(i) == 17 % cl 
        plot(tc, data(:,pos(i)), 'DisplayName', 'actual');
        ylabel(name(i));
        xlabel('t/c [\%]');
        hold on
        grid on
        plot(tc, cl, 'DisplayName', 'design');
        subplot(2,1,2)
        plot(r, data(:,pos(i)), 'DisplayName', 'actual');
        hold on
        plot(r, cl, 'DisplayName', 'design');
    end
    

    grid on
    xlabel('r [m]');
    ylabel(name(i));
    legend;
%     xlim([min(tcratio), max(tcratio) ])
end



%% a, a', alpha, Cl, Cl/Cd, CP, CT  vs r for different TSRs

name = {'a [-]', 'a''[-]', '$\alpha$ [deg]' '$C_l$ [-]', '$C_l/C_d$ [-]', '$C_P$ [-]', '$C_T$ [-]'};
pos = [2, 3, 5, 17, 18, 34, 33];
TSRs = [6.75]; %Choose desired TSR values (set equal to TSR for all)
for i=1:length(pos)
    figure;
    for j=1:length(TSRs)
    [~, idx] = min(abs(TSR - TSRs(j)));
    filename = append(path,'_u',erase(num2str(WSP(idx)),'.'),'.ind');
    data = readtable(filename,'Filetype', 'text');
    data = table2array(data);
    if pos(i) == 5
        plot(data(:,1), rad2deg(data(:,pos(i))), 'DisplayName', num2str(TSRs(j)));
    elseif pos(i) == 18
        plot(data(:,1), data(:,pos(i)-1)./data(:,pos(i)), 'DisplayName', num2str(TSRs(j)));
    else
        plot(data(:,1), data(:,pos(i)), 'DisplayName', num2str(TSRs(j)));
    end
    
    hold on
    end
    grid on
    xlabel('r [m]');
    ylabel(name(i));
    leg = legend;
    title(leg,'TSR')
    
end

%% P, T, CP, CT, pitch, omega, tip deflect. vs WSP
data_redesign_flex =  readtable(append(path1,'.pwr'),'Filetype', 'text' );
data_redesign_flex = table2array(data_redesign_flex);
data_redesign_rig =  readtable(append(path2,'.pwr'),'Filetype', 'text' );
data_redesign_rig = table2array(data_redesign_rig);

data_original_flex =  readtable('HAWC_inputs/original_design/DTU_10MW_flexible_hawc2s.pwr','Filetype', 'text' );
data_original_flex = table2array(data_original_flex);
data_original_rig =  readtable('HAWC_inputs/original_design/DTU_10MW_rigid_hawc2s.pwr','Filetype', 'text' );
data_original_rig = table2array(data_original_rig);

f1 = figure('Name', 'Power and Thrust Curves');
f2 = figure('Name', 'Twist and Rotational Speed');
f3 = figure('Name', 'Tip x');
f4 = figure('Name', 'Tip y');
f5 = figure('Name', 'Tip z');
figs = [f3,f4,f5]; 
leg = ['Redesign', 'original'];
for i=1:3
    if i==1
        data_flex = data_redesign_flex;
        data_rig = data_redesign_rig;
        linestyle = '-';
    elseif i==2
        data_flex = data_original_flex;
        data_rig = data_original_rig;
        linestyle = '--';
    end
    
    WSP = data_flex(:,1);
    P = data_flex(:,2);
    CP = data_flex(:,4);
    T = data_flex(:,3);
    CT = data_flex(:,5);
    omega = data_flex(:,9);
    pitch = data_flex(:,10);

    set(0, 'currentfigure', f1)
    subplot(2,1,1);
    yyaxis left
    plot(WSP,P, linestyle);
    ylabel('P [kW]')
    yyaxis right
    plot(WSP,CP, linestyle);
    ylabel('$C_P$ [-]')
    xlabel('$U_{\infty}$ [m/s]')
    grid on
    hold on

    subplot(2,1,2);
    yyaxis left
    plot(WSP,T);
    ylabel('T [kN]')
    yyaxis right
    plot(WSP,CT);
    ylabel('$C_T$ [-]')
    xlabel('$U_{\infty}$ [m/s]')
    grid on
    hold on

    set(0, 'currentfigure', f2)
    yyaxis left
    plot(WSP, pitch, linestyle);
    ylabel('$\theta$ [deg]')
    yyaxis right
    plot(WSP, omega, linestyle);
    ylabel('$\omega$ [rad/s]')
    xlabel('$U_{\infty}$ [m/s]')
    grid on
    hold on


    idx = [11, 12, 13];
    label = {'$x$', '$y$', '$z$'};
    for j=1:length(idx)
        deflection = data_flex(:,idx(j));% - data_rig(:,idx(j));
        set(0, 'currentfigure', figs(j))
        plot(WSP, deflection, linestyle);
        xlabel('$U_{\infty}$ [m/s]')
        ylabel(append(label(j), ' [m]'))
        grid on
        hold on
    end

end

%% Functions

function t = thickness(r,polynomials)
    %Absolute thickness [m] as a function of radius [m] for 35-m blade
    if r<=5
        r = 5;
    end
    p = polynomials.t;
    t = polyval(p,r);
end


function alpha = alpha_des(tcratio,polynomials)
    %Design AOA [deg] as a function of t/c [%]

    if tcratio<24
        tcratio = 24;
    elseif tcratio>100
        tcratio = 100;
    else
    end

    if tcratio < 48
        p = polynomials.alpha1;
        alpha = polyval(p,tcratio);        
    elseif 36 <= tcratio
        p = polynomials.alpha2;
        alpha = polyval(p,tcratio);
    else
        print('Invalid t/c input')
    end
end


function cl = cl_des(tcratio,polynomials)
    %Design cl [-] as a function of t/c [%]

    if tcratio<24
        tcratio = 24;
    elseif tcratio>100
        tcratio = 100;
    else
    end

    if tcratio < 48
        p = polynomials.cl1;
        cl = polyval(p,tcratio);
    elseif 48 <= tcratio
        p = polynomials.cl2;
        cl = polyval(p,tcratio);
    else
        print('Invalid t/c range')
    end

end

function clcd = clcd_des(tcratio,polynomials)
    %Design cl/cd [-] as a function of t/c [%]

    if tcratio<24
        tcratio = 24;
    elseif tcratio>100
        tcratio = 100;
    else
    end
    
    if tcratio < 48
        p = polynomials.clcd1;
        clcd = polyval(p,tcratio);
    elseif 48 <= tcratio
        p = polynomials.clcd2;
        clcd = polyval(p,tcratio);         
    else
        print('Invalid t/c range')
    end
end