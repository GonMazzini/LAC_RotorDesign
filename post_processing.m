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
path1 = 'HAWC_inputs/DTU_10MW_final_ASIER_flexible_hawc2s'; % Main file 
path2 = 'HAWC_inputs/DTU_10MW_final_ASIER_flexible_hawc2s'; % rigid to correct deflections 

%Operational file path:
operational_file = 'HAWC_inputs/data/operation_final.opt';
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

f1 = figure('Name', 'Power Curve');
f2 = figure('Name', 'Thrust Curves');
f3 = figure('Name', 'Twist');
f4 = figure('Name', 'Rotational Speed');
f5 = figure('Name', 'Tip x');
f6 = figure('Name', 'Tip y');
f7 = figure('Name', 'Tip z');
figs = [f5,f6,f7]; 
%leg = ['Redesign', 'original'];
for i=1:2
    if i==1
        data_flex = data_redesign_flex;
        data_rig = data_redesign_rig;
        linestyle = '-';
        colour = [0.8500 0.3250 0.0980];
        lab = 'Redesign';
        R = 97.77;
    elseif i==2
        data_flex = data_original_flex;
        data_rig = data_original_rig;
        linestyle = '--k';
        colour = 'black';
        lab = 'DTU 10MW';
        R = 178.3/2;
    end
    
    WSP = data_flex(:,1);
    P = data_flex(:,2);
    CP = data_flex(:,4);
    T = data_flex(:,3);
    CT = data_flex(:,5);
    omega = data_flex(:,10);
    pitch = data_flex(:,9);

    set(0, 'currentfigure', f1)
    plot(WSP,P, linestyle, 'Color', colour, 'DisplayName', lab);
    ylabel('P [kW]')
    xlabel('$U_{\infty}$ [m/s]')
    grid on
    hold on
    legend('Location', 'Best')
    
    set(0, 'currentfigure', f2)
    plot(WSP,T, linestyle, 'Color', colour, 'DisplayName', lab);
    ylabel('T [kN]')
    xlabel('$U_{\infty}$ [m/s]')
    grid on
    hold on
    legend
    
    set(0, 'currentfigure', f3)
    plot(WSP, pitch, linestyle, 'Color', colour, 'DisplayName', lab);
    ylabel('$\theta$ [deg]')
    grid on
    hold on
    legend('Location', 'Best')
    
    set(0, 'currentfigure', f4) 
    plot(WSP, omega, linestyle, 'Color', colour, 'DisplayName', lab);
    ylabel('$\omega$ [rad/s]')
    xlabel('$U_{\infty}$ [m/s]')
    grid on
    hold on
    legend

    idx = [11, 12, 13];
    label = {'$\Delta x$ [m]', '$y$ [m]', '$\Delta z$ [m]'};
    for j=1:length(idx)
        deflection = data_flex(:,idx(j));
        if idx(j) == 13
            deflection = data_flex(:,idx(j)) - data_rig(:,idx(j));
        elseif idx(j) == 12
            deflection = data_flex(:,idx(j)) - data_rig(:,idx(j));
        elseif idx(j) == 11
            deflection = data_flex(:,idx(j));
        end
        set(0, 'currentfigure', figs(j))
        plot(WSP, deflection, linestyle, 'Color', colour, 'DisplayName', lab);
        xlabel('$U_{\infty}$ [m/s]')
        ylabel(label(j))
        grid on
        hold on
        legend
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