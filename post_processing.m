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


%% 

TSR = 6:0.25:9;



fileinfo = dir('HAWC_inputs');
filenames = {fileinfo.name};
filenames = filenames([fileinfo.bytes]>0);

for i=1:length(filenames)
    if strcmp(filenames{i}(end-2:end), 'pwr')
        prefix = filenames{i}(1:end-4);
        break
    end
end



prefix = append('HAWC_inputs/', prefix);




%% P, T, CP, CT vs TSR

data = readtable(append(prefix,'.pwr'), 'Filetype', 'text');
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

load('polynomials.mat');

name = {'$\alpha$ [deg]' '$C_l$ [-]', '$C_l/C_d$ [-]'};
pos = [5, 17, 18];
desired_TSR = 6.75; %Choose desired TSR value

% tcratio = 14:100;
% for i=1:length(tcratio)
%     alpha(i,1) = alpha_des(tcratio(i), polynomials);
%     cl(i,1) = cl_des(tcratio(i), polynomials);
%     clcd(i,1) = clcd_des(tcratio(i), polynomials);
% end
[tc, alpha, cl, clcd] = deal(zeros(length(data(:,1)), 1));


idx = find(TSR == desired_TSR);
filename = sprintf(append(prefix,'_u800%d.ind'), idx-1);
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
        hold on
        plot(tc, alpha, 'DisplayName', 'design');
        subplot(2,1,2)
        plot(r, rad2deg(data(:,pos(i))), 'DisplayName', 'actual');
        hold on
        plot(r, alpha, 'DisplayName', 'design');
        
    elseif pos(i) == 18 % cl/cd
        plot(tc, data(:,pos(i)-1)./data(:,pos(i)), 'DisplayName', 'actual');
        hold on
        plot(tc, clcd, 'DisplayName', 'design');
        subplot(2,1,2)
        plot(r, data(:,pos(i)-1)./data(:,pos(i)), 'DisplayName', 'actual');
        hold on
        plot(r, clcd, 'DisplayName', 'design');
    elseif pos(i) == 17 % cl 
        plot(tc, data(:,pos(i)), 'DisplayName', 'actual');
        hold on
        plot(tc, cl, 'DisplayName', 'design');
        subplot(2,1,2)
        plot(r, data(:,pos(i)), 'DisplayName', 'actual');
        hold on
        plot(r, cl, 'DisplayName', 'design');
    end
    

    grid on
    xlabel('t/c [\%]');
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
    idx = find(TSR == TSRs(j));
    filename = sprintf(append(prefix,'_u800%d.ind'),idx-1);
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

%% P, T, CP, CT vs V
data =  readtable('HAWC_inputs/DTU_10MW_rigid_MWS_hawc2s.pwr','Filetype', 'text' );
data = table2array(data);

v =  data(:,1);
P = data(:,2);
CP = data(:,4);
T = data(:,3);
CT = data(:,5);

subplot(2,1,1);
yyaxis left
plot(v,P);
ylabel('P [kW]')
yyaxis right
plot(v,CP);
ylabel('$C_P$ [-]')
xlabel('$U_{\infty}$ [m/s]')
grid on

subplot(2,1,2);
yyaxis left
plot(v,T);
ylabel('T [kN]')
yyaxis right
plot(v,CT);
ylabel('$C_T$ [-]')
xlabel('$U_{\infty}$ [m/s]')
grid on



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