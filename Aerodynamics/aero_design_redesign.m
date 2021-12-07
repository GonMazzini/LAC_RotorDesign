clear all; close all; clc

% This code outputs a struct with the optimized design minimizing the
% residuals
% This code also provides: 

load('polynomials.mat');
% change the name and clear 
polynomialsfinal = polynomials; clear polynomials;

load('polynomials.mat');
[dir,~,~]=fileparts(pwd);
load(append(dir, '\blade_original\bladedat.txt'));
olddesign.r = bladedat(:,1);
olddesign.beta = bladedat(:,2);
olddesign.c = bladedat(:,3);
olddesign.tc = bladedat(:,4);clear bladedat
olddesign.t = olddesign.tc.*olddesign.c/100;

%% First redesign 
First_rotor.R = 97.7; % 
First_rotor.tsr = 6.75; % 
First_rotor.B = 3;
First_rotor.a = 1/3;

First_rotor.r_lst = 5:0.5:First_rotor.R;
[First_rotor.t, First_rotor.c, First_rotor.phi, First_rotor.alpha, First_rotor.beta, First_rotor.cl, First_rotor.cd, ...
    First_rotor.ap, First_rotor.cp, First_rotor.ct] = deal(NaN(length(First_rotor.r_lst), 1));


x0 = [8.0, 0.0001]; % initial guesses
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds


for i = 1:length(First_rotor.r_lst)
    res = lsqnonlin(@(x)residuals(x, First_rotor,polynomials, i, 0, olddesign), x0, lb, ub);
    x0 = [res(1), res(2)];
    First_rotor = residuals(x0, First_rotor,polynomials, i, 1, olddesign);
end
 
First_rotor.c(1:5)

First_rotor.CP = trapz(First_rotor.r_lst, First_rotor.cp.*First_rotor.r_lst')*2/First_rotor.R^2;
First_rotor.CT = trapz(First_rotor.r_lst, First_rotor.ct.*First_rotor.r_lst')*2/First_rotor.R^2;



%% Final design
loop = true; % true
if loop
    tsr_vec = 6:0.25:8;
    cp_vec = zeros(1,length(tsr_vec));
    ct_vec = zeros(1,length(tsr_vec));
else    
    tsr_vec = 7.25;
end 

for J=1:length(tsr_vec)
    rotor.R = 97.7; % 92.8831
    rotor.tsr = tsr_vec(J); % 6.75
    rotor.B = 3;
    rotor.a = 1/3;
    rotor.r_lst = 5:0.5:rotor.R;
    [rotor.t, rotor.c, rotor.phi, rotor.alpha, rotor.beta, rotor.cl, rotor.cd, ...
        rotor.ap, rotor.cp, rotor.ct] = deal(NaN(length(rotor.r_lst), 1));

    x0 = [8.0, 0.0001]; % initial guesses
    lb = [0, 0]; % lower bounds
    ub = [inf, 1]; % upper bounds


    for i = 1:length(rotor.r_lst)
        res = lsqnonlin(@(x)residuals(x, rotor, polynomialsfinal, i, 0, olddesign), x0, lb, ub);
        x0 = [res(1), res(2)];
        rotor = residuals(x0, rotor,polynomialsfinal, i, 1, olddesign);
    end

    rotor.c(1:5)

    rotor.CP = trapz(rotor.r_lst, rotor.cp.*rotor.r_lst')*2/rotor.R^2;
    rotor.CT = trapz(rotor.r_lst, rotor.ct.*rotor.r_lst')*2/rotor.R^2;
    cp_vec(1,J) = rotor.CP;
    ct_vec(1,J) = rotor.CT;
end

if loop
    figure(1)
    plot(tsr_vec  , cp_vec )
    [argvalue, argmax] = max(cp_vec);
    xline(tsr_vec(argmax), 'LineWidth', 2, 'LineStyle','--')
    grid on;
    xlabel('TSR')
    ylabel('Cp (using residuals)')
else   
    disp(rotor.tsr)
    disp(rotor.CP)
    disp(rotor.CT)
    disp("------")
end
%% Figure 1: Chord distribution

figure(2)
plot(olddesign.r,olddesign.c ,rotor.r_lst, rotor.c, First_rotor.r_lst, First_rotor.c)
legend('DTU-10MW RWT','Final Rotor. Des.', '1st Rotor. Des.')
title('Chord distribution along blade length')
xlabel('r [m]')
ylabel('Chord [m]')
grid on
rotor.c_raw = rotor.c;
% plot(rotor.r_lst, rotor.c); hold on;
% figure(2)
% plot(First_rotor.r_lst, First_rotor.c)
% title('Chord distribution along blade length')
% xlabel('r [m]')
% ylabel('Chord [m]')
% grid on

%% Figure 2: Cp (WARNING: IF LOOP=True, the last Cp will be consider and not the optimal)

figure(3)
plot(rotor.r_lst, rotor.cp ,First_rotor.r_lst, First_rotor.cp)
legend('Final Rotor. Des.', '1st Rotor. Des.')
title('Power Coefficient along blade length')
xlabel('r [m]')
ylabel('$C_p$')
grid on

%% Adjust Chord for final design 
for i=1:length(rotor.beta)
    % this is limitng twist
    if rotor.beta(i)>25
        rotor.beta(i) = 25;
    end
    tcratio = rotor.t(i)/rotor.c(i)*100;
    % limiting the thickes airfoi, this deviate from original shape near tip
    if tcratio < 24.1
        tcratio = 24.1;
        rotor.c(i) = rotor.t(i)/(tcratio/100); 
    end
    
    
    a = olddesign.r(1:6);
    b = olddesign.c(1:6);
    b(5) = b(5)*1.02;  %1.02
    b(6) = b(6)*1.03;  %1.03
    
    if rotor.r_lst(i)<32 && rotor.r_lst(i)>9.2
    rotor.c(i,1) = spline(a, b , rotor.r_lst(i))*rotor.R/olddesign.r(end);
    elseif rotor.r_lst(i)<10 % linear interp
        rotor.c(i,1) = interp1([rotor.r_lst(1), 10],...
                              [olddesign.c(1), spline(a, b , 10)*rotor.R/olddesign.r(end)], rotor.r_lst(i));
    end 
    
    % lo siguiente hace un ajuste continuo entre spline y ultimo tramo
    if rotor.c(i)>rotor.R/olddesign.r(end)*max(olddesign.c)
        rotor.c(i) = rotor.R/olddesign.r(end)*max(olddesign.c);
    end
       
end


rotor.r_lst(end+1) = rotor.R;
rotor.beta(end+1) = rotor.beta(end);
rotor.t(end+1) = rotor.t(end)*0.55;
rotor.c(end+1) = rotor.c(end)*0.55;

rotor.c(1:end/2) = smooth(rotor.c(1:end/2));

%% Adjust chord for First design 
for i=1:length(First_rotor.beta)
    % this is limitng twist
    if First_rotor.beta(i)>25
        First_rotor.beta(i) = 25;
    end
    tcratio = First_rotor.t(i)/First_rotor.c(i)*100;
    % limiting to the thickes airfoil
    if tcratio < 24.1
        tcratio = 24.1;
        First_rotor.c(i) = First_rotor.t(i)/(tcratio/100); 
    end
    
    
    a = olddesign.r(1:6);
    b = olddesign.c(1:6);
    b(5) = b(5)*1.02;  %1.02
    b(6) = b(6)*1.03;  %1.03
    
    if First_rotor.r_lst(i)<43 && First_rotor.r_lst(i)>9.2
    First_rotor.c(i,1) = spline(a, b , First_rotor.r_lst(i))*First_rotor.R/olddesign.r(end);
    elseif First_rotor.r_lst(i)<10 % linear interp
        First_rotor.c(i,1) = interp1([First_rotor.r_lst(1), 10],...
                              [olddesign.c(1), spline(a, b , 10)*First_rotor.R/olddesign.r(end)], First_rotor.r_lst(i));
    end 

    if First_rotor.c(i)>First_rotor.R/olddesign.r(end)*max(olddesign.c) % entra con r = 23 m, i=37
        First_rotor.c(i) = First_rotor.R/olddesign.r(end)*max(olddesign.c);
    end
       
end


First_rotor.r_lst(end+1) = First_rotor.R;
First_rotor.beta(end+1) = First_rotor.beta(end);
First_rotor.t(end+1) = First_rotor.t(end)*0.55;
First_rotor.c(end+1) = First_rotor.c(end)*0.55;


%% Figure 3 (subplots) - Comparison of designs
% compare raw chord vs modified
figure(4)
plot(rotor.r_lst(1:end-1), rotor.c(1:end-1),rotor.r_lst(1:end-1), rotor.c_raw)
xline(9.2)
xline(32)
ylabel('chord'); xlabel('radio [m]');legend('Adjusted chord','Raw design');
%% Subplots for DTU, 1st Redesign and Final design
figure(5)
subplot(2,2,1)
plot(olddesign.r,olddesign.c);hold on; % DTU
plot(First_rotor.r_lst,First_rotor.c)
plot(rotor.r_lst, rotor.c);hold on; % final
title('Chord distribution')
xlabel('r [m]');ylabel('c [m]')
legend('DTU 10 MW RWT','1st Design','Final Design')
grid on; box on;

subplot(2,2,2) % beta
plot(olddesign.r,olddesign.beta);hold on; % DTU
plot(First_rotor.r_lst,First_rotor.beta);hold on; 
plot(rotor.r_lst, rotor.beta);hold on; % final
title('Twist distribution')
xlabel('r [m]');ylabel('$\beta [deg]$')
legend('DTU 10 MW RWT','1st Design','Final Design')
grid on; box on;

subplot(2,2,3) % Relative thicknes
plot(olddesign.r,olddesign.tc);hold on; % DTU
plot(First_rotor.r_lst, 100*First_rotor.t./First_rotor.c) 
plot(rotor.r_lst, 100*rotor.t./rotor.c);hold on; % final
title('Relative thickness distribution')
xlabel('r [m]');ylabel('$t/c [\%]$')
legend('DTU 10 MW RWT','1st Design','Final Design')
grid on; box on;

subplot(2,2,4) % Absolute thickness 
plot(olddesign.r,olddesign.t);hold on; % DTU
plot(First_rotor.r_lst, First_rotor.t);
plot(rotor.r_lst, rotor.t);hold on; % final
title('Absolute thickness distribution')
xlabel('r [m]');ylabel('t [m]')
legend('DTU 10 MW RWT','1st Design','Final Design')
grid on; box on;

%% Generate geometry file

%Uncomment to generate file with radius, chord, thickness, twist and t/c.

% R  = rotor.r_lst';
% c = rotor.c;
% t = rotor.t;
% beta = rotor.beta;
% tc = 100*rotor.t./rotor.c;
% 
% T = table(R, c, t, beta, tc, 'VariableNames', { 'R', 'c', 't', 'beta', 'tc'});
% 
% writetable(T, 'redesign_geometry.txt')


%%



%% F U N C T I O N S

function output = residuals(x, struct, struct2, idx, mode, olddesign)
    r = struct.r_lst(idx);

    c = x(1);
    ap = x(2);
    
    t = thickness(r,struct2);
    if t>max(olddesign.t)
        t = max(olddesign.t);
    end
    tcratio = t/c*100;
%     if tcratio < 24.1
%         tcratio = 24.1;
%         c = t/(tcratio/100); 
%         %t = c* (tcratio/100);
%     end
    clcd = clcd_des(tcratio,struct2);
    cl = cl_des(tcratio,struct2);
    
    % intermediate variables
    phi = atan((1-struct.a)*struct.R/((1+ap)*r*struct.tsr));
    cd = cl/clcd;
    cy = cl*cos(phi) + cd*sin(phi);
    cx = cl*sin(phi) - cd*cos(phi);
    f = struct.B/2*(struct.R-r)/(r*sin(phi));
    F = 2/pi*acos(exp(-f));
    sigma = struct.B*c/(2*pi*r);
    
    % residuals
    
    res_c = 4*pi*r*sin(phi)^2*F*2*struct.a/(cy*struct.B*(1-struct.a)) - c;
    res_ap = 1/(4*F*sin(phi)*cos(phi)/sigma/cx-1) - ap;
    res = [res_c, res_ap];

    if mode == 0
        output = res;
       
    else
        struct.t(idx) = t;
        struct.c(idx) = c;
        struct.phi(idx) = phi;
        struct.alpha(idx) = alpha_des(tcratio,struct2);
        struct.beta(idx) = rad2deg(phi)-alpha_des(tcratio,struct2);
        struct.cl(idx) = cl;
        struct.cd(idx) = cd;
        struct.ap(idx) = ap;
        struct.cp(idx) = ((1-struct.a)^2 + (struct.tsr*r/struct.R)^2*(1+ap)^2)*struct.tsr*r/struct.R*sigma*cx;
        struct.ct(idx) = ((1-struct.a)^2 + (struct.tsr*r/struct.R)^2*(1+ap)^2)*sigma*cy;
    
        output = struct;
    end
end

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


