clear all; close all; clc

load('polynomials.mat');
load('bladedat.txt');
olddesign.r = bladedat(:,1);
olddesign.beta = bladedat(:,2);
olddesign.c = bladedat(:,3);
olddesign.tc = bladedat(:,4);clear bladedat
olddesign.t = olddesign.tc.*olddesign.c/100;

rotor.R = 97.7;
rotor.tsr = 8.0;
rotor.B = 3;
rotor.a = 1/3;

rotor.r_lst = 5:0.5:rotor.R*0.98;
[rotor.t, rotor.c, rotor.phi, rotor.alpha, rotor.beta, rotor.cl, rotor.cd, ...
    rotor.ap, rotor.cp, rotor.ct] = deal(NaN(length(rotor.r_lst), 1));


x0 = [8, 0.0001]; % initial guesses
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds

for i = 1:length(rotor.r_lst)
    res = lsqnonlin(@(x)residuals(x, rotor,polynomials, i, 0), x0, lb, ub);
    x0 = [res(1), res(2)];
    rotor = residuals(x0, rotor,polynomials, i, 1);
end
 
rotor.c(1:5)

rotor.CP = trapz(rotor.r_lst, rotor.cp.*rotor.r_lst')*2/rotor.R^2;
rotor.CT = trapz(rotor.r_lst, rotor.ct.*rotor.r_lst')*2/rotor.R^2;
%%
% x = linspace(0,100,100);
% y = zeros(length(x),1);
% for i=1:length(x)
%     y(i,1) = cl_des(x(i));
% end
% figure
% plot(x,y)
% grid on

figure
plot(rotor.r_lst, rotor.c)
title('Chord distribution along blade length')
xlabel('r [m]')
ylabel('Chord [m]')
grid on

figure
plot(rotor.r_lst, rotor.cp)
title('Power Coefficient along blade length')
xlabel('r [m]')
ylabel('$C_p$')
grid on

%% Comparison of designs

figure
subplot(2,2,1)
plot(rotor.r_lst, rotor.c);hold on;
plot(olddesign.r,olddesign.c)
title('Chord distribution')
xlabel('r [m]');ylabel('c [m]')
legend('Redesigned','Old')
grid on; box on;

subplot(2,2,2)
plot(rotor.r_lst, rotor.beta);hold on;
plot(olddesign.r,olddesign.beta)
title('Twist distribution')
xlabel('r [m]');ylabel('$\beta [deg]$')
legend('Redesigned','Old')
grid on; box on;

subplot(2,2,3)
plot(rotor.r_lst, 100*rotor.t./rotor.c);hold on;
plot(olddesign.r,olddesign.tc)
title('Relative thickness distribution')
xlabel('r [m]');ylabel('$t/c [\%]$')
legend('Redesigned','Old')
grid on; box on;

subplot(2,2,4)
plot(rotor.r_lst, rotor.t);hold on;
plot(olddesign.r,olddesign.t)
title('Absolute thickness distribution')
xlabel('r [m]');ylabel('t [m]')
legend('Redesigned','Old')
grid on; box on;

%% F U N C T I O N S

function output = residuals(x, struct,struct2, idx, mode)

    c = x(1);
    ap = x(2);
    r = struct.r_lst(idx);

    t = thickness(r,struct2);
    tcratio = t/c*100;
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
