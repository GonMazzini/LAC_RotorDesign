clear all; close all; clc

load('polynomials.mat');
load('blade_original\bladedat.txt');
olddesign.r = bladedat(:,1);
olddesign.beta = bladedat(:,2);
olddesign.c = bladedat(:,3);
olddesign.tc = bladedat(:,4);clear bladedat
olddesign.t = olddesign.tc.*olddesign.c/100;

rotor.R = 97.7;
rotor.tsr = 6.75;
rotor.B = 3;
rotor.a = 1/3;

rotor.r_lst = 5:0.5:rotor.R*0.95;
[rotor.t, rotor.c, rotor.phi, rotor.alpha, rotor.beta, rotor.cl, rotor.cd, ...
    rotor.ap, rotor.cp, rotor.ct] = deal(NaN(length(rotor.r_lst), 1));


x0 = [8.0, 0.0001]; % initial guesses
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds


for i = 1:length(rotor.r_lst)
    res = lsqnonlin(@(x)residuals(x, rotor,polynomials, i, 0, olddesign), x0, lb, ub);
    x0 = [res(1), res(2)];
    rotor = residuals(x0, rotor,polynomials, i, 1, olddesign);
end
 
rotor.c(1:5)

rotor.CP = trapz(rotor.r_lst, rotor.cp.*rotor.r_lst')*2/rotor.R^2;
rotor.CT = trapz(rotor.r_lst, rotor.ct.*rotor.r_lst')*2/rotor.R^2;
%%

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

for i=1:length(rotor.beta)
    if rotor.beta(i)>25
        rotor.beta(i) = 25;
    end
    tcratio = rotor.t(i)/rotor.c(i)*100;
    if tcratio < 24.1
        tcratio = 24.1;
        rotor.c(i) = rotor.t(i)/(tcratio/100); 
    end

    
    a = olddesign.r(1:6);
    b = olddesign.c(1:6);
    b(5) = b(5)*1.02;
    b(6) = b(6)*1.03;

    if rotor.r_lst(i)<43 && rotor.r_lst(i)>9.2
    rotor.c(i,1) = spline(a, b , rotor.r_lst(i))*rotor.R/olddesign.r(end);
    elseif rotor.r_lst(i)<10
        rotor.c(i,1) = interp1([rotor.r_lst(1), 10],...
                              [olddesign.c(1), spline(a, b , 10)*rotor.R/olddesign.r(end)], rotor.r_lst(i));
    end 

     if rotor.c(i)>rotor.R/olddesign.r(end)*max(olddesign.c)
        rotor.c(i) = rotor.R/olddesign.r(end)*max(olddesign.c);
    end
       
end


rotor.r_lst(end+1) = rotor.R;
rotor.beta(end+1) = rotor.beta(end);
rotor.t(end+1) = rotor.t(end)*0.55;
rotor.c(end+1) = rotor.c(end)*0.55;



%% Comparison of designs

figure
subplot(2,2,1)
plot(rotor.r_lst, rotor.c);hold on;
plot(olddesign.r,olddesign.c)
title('Chord distribution')
xlabel('r [m]');ylabel('c [m]')
legend('Redesigned','DTU 10 MW RWT')
grid on; box on;

subplot(2,2,2)
plot(rotor.r_lst, rotor.beta);hold on;
plot(olddesign.r,olddesign.beta)
title('Twist distribution')
xlabel('r [m]');ylabel('$\beta [deg]$')
legend('Redesigned','DTU 10 MW RWT')
grid on; box on;

subplot(2,2,3)
plot(rotor.r_lst, 100*rotor.t./rotor.c);hold on;
plot(olddesign.r,olddesign.tc)
title('Relative thickness distribution')
xlabel('r [m]');ylabel('$t/c [\%]$')
legend('Redesigned','DTU 10 MW RWT')
grid on; box on;

subplot(2,2,4)
plot(rotor.r_lst, rotor.t);hold on;
plot(olddesign.r,olddesign.t)
title('Absolute thickness distribution')
xlabel('r [m]');ylabel('t [m]')
legend('Redesigned','DTU 10 MW RWT')
grid on; box on;


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
