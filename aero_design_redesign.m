clear all; close all; clc

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
    res = lsqnonlin(@(x)residuals(x, rotor, i, 0), x0, lb, ub);
    x0 = [res(1), res(2)];
    rotor = residuals(x0, rotor, i, 1);
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
grid on

figure
plot(rotor.r_lst, rotor.cp)
grid on


%% F U N C T I O N S

function output = residuals(x, struct, idx, mode)

    c = x(1);
    ap = x(2);
    r = struct.r_lst(idx);

    t = thickness(r);
    tcratio = t/c*100;
    clcd = clcd_des(tcratio);
    cl = cl_des(tcratio);
    
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
        struct.alpha(idx) = alpha_des(tcratio);
        struct.beta(idx) = phi-alpha_des(tcratio);
        struct.cl(idx) = cl;
        struct.cd(idx) = cd;
        struct.ap(idx) = ap;
        struct.cp(idx) = ((1-struct.a)^2 + (struct.tsr*r/struct.R)^2*(1+ap)^2)*struct.tsr*r/struct.R*sigma*cx;
        struct.ct(idx) = ((1-struct.a)^2 + (struct.tsr*r/struct.R)^2*(1+ap)^2)*sigma*cy;
    
        output = struct;
    end
end

function t = thickness(r)
    %Absolute thickness [m] as a function of radius [m] for 35-m blade
    if r<=5
        r = 5;
    end
    t = -1.000521397083093e-10*r^6 + 3.827600673875859e-08*r^5 - 5.799931472937034e-06*r^4 ...
        + 4.275676654471356e-04*r^3 - 0.014369744235553*r^2 + 0.077853410396680*r + 5.796361265081653;
end


function alpha = alpha_des(tcratio)
    %Design AOA [deg] as a function of t/c [%]

    if tcratio<24
        tcratio = 24;
    elseif tcratio>100
        tcratio = 100;
    else
    end

    if tcratio < 48
        alpha = 0.002634259259259*tcratio^3 - 0.281034722222214*tcratio^2 ...
            + 9.382141666666378*tcratio -90.387299999996500;
    elseif 36 <= tcratio
        alpha = -0.072682692307692*tcratio + 7.268269230769232;
    else
        print('Invalid t/c input')
    end
end


function cl = cl_des(tcratio)
    %Design cl [-] as a function of t/c [%]

    if tcratio<24
        tcratio = 24;
    elseif tcratio>100
        tcratio = 100;
    else
    end

    if tcratio < 48
        cl = -1.013773148148117e-04*tcratio^3 + 0.009146458333333*tcratio^2 ...
            -0.274282499999987*tcratio + 4.146159999999838;
    elseif 48 <= tcratio
        cl = -0.016202307692308*tcratio + 1.620230769230769;
    else
        print('Invalid t/c range')
    end

end

function clcd = clcd_des(tcratio)
    %Design cl/cd [-] as a function of t/c [%]

    if tcratio<24
        tcratio = 24;
    elseif tcratio>100
        tcratio = 100;
    else
    end
    
    if tcratio < 48
        clcd = -0.002717789431521*tcratio^3 + 0.287987668161641*tcratio^2 ...
            - 13.087292916483413*tcratio + 2.900341253728703e+02;
    elseif 48 <= tcratio
        clcd = -0.476959307986685*tcratio + 47.695930798668506;
    else
        print('Invalid t/c range')
    end
end
