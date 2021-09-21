clear all; close all; clc

rotor.R = 97.77;
rotor.tsr = 9.0;
rotor.B = 3;
rotor.a = 1/3;

rotor.r_lst = 5:0.5:0.98*rotor.R ;
[rotor.t, rotor.c, rotor.phi, rotor.alpha, rotor.beta, rotor.cl, rotor.cd, ...
    rotor.ap, rotor.cp, rotor.ct] = deal(NaN(length(rotor.r_lst),1));


x0 = [3, 0.0001]; % initial guesses
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds

for i = 1:length(rotor.r_lst)
    res = lsqnonlin(@(x)residuals(x, rotor, i, 0),x0,lb,ub);
    x0 = [res(1), res(2)];
    rotor = residuals(x0, rotor, i, 1);
end
 
rotor.c(1:5)
% output = residuals([3.0,0.001], rotor, 1, 1)

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
    t = 9.35996E-08*r^6 - 1.2911E-05*r^5 + 7.15038E-04*r^4 - 2.03735E-02*r^3 ...
        + 3.17726E-01*r^2 - 2.65357E+00*r + 10.2616;
end


function alpha = alpha_des(tcratio)
    %Design AOA [deg] as a function of t/c [%]

    if tcratio<15
        tcratio = 15;
    elseif tcratio>100
        tcratio = 100;
    else
    end

    if tcratio < 36
        alpha = 2.32706E-05*tcratio^5 - 2.87331E-03*tcratio^4 + 1.36343E-01*tcratio^3 ...
                    -3.10470E+00*tcratio^2 + 3.38460E+01*tcratio - 1.36500E+02;
    elseif 36 <= tcratio
        alpha = -0.0078*tcratio + 0.7813;
    else
        print('Invalid t/c input')
    end
end


function cl = cl_des(tcratio)
    %Design cl [-] as a function of t/c [%]

    if tcratio<15
        tcratio = 15;
    elseif tcratio>100
        tcratio = 100;
    else
    end

    if tcratio < 36
        cl = -7.34862E-07*tcratio^5 + 1.10229E-04*tcratio^4 - 6.40432E-03*tcratio^3 ...
                    + 1.79563E-01*tcratio^2 - 2.43397E+00*tcratio + 1.36000E+01;
    elseif 36 <= tcratio
        cl = -0.0094*tcratio + 0.9375;
    else
        print('Invalid t/c range')
    end

end

function clcd = clcd_des(tcratio)
    %Design cl/cd [-] as a function of t/c [%]

    if tcratio<15
        tcratio = 15;
    elseif tcratio>100
        tcratio = 100;
    else
    end
    
    if tcratio < 36
        clcd = -8.10212E-03*tcratio^4 + 8.73313E-01*tcratio^3 - 3.41180E+01*tcratio^2 ...
                + 5.66297E+02*tcratio -3.24932E+03;
    elseif 36 <= tcratio
        clcd = -0.8906*tcratio + 89.063;
    else
        print('Invalid t/c range')
    end
end
