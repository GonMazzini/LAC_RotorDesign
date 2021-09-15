clear all
close all
clc

%% Tip radius

R1 = 178.3/2;
V1 = 11.4;
I1 = 0.16;
V1_max = V1*(1+2*I1);
I2 = 0.14;

R2_guess = 90;

dif = 1e12;
it = 0;
while dif>1e-6
    
    V2 = (V1^3*R1^2/R2_guess^2)^(1/3);
    V2_max = V2*(1+2*I2);

    R2 = R1*V1_max/V2_max;
    dif = abs(R2-R2_guess);
    R2_guess = R2;
    it = it + 1;
end