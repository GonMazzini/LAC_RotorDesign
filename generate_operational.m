TSR = 6:0.25:9;
R = 97.77;
U_inf = 8;
omega = TSR*U_inf/R*30/pi;
v = U_inf+ (1:length(TSR))/1000;


T = table(v', zeros(size(v')), omega');
writetable(T,'HAWC_inputs\data\operation_7pt.dat','Delimiter','\t', 'WriteVariableNames', false);
type HAWC_inputs\data\operation_7pt.dat