%% USE THIS SECTION FOR GENERATING MULTIPLE TSR WITH SAME WS
TSR = 6:0.25:9;
R = 97.77*0.95;
U_inf = 8;
omega = TSR*U_inf/R*30/pi;
v = U_inf+ (1:length(TSR))/1000;


T = table(v', zeros(size(v')), omega');
writetable(T,'HAWC_inputs\data\operation_final.dat','Delimiter','\t', 'WriteVariableNames', false);
type HAWC_inputs\data\operation_final.dat

%% USE THIS SECTION FOR GENERATING MULTIPLE SPEED POINTS
%Design TSR 7.25
%Rated WSP 10.93 m/s
%Rated omega 7.7397 RPM

TSR = 7.25;
R = 97.77;
U_inf = 4:0.5:10.5;
omega = TSR*U_inf/R*30/pi;

T = table(U_inf', zeros(size(U_inf')), omega');
writetable(T,'HAWC_inputs\data\operation_mult_WS_wPitch.dat','Delimiter','\t', 'WriteVariableNames', false);
type HAWC_inputs\data\operation_mult_WS.dat

%% USE THIS SECTION FOR CHANGING PITCH
% TODO 
