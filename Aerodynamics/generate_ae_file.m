%%
RforDAT = [0; rotor.r_lst(1:end-1)'];  % put a zero for the first row
TforDAT = [rotor.t(1) ; rotor.t(1:end-1)];
CforDAT = [rotor.c(1) ; rotor.c(1:end-1)];

radial =  RforDAT; % round(rotor.r_lst', 3, 'significant') ; %

t_c =  100 * TforDAT ./ CforDAT;

rigid = ones(size(RforDAT'))'; % just a column of ones

chord = CforDAT;


radial = [1; 1; radial];
chord = [cell(1,1); num2cell(length(t_c)); num2cell(chord)];
t_c = [cell(1,1); cell(1,1); num2cell(t_c)];
rigid = [cell(1,1); cell(1,1); num2cell(rigid)];

disp(size(radial))
disp(size(chord))
disp(size(t_c))
disp(size(rigid))

% T.Properties.VariableNames = {'a',    'b',    'c',    'd'};
T = table(radial, chord, t_c, rigid);
T.semicolon = [cell(1,1); cell(1,1) ;strrep(num2cell( num2str(ones(size(RforDAT))) ), '1', ';')];

[dir,~,~]=fileparts(pwd);
writetable(T,append(dir,'\HAWC_inputs\data\DTU_10MW_RWT_final3_ae.dat'),'Delimiter','\t', 'WriteVariableNames', false);
% type HAWC_inputs\data\DTU_10MW_RWT_final_ae.dat

%warning('Remember to add the number of lines in the _ae file');
%fid = fopen('tabledata.txt','wt');

%fgetl(fid)