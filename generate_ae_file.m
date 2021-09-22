%%
radial =  rotor.r_lst'; % round(rotor.r_lst', 3, 'significant') ; %
disp(size(radial))
t_c =  100 * rotor.t ./ rotor.c;
rigid = ones(size(rotor.r_lst')); % just a column of ones

chord = rotor.c;
disp(size(chord))
disp(size(t_c))
disp(size(rigid))

% T.Properties.VariableNames = {'a',    'b',    'c',    'd'};
T = table(radial, chord, t_c, rigid);


T.semicolon = strrep(num2cell( num2str(ones(size(rotor.r_lst'))) ), '1', ';');
T.Properties.CustomProperties

writetable(T,'tabledata.dat','Delimiter','\t');
type tabledata.dat

%fid = fopen('tabledata.txt','wt');

%fgetl(fid)