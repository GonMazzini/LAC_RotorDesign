%% Preparation of the htc file

%Inputs
R1 = 178.3/2-2.8; %original radius

sec = {};
for i=1:27
    sec{i,1} = append('sec ', num2str(i));
end

x_pos = [0 				
-2.06477e-05	
-0.0072881		
-0.0189235		
-0.0541282		
-0.126633		
-0.225666		
-0.288563		
-0.399194		
-0.576634		
-0.707136		
-0.791081		
-0.837195		
-0.853948		
-0.849367		
-0.79392 		
-0.716284		
-0.634358		
-0.553179		
-0.475422		
-0.40318 		
-0.330085		
-0.31014 		
-0.286719		
-0.255823		
-0.207891		
-0.089894]*(rotor.R-2.8)/R1;

y_pos = [ 7.006e-05 		
 -0.0122119		
 -0.0249251		
 -0.0273351		
 -0.0282163		
 -0.021321 		
 -0.0128378 	
 -0.00770659	
 -0.00488317	
 -0.0180296		
 -0.0501772		
 -0.0941228		
 -0.14888 		
 -0.214514		
 -0.290618		
-0.462574 		
 -0.688437		
 -0.960017		
 -1.28424		
 -1.66402		
-2.10743 		
 -2.6563 		
-2.78882 		
 -2.92517		
 -3.06577		
 -3.20952		
 -3.33685]*(rotor.R-2.8)/R1; 

z_pos = linspace(4.44089e-16, rotor.R - 2.8, 27).';

twist = -interp1(rotor.r_lst, rotor.beta, z_pos);
twist(isnan(twist)) = -25;

last = repmat(';',[27,1]);
%Table creation
T = table(sec, x_pos, y_pos, z_pos, twist, last);
[dir,~,~]=fileparts(pwd);
writetable(T,append(dir,'\HAWC_inputs/blade_coords_htc.txt'), 'Delimiter', ' ', 'QuoteStrings', false );