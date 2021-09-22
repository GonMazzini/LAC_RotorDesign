clear ans; close all; clc;

load('data.mat')

DTU10MWrigidhawc2su8000(1,:)=[];
DTU10MWrigidhawc2su8001(1,:)=[];
DTU10MWrigidhawc2su8002(1,:)=[];
DTU10MWrigidhawc2su8003(1,:)=[];
DTU10MWrigidhawc2su8004(1,:)=[];
DTU10MWrigidhawc2su8005(1,:)=[];
DTU10MWrigidhawc2su8006(1,:)=[];

figure()
colormap winter
subplot(2,3,1)
    plot(DTU10MWrigidhawc2su8000.Sm1, DTU10MWrigidhawc2su8000.A2,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8001.Sm1, DTU10MWrigidhawc2su8001.A2,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8002.Sm1, DTU10MWrigidhawc2su8002.A2,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8003.Sm1, DTU10MWrigidhawc2su8003.A2,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8004.Sm1, DTU10MWrigidhawc2su8004.A2,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8005.Sm1, DTU10MWrigidhawc2su8005.A2,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8006.Sm1, DTU10MWrigidhawc2su8006.A2,'linewidth',1);hold on;
    set(gca, 'FontName', 'Times New Roman','FontSize',11)
    xlabel('R [m] ','Interpreter','LaTex')
    ylabel(' a [-]','Interpreter','LaTex'); 
    grid on; box on;
    title('Axial Induction factor','FontSize',13)

subplot(2,3,4)
    plot(DTU10MWrigidhawc2su8000.Sm1, DTU10MWrigidhawc2su8000.AP3,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8001.Sm1, DTU10MWrigidhawc2su8001.AP3,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8002.Sm1, DTU10MWrigidhawc2su8002.AP3,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8003.Sm1, DTU10MWrigidhawc2su8003.AP3,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8004.Sm1, DTU10MWrigidhawc2su8004.AP3,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8005.Sm1, DTU10MWrigidhawc2su8005.AP3,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8006.Sm1, DTU10MWrigidhawc2su8006.AP3,'linewidth',1);hold on;    set(gca, 'FontName', 'Times New Roman','FontSize',11)
    xlabel('R [m] ','Interpreter','LaTex')
    ylabel(' a_p [-]');
    grid on; box on;
    title('Tangential Induction factor','FontSize',13)
subplot(2,3,2)
    plot(DTU10MWrigidhawc2su8000.Sm1, DTU10MWrigidhawc2su8000.CL017,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8001.Sm1, DTU10MWrigidhawc2su8001.CL017,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8002.Sm1, DTU10MWrigidhawc2su8002.CL017,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8003.Sm1, DTU10MWrigidhawc2su8003.CL017,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8004.Sm1, DTU10MWrigidhawc2su8004.CL017,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8005.Sm1, DTU10MWrigidhawc2su8005.CL017,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8006.Sm1, DTU10MWrigidhawc2su8006.CL017,'linewidth',1);hold on;    
    set(gca, 'FontName', 'Times New Roman','FontSize',11)
    xlabel('R [m] ','Interpreter','LaTex')
    ylabel(' C_L [-]'); 
    grid on; box on;
    title('Lift Coefficient','FontSize',13)
subplot(2,3,5)
    plot(DTU10MWrigidhawc2su8000.Sm1, DTU10MWrigidhawc2su8000.CD018,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8001.Sm1, DTU10MWrigidhawc2su8001.CD018,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8002.Sm1, DTU10MWrigidhawc2su8002.CD018,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8003.Sm1, DTU10MWrigidhawc2su8003.CD018,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8004.Sm1, DTU10MWrigidhawc2su8004.CD018,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8005.Sm1, DTU10MWrigidhawc2su8005.CD018,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8006.Sm1, DTU10MWrigidhawc2su8006.CD018,'linewidth',1);hold on;    
    set(gca, 'FontName', 'Times New Roman','FontSize',11)
    xlabel('R [m] ','Interpreter','LaTex')
    ylabel(' C_D [-]'); 
    grid on; box on;
    title('Drag Coefficient','FontSize',13)
subplot(2,3,3)
    plot(DTU10MWrigidhawc2su8000.Sm1, DTU10MWrigidhawc2su8000.CT33,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8001.Sm1, DTU10MWrigidhawc2su8001.CT33,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8002.Sm1, DTU10MWrigidhawc2su8002.CT33,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8003.Sm1, DTU10MWrigidhawc2su8003.CT33,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8004.Sm1, DTU10MWrigidhawc2su8004.CT33,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8005.Sm1, DTU10MWrigidhawc2su8005.CT33,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8006.Sm1, DTU10MWrigidhawc2su8006.CT33,'linewidth',1);hold on;        
    set(gca, 'FontName', 'Times New Roman','FontSize',11)
    xlabel('R [m] ','Interpreter','LaTex')
    ylabel(' C_T [-]');
    grid on; box on;
    title('Thrust Coefficient','FontSize',13)
subplot(2,3,6)
    plot(DTU10MWrigidhawc2su8000.Sm1, DTU10MWrigidhawc2su8000.CP34,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8001.Sm1, DTU10MWrigidhawc2su8001.CP34,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8002.Sm1, DTU10MWrigidhawc2su8002.CP34,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8003.Sm1, DTU10MWrigidhawc2su8003.CP34,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8004.Sm1, DTU10MWrigidhawc2su8004.CP34,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8005.Sm1, DTU10MWrigidhawc2su8005.CP34,'linewidth',1);hold on;
    plot(DTU10MWrigidhawc2su8006.Sm1, DTU10MWrigidhawc2su8006.CP34,'linewidth',1);hold on;        
    set(gca, 'FontName', 'Times New Roman','FontSize',11)
    xlabel('R [m] ','Interpreter','LaTex')
    ylabel(' C_p [-]');
    
    grid on; box on;
    title('Power Coefficient','FontSize',13)
    legend('TSR 6','TSR 6.5','TSR 7','TSR 7.5','TSR 8','TSR 8.5','TSR 9','Interpreter','LaTex')
    
  
set(gcf,'PaperUnits','centimeters','PaperSize',[20 13])
fig = gcf;fig.PaperUnits = 'centimeters'; 
fig.PaperPosition = [0 0 20 13];fig.Units = 'centimeters';
fig.PaperSize=[20 13];fig.Units = 'centimeters';
print(fig,'VANPLOT','-dpdf','-r200')





