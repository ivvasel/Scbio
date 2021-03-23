function electro = calculaelectromiogramas_45913106(r,activaciones)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
load('AP_ventriculo.mat');
% r=30; %retrasito
% activaciones=5;
AP_sano_aux=[];
AP_isch_aux=[];
t=length(AP_sano)*activaciones;


for i=1:activaciones+1 %Concatenacion extra para cortar con el r necesario
    AP_sano_aux = [AP_sano_aux AP_sano];
    AP_isch_aux = [AP_isch_aux AP_isch];
end
AP_sano_r = AP_sano_aux(r:activaciones*length(AP_sano)+r-1);
AP_sano = AP_sano_aux(1:activaciones*length(AP_sano));

AP_isch_r = AP_isch_aux(r:activaciones*length(AP_isch)+r-1);
AP_isch = AP_isch_aux(1:activaciones*length(AP_isch));

EGM_sano= AP_sano-AP_sano_r;
EGM_ischemico= AP_isch-AP_isch_r;
EGM_transicion= AP_sano-AP_isch_r;

%%OUTPUTS%%
electro.AP_sano = AP_sano;
electro.AP_sano_r = AP_sano_r;

electro.AP_isch = AP_isch;
electro.AP_isch_r = AP_isch_r;

electro.EGM_sano = EGM_sano;
electro.EGM_isch = EGM_ischemico;
electro.EGM_trans = EGM_transicion;

electro.t = t;


% close all
% figure
% subplot(2,1,1)
% plot(1:t,AP_sano,'LineWidth',2),hold on
% plot(1:t,AP_sano_r,'LineWidth',2),hold on
% legend('Celula 1','Celula 2');
% xlabel('Time (ms)','FontWeight','bold') 
% ylabel('Voltage (mV)','FontWeight','bold') 
% title('Potenciales de acción','FontWeight','bold') 
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold') 


% subplot(2,1,2)
% plot(EGM_sano,'LineWidth',2)
% xlabel('Time (ms)','FontWeight','bold') 
% ylabel('Voltage (mV)','FontWeight','bold') 
% title('Electromiograma','FontWeight','bold') 
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold') 

% figure
% subplot(2,1,1)
% plot(AP_sano,'LineWidth',2),hold on
% plot(AP_sano_r,'LineWidth',2),hold on
% legend('Celula 1','Celula 2');
% xlabel('Time (ms)','FontWeight','bold') 
% ylabel('Voltage (mV)','FontWeight','bold') 
% title('Potenciales de acción','FontWeight','bold') 
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold') 
% 
% subplot(2,1,2)
% plot(EGM_sano,'LineWidth',2)
% xlabel('Time (ms)','FontWeight','bold') 
% ylabel('Voltage (mV)','FontWeight','bold') 
% title('Electromiograma','FontWeight','bold') 
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')

end