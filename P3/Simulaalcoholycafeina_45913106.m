function [cont,alc,caf]=Simulaalcoholycafeina_45913106
gbarNa = 1.2;
gbarNaAlc = 0.8*gbarNa;
gbarNaCaf = 1.2*gbarNa;

cont = mainHH_Euler_45913106(gbarNa);
alc = mainHH_Euler_45913106(gbarNaAlc);
caf = mainHH_Euler_45913106(gbarNaCaf);

% figure
% plot(alc.t,alc.V,'LineWidth',2),hold on
% plot(cont.t,cont.V,'LineWidth',2),hold on
% plot(caf.t,caf.V,'LineWidth',2),hold on
% legend('Alcohol','Control','Cafeina');
% xlabel('Time (ms)','FontWeight','bold') 
% ylabel('Voltage (mV)','FontWeight','bold') 
% title('Voltage Change for Hodgkin-Huxley Model','FontWeight','bold') 
% % set(gca,'FontSize',12)
% % set(gca,'FontWeight','bold') 

end