function mainHH=mainHH_Euler_45913106(gbarNa)
%==============================
%% HODKING-HUXLEY MODEL NEURAL CELL MODEL
%============================== 

if nargin < 1
    gbarNa = 1.2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constants set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cm=0.01; % Membrane Capcitance mF/cm^2
dt=0.05; % Time Step ms
t=0:dt:100; %Time Array ms
Iap=0.1; %External Current Applied
ENa=55.17; % mv Na reversal potential
EK=-72.14; % mv K reversal potential
El=-49.42; % mv Leakage reversal potential
% gbarNa=1.2; % mS/cm^2 Na conductance
gbarK=0.36; % mS/cm^2 K conductance
gbarl=0.003; % mS/cm^2 Leakage conductance
%Potenciales de Nerst
VNa=55.17;
Vk=-72.14;
Vl=-49.42;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=zeros(1,length(t));
    n=zeros(1,length(t));
    h=zeros(1,length(t));
    gNa=zeros(1,length(t));
    gK=zeros(1,length(t));
    gl=zeros(1,length(t));
    INa=zeros(1,length(t));
    IK=zeros(1,length(t));
    Il=zeros(1,length(t));
    V=zeros(1,length(t));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1)=-60; % Initial Membrane voltage
v=V(1);
am=0.1*(v+35)/(1-exp(-(v+35)/10));
bm=4.0*exp(-0.0556*(v+60));
m(1)=am/(am+bm); % Initial m-value
an=0.01*(v+50)/(1-exp(-(v+50)/10));
bn=0.125*exp(-(v+60)/80);
n(1)=an/(an+bn); % Initial n-value
ah=0.07*exp(-0.05*(v+60));
bh=1/(1+exp(-(0.1)*(v+30)));
h(1)=ah/(ah+bh); % Initial h-value
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP CALCULATING VM FOR EACH INSTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for i=1:length(t)-1
  
    
%En este bucle se han de calcular
% 1) A partir del valor de las compuertas (i.e. m,n y h) en el intante 
%    i calculamos el valor de las conductancias en dicho instante (i.e.
%     (gNa(i) , gK(i) , gl(i))
    gNa(i) = gbarNa*m(i)^3*h(i);
    gK(i) = gbarK*n(i)^4;
    gl(i) = gbarl;

    
% 2) Una vez conocida la conductancia se calcula el valor de cada
%    corriente en el instante i (INa(i), IK(i) e Il(i))

    INa(i) = -gNa(i)*(V(i)-VNa);
    IK(i) = -gK(i)*(V(i)-Vk);
    Il(i) = -gl(i)*(V(i)-Vl);

    
% 3) Mediante el método de Euler y sabiendo el valor de cada corriente 
%    en el instante i, así como el valor de la tensión en el instante i
%    podemos predecir cual será el valor de la tensión en el instante
 V(i+1)=V(i)+(dt)*((1/Cm)*(Iap+INa(i)+IK(i)+Il(i))); 

    
 % 4) También mediante el método de Euler podemos predecir cual será el
 %    valor de las compuertas en el instante i+1 (i.e. m(i+1), n(i+1) y
 %    h(i+1)     
 m(i+1)=m(i)+dt*(am*(1-m(i))-bm*m(i)); 
 n(i+1)=n(i)+dt*(an*(1-n(i))-bn*n(i)); 
 h(i+1)=h(i)+dt*(ah*(1-h(i))-bh*h(i)); 

    
% Ya disponemos de todo lo necesario para pasar al siguiente instante
% del bucle

v=V(i);
am=0.1*(v+35)/(1-exp(-(v+35)/10));
bm=4.0*exp(-0.0556*(v+60));
m(i)=am/(am+bm); % Initial m-value
an=0.01*(v+50)/(1-exp(-(v+50)/10));
bn=0.125*exp(-(v+60)/80);
n(i)=an/(an+bn); % Initial n-value
ah=0.07*exp(-0.05*(v+60));
bh=1/(1+exp(-(0.1)*(v+30)));
h(i)=ah/(ah+bh);

end

mainHH.V = V;
mainHH.dt = dt;
mainHH.t = t;
mainHH.INa = INa;
mainHH.IK = IK;
mainHH.Il = Il;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ACTION POTENTIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
 
% plot(t,V,'LineWidth',2),hold on
% legend('Action Potential');
% xlabel('Time (ms)','FontWeight','bold') 
% ylabel('Voltage (mV)','FontWeight','bold') 
% title('Voltage Change for Hodgkin-Huxley Model','FontWeight','bold') 
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')  
%  
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT ION CHANNEL CURRENTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(t(1:length(INa)),INa,'g',t(1:length(INa)),IK,'r',t(1:length(INa)),Il,'k','LineWidth',2),hold on
% ylabel('Ion channel currents','FontWeight','bold') 
% xlabel('Time (ms)','FontWeight','bold') 
% %axis([0 5 0 1])
% legend('INa','IK','Ir');
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')  
end



