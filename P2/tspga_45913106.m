% function 1.	tspga_DNI (fmode,model)
function tspga_45913106(fmode,model)
% Función que utiliza algoritmos genéticos para resolver el problema del
% viajante (Travel Sales Problem TSP)
% INPUTS
    % mode: indica el modo de generación los nodos del modelo
    %       está función llamará a la función generamodelo_DNI creada en la
    %       primera práctica
    %    modo=0; % los nodos se genera aleatoriamente
    %    modo=1; % los nodos los indica el usuario manualmente
    %    modo=2; % utilizar los nodos almacenados
    %    mode=3; % la estructura del modelo se introduce como una variable
    % model: estructura generada previamente con generamodelo_DNI.m
% OUPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERACION DEL MODELO EN FUNCIÓN DE LAS INSTRUCCIONES DEL USUARIO %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
     fmode=0;
end
numberofnodes=20;
tam=10;
if fmode==3
    if nargin<2
         disp('Error: El modo 3 requiere introducir el modelo')
    end
end
if fmode<3
    model=generamodelo_45913106(numberofnodes,tam,fmode);
end
x=model.x;
y=model.y;
numberofnodes=length(x);
save('model','model');
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETROS DEL ALGORITMO GENÉTICO %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npar=numberofnodes; % # of optimization variables
Nt=npar; % # of columns in population matrix
maxit=5000; % max number of iterations
popsize=20; % set population size / miembros de la población
mutrate=.05; % set mutation rate
selection=0.5; % fraction of population kept / fracción de miembros sobreviven
iga=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INICIALIZA LA POBLACIÓN Y VECTOR DE COSTES%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pop= ones(popsize,npar);
for i=1:popsize
    pop(i,:)=randperm(npar);
end

cost=tspfun(pop,model); %función que evalúa el fitness de cada miembro de la poblacion
[cost,ind]=sort(cost); %devuelve el vector ordenado y los indices del antiguo vector.
pop=pop(ind,:); %reordena la población de acuerdo al coste



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EN CADA ITERACIÓN VAMOS A IR GUARDANDO EL VALOR DE COSTE DEL MEJOR INDIVIDUO Y DE LA POBLACIÓN MEDIA PARA POSTERIORMENTE PODER VER EL PROGRESO%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 minc(1)=min(cost); % minc contains min of population
 meanc(1)=mean(cost); % meanc contains mean of population
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERO EL VECTOR DE PROBABILIDADES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g. si tengo una población de 10 individuos y solo la mitad sobreviven
% para reproducirse entonces el vector de probabilidades sería.
% Prob=[5 4 4 3 3 3 2 2 2 2 1 1 1 1 1];
% Seleccionamos a la mejor mitad de la población.

        odds=[];
        k=popsize*selection;
    for i=1:popsize*selection
        for j=1:k
            odds=[i odds];        
        end
            k=k-1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PONEMOS EL MUNDO VIRTUAL A FUNCIONAR (MAIN LOOP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while iga<maxit
    iga=iga+1; % increments generation counter
    
    pop=pop([1:popsize*selection],:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ELIJO A LOS PADRES Y MADRES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % De entre el vector de probabilidades elijo aleatoriamente
    % a los padres y a las madres.
        N=round(popsize*selection/2);
        indicespadres=ceil(length(odds)*rand(1,N));
        indicesmadres=ceil(length(odds)*rand(1,N));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Emparejamientos y generación de los hijos %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hijos=[];
    for i=1:N %Bucle de emparejamiento e hijos
        
        %Obtenemos quien son los padres
        padre=odds(indicespadres(i)); 
        madre=odds(indicesmadres(i));
        %Cogemos los genes de la matriz pop
        genespadre = pop(padre,:);
        genesmadre = pop(madre,:);      

        at=randi(npar);
        %Obtenemos los dos hijos de los dos padres seleccionamos
        hijo1=op_cruce(genespadre,genesmadre,npar,at);
        hijo2=op_cruce(genesmadre,genespadre,npar,at);   
        
        
        
        hijos=[hijos;hijo1;hijo2];
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate the population
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        col1=0;
        col2=0;
        while(col1==col2) %Evitar que los genes que mutan coincidan
            col1=randi(Nt);
            col2=randi(Nt);
        end
        hijos(:,[col1 col2])=hijos(:,[col2 col1]); %Mutacion
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    pop=[pop;hijos];
    cost=tspfun(pop,model);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

%     Sort the costs and associated parameters
    [cost,ind]=sort(cost);
    pop=pop(ind,:);
    
%     Do statistics
    minc(iga)=min(cost);
    meanc(iga)=mean(cost);
        
    iga
end %iga
 


% Displays the output
day=clock;
disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5),day(6)),0))
%disp(['optimized function is ' ff])
format short g
disp(['popsize=' num2str(popsize) 'mutrate=' num2str(mutrate) '# par=' num2str(npar)])
disp([' best cost=' num2str(cost(1))])
disp(['best solution']); disp([num2str(pop(1,:))])
 
figure(2)
iters=1:maxit;
plot(iters,minc,iters,meanc,'--','LineWidth',2);
hold on
xlabel('generation','FontWeight','bold','FontSize',12);
ylabel('cost','FontWeight','bold','FontSize',12);
legend('Coste Mínimo', 'Coste Medio')
box on
set(gca,'FontSize',12)
set(gca,'FontWeight','bold')
plotsolution(x,y,pop)


end

function cost=tspfun(pop,model)
% Cost function for the traveling salesman problem.
[popsize,npar]=size(pop);
cost=[];
for i=1 :popsize
coste=[];
genes=pop(i,:);
    for l=2 :npar
        nodoi=genes(l-1);%Obtenemos el nodo i
        nodoj=genes(l);%Obtenemos el nodo j
        coste=[coste model.D(nodoi,nodoj)]; %Añadimos la distancia entre cada salto del camino seguido por la hormiga
                
    end
    coste=sum(coste);
    cost=[cost coste];
end
end
function plotsolution(x,y,pop)
 
%%%%%%%%%%%%%%%%%%%
% PINTA LOS NODOS %
%%%%%%%%%%%%%%%%%%%
figure(1);
clf %clear figure
tam=max([x; y])+1; % tamaño máximo del modelo
plot(x,y,'x') %pinta los nodos
axis([0 tam 0 tam]) %fija los ejes
 
%%%%%%%%%%%%%%%%%%%%%
% PINTA LA SOLUCIÓN %
%%%%%%%%%%%%%%%%%%%%%
hold on
sol=pop(1,:);
for i=1:length(sol)-1
    p1x=x(sol(i));
    p2x=x(sol(i+1));
    p1y=y(sol(i));
    p2y=y(sol(i+1));
    line([p1x p2x], [p1y p2y],'LineWidth',2)         
end
    p1x=x(sol(1));
    p2x=x(sol(end));
    p1y=y(sol(1));
    p2y=y(sol(end));
    line([p1x p2x], [p1y p2y],'LineWidth',2)      
axis square
 
hold on
%Plot the number of each node
for i=1:length(x)    
    text(x(i),y(i),num2str(i),'FontWeight','bold','FontSize',12)
end
box on
xlabel('X coordenate (m)','FontWeight','bold')
ylabel('Y coordenate (m)','FontWeight','bold')
set(gca,'FontSize',12)
set(gca,'FontWeight','bold')
end
function [hijo] = op_cruce(genespadre,genesmadre,npar,at)
%Metodo para realizar la operacion de cruce
%npar = numero de genes
i=at-1;
hijo=zeros(1,npar);
auxmadre=genesmadre;
counter1=npar/2;
        while(counter1~=0)
            i = i +1;
            if i > npar
                i=1;
            end  
            hijo(i)=genespadre(i);
            counter1=counter1 - 1;
            indi=auxmadre==hijo(i); %Indice del gen utilizado por el otro padre
            auxmadre(:,indi)=[]; %Contiene los genes restantes que no han sido utilizado por el otro padre
            
        end
        
        j=1;
        counter2=npar/2;
     
        %Genes de la madre
        while(counter2~=0)
            i=i+1;
            if i > npar
                i=1;                              
            end
            hijo(i)=auxmadre(j); 
            counter2=counter2 - 1;
            
            j=j+1;
        end
end
function model=generamodelo_45913106(numberofnodes,tam,mode)
    %% Funcion que genera un modelo plano de triángulos
    %  Contact info: Rubén Molero Alabau (rumoal1@itaca.upv.es)
    %Andreu M. Climent (andreu.climent@gmail.com)
    
    % Este modelo tendrá N nodos y M triangulos en una superficie cuadrada
    % de lado T.
    % model=generamodelo_profesor(numberofnodes,tam)
    %% INPUTS:
    % numberofnodes: number of nodes of the model
    % tam: size of the map (tam=10 genera el modelo sobre superficie 10x10)
    % mode: indica el modo de generación los nodos del modelo
    %       modo=0; % los nodos se genera aleatoriamente
    %       modo=1; % los nodos los indica el usuario manualmente
    %       modo=2; % utilizar los nodos almacenados
    %% OUTPUTS
    % model    estructura con la información del modelo
    % model.dt estructura del modelo en triangulos
    % model.dt.X: vector de dimension Nx2 indicando las coordinadas X e Y de 
    %           cada nodo del modelo
    % model.dt.Triangulation: vector de dimension [Mx3] indicando el numero
    %           de los tres triangulos vertices de cada triangulo
    % model.x  coordenadas x de los puntos del modelo
    % model.y  coordenadas y de los puntos del modelo
    % model.D  matriz de NxN dimensiones con la distancia entre cada nodo
    % model.nei estructura donde se almacenan los vecinos de cada nodo
    %           e.g. model.nei{1}= [2 3 5] indica que el nodo 1 es vecino
    %           de los nodos 2, 3 y 5.

close all % Close all opened figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF INPUT VARIABLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    mode=0;
end
if nargin<2
    tam=100;
end
if nargin<1 
    numberofnodes=10; 
end

%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF NODES %
%%%%%%%%%%%%%%%%%%%%%%
%% MODO 0: Aleatoria
if mode==0
x=tam+(-tam)*rand(numberofnodes,1);
y=tam+(-tam)*rand(numberofnodes,1);
end

%% MODO1: Manualmente
if mode==1
    [x,y]=ginput(numberofnodes);
    x=tam+(-tam)*x;
    y=tam+(-tam)*y;
    close
end
%% MODO2: Almacenado
if mode==2
    auxtam=tam;%almacenamos variable introducida
    load('var','x','y','tam');
    x=x/tam;
    y=y/tam;
    tam=auxtam;%devolvemos el valor introducido por el usuario
    x=tam*x;
    y=tam*y;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF TRIANGLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('var.mat','x','y','tam');
dt=delaunayTriangulation(x,y);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION OF MATRIX OF DISTANCES BETWEEN NODES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = ones(numberofnodes,numberofnodes);
for i=1:numberofnodes
    Xi=dt.Points(i,1);
    Yi=dt.Points(i,2);
    
    for j=1:numberofnodes
        Xj=dt.Points(j,1);
        Yj=dt.Points(j,2);
        dist=sqrt(((Xi-Xj)^2)+((Yi-Yj)^2));
        D(i,j)=dist*D(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION OF NEIGHBORS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
nei={};
for i=1:numberofnodes
%   Obtengo los índices de fila y columna donde se encuentra el nodo i
    [fila,columna]=find(dt.ConnectivityList==i);
%   Crea el vector con todas las fila(puntos que forman cada triangulo) en
%   el nodo i
    nei=[nei dt.ConnectivityList(fila,:)];
%   Quita nodos repetidos en el vector nei
    aux=unique(nei{i}); 
    index=find(aux==i); %Devuelve el índice que apunta al propio nodo
    aux=aux(:); %Convierte en vector columna para trabajar siempre igual
    aux(index,:)=[]; %Elimina el propio nodo del vector aux de vecinos
    nei{i}=aux; %Almacena los vecinos
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE OUTPUT VARIABLE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.dt=dt; 
model.x=x;
model.y=y;
model.D=D; %distance between each node
model.nei=nei; %estructure with the neighbors of each node
 
%%%%%%%%%%%%%%%
% PLOT OUTPUT %
%%%%%%%%%%%%%%%
figure
set(gca,'FontSize',12) %# Fix font size of the text in the current axes 
set(gca,'FontWeight','bold')  %# Fix Bold text in the current axes 
plot(model.x,model.y,'x') 
triplot(model.dt,'r');  %# Plot the Delaunay triangulation
for i=1:length(model.dt.Points) %# Plot the number of each node
    text(model.dt.Points(i,1),model.dt.Points(i,2),num2str(i),'FontWeight','bold')
end
axis([0 tam 0 tam]) %# Fix axes representation size
box on %# Plot all lines in the box of the figure
xlabel('X coordenate (m)') %# Label in the x axes
ylabel('Y coordenate (m)') %# Label in the y axes
title('Generated model') %# Title
end
