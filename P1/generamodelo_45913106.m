function model=generamodelo_45913106(numberofnodes,tam,mode)
    %% Funcion que genera un modelo plano de tri�ngulos
    %  Contact info: Rub�n Molero Alabau (rumoal1@itaca.upv.es)
    %Andreu M. Climent (andreu.climent@gmail.com)
    
    % Este modelo tendr� N nodos y M triangulos en una superficie cuadrada
    % de lado T.
    % model=generamodelo_profesor(numberofnodes,tam)
    %% INPUTS:
    % numberofnodes: number of nodes of the model
    % tam: size of the map (tam=10 genera el modelo sobre superficie 10x10)
    % mode: indica el modo de generaci�n los nodos del modelo
    %       modo=0; % los nodos se genera aleatoriamente
    %       modo=1; % los nodos los indica el usuario manualmente
    %       modo=2; % utilizar los nodos almacenados
    %% OUTPUTS
    % model    estructura con la informaci�n del modelo
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
%   Obtengo los �ndices de fila y columna donde se encuentra el nodo i
    [fila,columna]=find(dt.ConnectivityList==i);
%   Crea el vector con todas las fila(puntos que forman cada triangulo) en
%   el nodo i
    nei=[nei dt.ConnectivityList(fila,:)];
%   Quita nodos repetidos en el vector nei
    aux=unique(nei{i}); 
    index=find(aux==i); %Devuelve el �ndice que apunta al propio nodo
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

