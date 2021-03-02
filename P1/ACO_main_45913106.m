%% PRACTICA 1 SISTEMAS COMPLEJOS BIOINSPIRADOS
% Este archivo contine 3 funciones principales:
% 1. generamodelo_45913106
% 2. aco_enrutamiento_45913106
% 3. aco_TSP_45913106
clear all
clc

%% Interfaz de usuario
disp('Elige una función:')
disp('1: Genera modelo')
disp('2: Aco_enrutamiento')
disp('3: Aco_TSP')

funcion=input('Elección: ');

switch funcion
    case 1
        disp('Modo de generación del modelo')
        disp('0: Aleatoria')
        disp('1: Manual')
        disp('2: Almacenado')
        modo=input('Modo: ');
        switch modo
            case 0
                clc
                disp('Introduce: generamodelo_45913106(numberofnodes,tam,0)')
                disp('numberofnodes: número de nodos')
                disp('tam: tamaño máximo de la cuadrícula')
            case 1
                clc
                disp('Introduce: generamodelo_45913106(numberofnodes,tam,1)')
                disp('numberofnodes: número de nodos')
                disp('tam: tamaño máximo de la cuadrícula')
                disp('Deberás clickar en la figura la cantidad de numberofnodes veces') 
            case 2
                clc
                disp('Introduce: generamodelo_45913106(numberofnodes,tam,2)')
                disp('numberofnodes: número de nodos')
                disp('tam: tamaño máximo de la cuadrícula') 
                disp('Si deseas utilizar nodos almacenados propios debes de crear el archivo con save(`var.mat´,´x ´,`y`,`tam`)')
                disp('Si no deberás crear primero un modelo con alguno de los modos anteriores y después ejectutar este modo')
        end
        
    case 2
        disp('Elige el modo de funcionamiento')
        disp('0: Genera un modelo aleatorio')
        disp('1: Utilizar un modelo ya generado')
        modo=input('Modo: ');
        switch modo
            case 0
                clc
                disp('Introduce: aco_enrutamiento_45913106()')
                disp('Generará un nuevo modelo con los parámetros predefinidos')
                disp('Este modo guarda en un arvhico modelo.mat el modelo que ha generado')
            case 1
                clc
                disp('Introduce: aco_enrutamiento_45913106(model,nest,food,nants,niter)')
                disp('model: modelo almacenado')
                disp('nest: nodo del hormiguero')
                disp('food: nodo de la comida')
                disp('nants: número de hormigas')
                disp('niter: número de veces que una hormiga sale a buscar comida')
        end
        
    case 3
        disp('Elige el modo de funcionamiento')
        disp('0: Genera un modelo aleatorio')
        disp('1: Utilizar un modelo ya generado')
        modo=input('Modo: ');
        switch modo
            case 0
                clc
                disp('Introduce: aco_TSP_45913106()')
                disp('Generará un nuevo modelo con los parámetros predefinidos')
                disp('Este modo guarda en un arvhico modelo.mat el modelo que ha generado')
            case 1
                clc
                disp('Introduce: aco_TSP_45913106(model,nants,niter)')
                disp('model: modelo almacenado')
                disp('nants: número de hormigas')
                disp('niter: número de veces que una hormiga sale a buscar el camino')
        end
        
end

%% Funciones utilizadas
%% GENERADORA DE MODELOS
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

%% ACO_ENRUTAMIENTO
function [BestSol]=aco_enrutamiento_45913106(model,nest,food,nants,niter)
%% Funcion que resuelve el problema de enrutamiento mediante optimización 
% por colonias de hormigas
    % Contact info: Rubén Molero Alabau (rumoal1@itaca.upv.es)
% Andreu M. Climent (andreu.climent@gmail.com)
    % [BestSol]=aco_enrutamiento_profesor(model,nest,food,nants,niter)
    %% INPUTS:
    % model    estructura con la información del modelo generada con
    %           generamodelo_XXX.m 
    %nest:     Número del nodo que se considerará origen
    %food:     Número del nodo que se considerará destino
    %nants:    Número de hormigas
    %niter:    Número de iteraciones máximas
    % OUTPUT:
    % BestSol:  Estructura con la mejor ant
    %   BestSol.Tabu :vector con la secuencia del camino optimo detectado
    %   BestSol.Cost :distancia del camino elegido.   %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF EMPTY INPUT VARIABLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<1, model=generamodelo_45913106(30,10,0); end %generate a random model
    if nargin<5, niter=100; end
    if nargin<4, nants=40; end    
    if nargin<3, 
        [val,pos]=max(model.x);
        food=pos; 
    end %food will be the last node
    if nargin<2, 
        [val,pos]=min(model.x);
        nest=pos; 
    end %nest will be the first node 
    
    save('modelo','model');
%%%%%%%%%%%%%%%%%
% ACO PARAMETERS %
%%%%%%%%%%%%%%%%%
MaxIt=niter;      % Maximum Number of Iterations
nAnt=nants;        % Number of Ants (Population Size)
alpha=1;        % Phromone Exponential Weight
beta=1;         % Heuristic Exponential Weight
rho=0.05;       % Evaporation Rate
Q=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE eta (visibility) and tau (PHEROMONE) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=1./model.D; %Matriz de visibilidad
nVar=length(model.x); %Número de nodos
tau0=10*Q/(nVar*mean(model.D(:))); %Feromona inicial
tau=tau0*ones(nVar,nVar); %Matriz de todas las feromonas


%%%%%%%%%%%%%%%%%
% ACO Main Loop %
%%%%%%%%%%%%%%%%%
for k=1:MaxIt
% Dentro del este bucle principal probablemente necesitarás las siguientes secciones:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize each Ants at nest %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nants
   ant(n).Tabu=nest;
   ant(n).Cost=[];
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct a solution for each ANT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exitosas=0;
    exito = false;
while(~exito)
for n=1:nants %Recorre todas las hormigas    
    antPath=ant(n).Tabu;
    pos=antPath(length(antPath)); %Posición actual obtenida del vector del camino   
    
    %%Comprueba si la posicion actual es en la comida y salta a la
    %%siguiente hormiga
    if pos==food
        continue     
    
    %%No ha llegado a la comida y sigue buscando
    else
        auxNei=model.nei{pos};%Los vecinos en la posicion de la hormiga
        checkNei=auxNei;%Almacena los vecinos en una variable variante
        allP=[];%Inicializa el vector de todas las prob entre la hormiga y los vecinos
        
        %%calcula tau y eta de la posicion con los posibles vecinos que no se hayan visitado ya. 
        %%Obtenemos el vector con las probabilidades de saltar a cada
        %%posible vecino
        
        for p=1:length(auxNei) %Indices de la matrix auxiliar de vecinos        
            posiblesalto=auxNei(p); %Obtiene un posible salto de model.nei{p}=auxNei(p)
            comprobar=find(antPath==posiblesalto); %Comprobar si posiblesalto ya se encuentra en el camino de la ant
            if isempty(comprobar) %Comprueba si se ha devuelto algún indice indicando que por ese nodo ya ha pasado la hormiga
                %Si no hay indice indica que la hormiga no ha pasado por
                %ese nodo
                %Calculamos probabilidad hacia el posiblesalto
                TauVecino=tau(pos,posiblesalto);%Obtenemos la feromona del camino hacia el nodo vecino
                EtaVecino=eta(pos,posiblesalto);%Obtenemos la visibilidad del nodo vecino         
                allP=[allP (TauVecino^alpha)*(EtaVecino^beta)]; %Vecotr que almacena los valores del sumatorio del denominador de la formula de pij       
        
            else %Si la hormiga ha pasado ya por el posible salto entonces:
                elimNei=find(checkNei==posiblesalto); %Obtiene el indice donde se encuentra el nodo que se desea eliminar de posibles saltos
                checkNei(elimNei,:)=[]; %Elimina el nodo del vecto. Se utiliza una variable extra para no interferir en el for que utiliza auxNei.            
            end 
        
        end
    
        sumaP = sum(allP); %Sumatorio del denominador de la formula pij
        P=[]; %Inicializa el vector P donde se almacenarán las probabilidades del cumsum
    
    
        %%Calculamos la probabilidad para cada caso y creamos el vector P
        for p=1:length(allP) %Recorre todo el vector sumaP que indica cuantas probabilidades se han calculado y por tanto cuantos posibles vecinos hay
        P=[P allP(p)/sumaP]; %Almacena las probabilidades pij
        end
    
    
        %%Elegimos de forma heuristica el nodo al que iremos
        r=rand;
        C=cumsum(P);
        j=find(r<=C,1,'first'); %selecciona el primero mayor que r y devuelve el indice
        %%Elige un posible salto entre todos y se añade al vector de nodos
        %%recorridos por la hormiga
        ant(n).Tabu=[antPath checkNei(j)];%introduce en el path el nodo seleccionado
    
        %%Comprueba si no hay vecinos a los que ir y reinicia el path
        %%Exactamente comprueba si se han eliminado todos los posibles
        %%saltos en el proceso anterior, esto indicará que el nodo no ha
        %%encontrado un camino por donde avanzar.
        if isempty(checkNei)
            antPath=[nest]; %Devuelve a la hormiga al nest
            ant(n).Tabu=antPath; %Reinicia la variable 
            continue %Salta a la siguiente hormiga
        end
        
        %%Comprueba si ha llegado la hormiga a la comida
        antPath=ant(n).Tabu; %Actualiza la varible con el avance de la hormiga incluido
        pos=antPath(length(antPath));%Actualiza la variable con el nodo actual de la hormiga
        if pos==food %Compara la el nodo de la hormiga con el de la comida
            exitosas=exitosas+1; %Si encuentra la comida suma +1 a las hormiga que han logrado llegar a la comida
            if exitosas==nants %Comprueba si todas las hormigas han llegado
                exito=true; %En de llegar todas se cambia el valor de la booleana exito
                
            end
        end 

    end    
end
end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UPDATE PHEROMONES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sumCoste=[];%Inicializamos el sumatorio de los costes
        for p=1:nants %Recorre todas las hormigas
            antPath=ant(p).Tabu;%Obtenemos el camino que ha llevado la hormiga
            coste=[ant(p).Cost];%Inicializamos el coste con el valor almacenado. Aquí se almacenarán los costes de cada salto para una única hormiga(inicializado en vacio)
            
            %For para almacenar la distancia entre todos los saltos
            for l=2:length(antPath)
                nodoi=antPath(l-1);%Obtenemos el nodo i
                nodoj=antPath(l);%Obtenemos el nodo j
                coste=[coste model.D(nodoi,nodoj)]; %Añadimos la distancia entre cada salto del camino seguido por la hormiga
                
            end 
            
            coste=sum(coste); %Sumamos las distancias del vector de una hormiga y obtenemos la distancia total recorrida
            ant(p).Cost=coste; %Almacenamos la distancia de una hormiga
            sumCoste=[sumCoste coste ];%Almacenamos en un vector aux el coste de todas las hormigas. Lo utilizaremos para encontrar el minimo valor entre todas las hormigas

        end
        
        betterant = find(sumCoste==min(sumCoste)); %Apunta a la hormiga con el camino mas corto
        betterpath = ant(betterant).Tabu; %Obtiene el camino de la hormiga con menor coste
        varFijatau=mrdivide(Q,ant(betterant(1)).Cost);%Calculo la variable que es constante. Apunta al índice 1 ya que puede ser que varias hormigas tengan el mismo coste, por ranto se elige la primera hormiga
        %x = mrdivide(B,A) forma alternativa de hacer div (utilizada porque
        %pensaba que la habitual fallaba :=) )
        for l=2:length(betterpath)%Recorre los indices del mejor camino
                nodoi=betterpath(l-1); %Nodosi
                nodoj=betterpath(l); %Nodosj
                tau(nodoi,nodoj)=tau(nodoi,nodoj)+varFijatau; %Actualizamos las distancias en la matriz de feromonas
                
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tau=(1-rho)*tau; %Evaporamos una parte de todos los caminos
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %GENERATE OUTPUT VARIABLE%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        BestSol.Tabu=betterpath;
        BestSol.Cost=min(sumCoste);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Show IterationInformation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['Iteration' num2str(k) ': Best Cost = ' num2str(BestSol.Cost)]);
        %figure 
        plotdistribution(model,tau) %plot map and distribution of pheromone
        title(['Generated model. Iteration: ' num2str(k)],'FontWeight','bold')
        pause(0.1) %stop during0.1 seconds to allow the user to see the progression

end

    
end

%% ACO_TSP
function [BestSol] = aco_TSP_45913106(model,nants,niter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF EMPTY INPUT VARIABLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<1,model=generamodelo_45913106(20,10,0); end %generate a random model
    if nargin<3,niter=50; end
    if nargin<2,nants=10; end      
nest=1;
save('modelo','model');
    
    

%%%%%%%%%%%%%%%%%
% ACO PARAMETERS %
%%%%%%%%%%%%%%%%%
MaxIt=niter;      % Maximum Number of Iterations
nAnt=nants;        % Number of Ants (Population Size)
alpha=1;        % Phromone Exponential Weight
beta=1;         % Heuristic Exponential Weight
rho=0.05;       % Evaporation Rate
Q=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE eta (visibility) and tau (PHEROMONE) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=1./model.D; %Matriz de visibilidad
nVar=length(model.x); %Número de nodos
tau0=10*Q/(nVar*mean(model.D(:))); %Feromona inicial
tau=tau0*ones(nVar,nVar); %Matriz de todas las feromonas


%%%%%%%%%%%%%%%%%
% ACO Main Loop %
%%%%%%%%%%%%%%%%%
for k=1:MaxIt
% Dentro del este bucle principal probablemente necesitarás las siguientes secciones:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize each Ants at nest %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nants
   ant(n).Tabu=nest;
   ant(n).Cost=[];
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct a solution for each ANT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exitosas=0;
    exito = false;
while(~exito)
for n=1:nants %Recorre todas las hormigas    
    antPath=ant(n).Tabu;
    pos=antPath(length(antPath)); %Posición actual obtenida del vector del camino   
    
    %%Comprueba si la posicion actual es en la comida y salta a la
    %%siguiente hormiga
    if length(antPath)==length(model.x)+1 && antPath(length(antPath))==nest %if es el ultimo y esta en el objetivo salta a la siguiente hormiga
        continue     
    
    %%No ha llegado a la comida y sigue buscando
    else
        auxNei=model.nei{pos};%Los vecinos en la posicion de la hormiga
        checkNei=auxNei;%Almacena los vecinos en una variable variante
        allP=[];%Inicializa el vector de todas las prob entre la hormiga y los vecinos
        
        %%calcula tau y eta de la posicion con los posibles vecinos que no se hayan visitado ya. 
        %%Obtenemos el vector con las probabilidades de saltar a cada
        %%posible vecino
        
        for p=1:length(auxNei) %Indices de la matrix auxiliar de vecinos        
            posiblesalto=auxNei(p); %Obtiene un posible salto de model.nei{p}=auxNei(p)
            comprobar=find(antPath==posiblesalto); %Comprobar si posiblesalto ya se encuentra en el camino de la ant
            if isempty(comprobar) %Comprueba si se ha devuelto algún indice indicando que por ese nodo ya ha pasado la hormiga
                %Si no hay indice indica que la hormiga no ha pasado por
                %ese nodo
                %Calculamos probabilidad hacia el posiblesalto
                TauVecino=tau(pos,posiblesalto);%Obtenemos la feromona del camino hacia el nodo vecino
                EtaVecino=eta(pos,posiblesalto);%Obtenemos la visibilidad del nodo vecino         
                allP=[allP (TauVecino^alpha)*(EtaVecino^beta)]; %Vecotr que almacena los valores del sumatorio del denominador de la formula de pij       
        
            else %Si la hormiga ha pasado ya por el posible salto entonces:
                if posiblesalto==nest && length(antPath)==length(model.x) %comprueba que el posible salto es el objetivo y que es el penultimo nodo posible
                    TauVecino=tau(pos,posiblesalto);%Obtenemos la feromona del camino hacia el nodo vecino
                    EtaVecino=eta(pos,posiblesalto);%Obtenemos la visibilidad del nodo vecino         
                    allP=[allP (TauVecino^alpha)*(EtaVecino^beta)];%Vecotr que almacena los valores del sumatorio del denominador de la formula de pij 
                else
                elimNei=find(checkNei==posiblesalto); %Obtiene el indice donde se encuentra el nodo que se desea eliminar de posibles saltos
                checkNei(elimNei,:)=[]; %Elimina el nodo del vecto. Se utiliza una variable extra para no interferir en el for que utiliza auxNei.  
                end
            end 
        
        end
    
        sumaP = sum(allP); %Sumatorio del denominador de la formula pij
        P=[]; %Inicializa el vector P donde se almacenarán las probabilidades del cumsum
    
    
        %%Calculamos la probabilidad para cada caso y creamos el vector P
        for p=1:length(allP) %Recorre todo el vector sumaP que indica cuantas probabilidades se han calculado y por tanto cuantos posibles vecinos hay
        P=[P allP(p)/sumaP]; %Almacena las probabilidades pij
        end
    
    
        %%Elegimos de forma heuristica el nodo al que iremos
        r=rand;
        C=cumsum(P);
        j=find(r<=C,1,'first'); %selecciona el primero mayor que r y devuelve el indice
        %%Elige un posible salto entre todos y se añade al vector de nodos
        %%recorridos por la hormiga
        ant(n).Tabu=[antPath checkNei(j)];%introduce en el path el nodo seleccionado
    
        %%Comprueba si no hay vecinos a los que ir y reinicia el path
        %%Exactamente comprueba si se han eliminado todos los posibles
        %%saltos en el proceso anterior, esto indicará que el nodo no ha
        %%encontrado un camino por donde avanzar.
        if isempty(checkNei)
            antPath=[nest]; %Devuelve a la hormiga al nest
            ant(n).Tabu=antPath; %Reinicia la variable 
            continue %Salta a la siguiente hormiga
        end
        
        %%Comprueba si ha llegado la hormiga a la comida
        antPath=ant(n).Tabu; %Actualiza la varible con el avance de la hormiga incluido
        pos=antPath(length(antPath));%Actualiza la variable con el nodo actual de la hormiga
        if length(antPath)==length(model.x)+1 && antPath(length(antPath))==nest %Compara la el nodo de la hormiga con el de la comida
            exitosas=exitosas+1; %Si encuentra la comida suma +1 a las hormiga que han logrado llegar a la comida
            if exitosas==nants %Comprueba si todas las hormigas han llegado
                exito=true; %En de llegar todas se cambia el valor de la booleana exito
                
            end
        end 

    end    
end
end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UPDATE PHEROMONES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sumCoste=[];%Inicializamos el sumatorio de los costes
        for p=1:nants %Recorre todas las hormigas
            antPath=ant(p).Tabu;%Obtenemos el camino que ha llevado la hormiga
            coste=[ant(p).Cost];%Inicializamos el coste con el valor almacenado. Aquí se almacenarán los costes de cada salto para una única hormiga(inicializado en vacio)
            
            %For para almacenar la distancia entre todos los saltos
            for l=2:length(antPath)
                nodoi=antPath(l-1);%Obtenemos el nodo i
                nodoj=antPath(l);%Obtenemos el nodo j
                coste=[coste model.D(nodoi,nodoj)]; %Añadimos la distancia entre cada salto del camino seguido por la hormiga
                
            end 
            
            coste=sum(coste); %Sumamos las distancias del vector de una hormiga y obtenemos la distancia total recorrida
            ant(p).Cost=coste; %Almacenamos la distancia de una hormiga
            sumCoste=[sumCoste coste ];%Almacenamos en un vector aux el coste de todas las hormigas. Lo utilizaremos para encontrar el minimo valor entre todas las hormigas

        end
        
        betterant = find(sumCoste==min(sumCoste)); %Apunta a la hormiga con el camino mas corto
        betterpath = ant(betterant).Tabu; %Obtiene el camino de la hormiga con menor coste
        varFijatau=mrdivide(Q,ant(betterant(1)).Cost);%Calculo la variable que es constante. Apunta al índice 1 ya que puede ser que varias hormigas tengan el mismo coste, por ranto se elige la primera hormiga
        %x = mrdivide(B,A) forma alternativa de hacer div (utilizada porque
        %pensaba que la habitual fallaba :=) )
        for l=2:length(betterpath)%Recorre los indices del mejor camino
                nodoi=betterpath(l-1); %Nodosi
                nodoj=betterpath(l); %Nodosj
                tau(nodoi,nodoj)=tau(nodoi,nodoj)+varFijatau; %Actualizamos las distancias en la matriz de feromonas
                
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tau=(1-rho)*tau; %Evaporamos una parte de todos los caminos
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %GENERATE OUTPUT VARIABLE%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        BestSol.Tabu=betterpath;
        BestSol.Cost=min(sumCoste);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Show IterationInformation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['Iteration' num2str(k) ': Best Cost = ' num2str(BestSol.Cost)]);
        %figure 
        plotdistribution2(model,tau) %plot map and distribution of pheromone
        title(['Generated model. Iteration: ' num2str(k)],'FontWeight','bold')
        pause(0.2) %stop during0.1 seconds to allow the user to see the progression

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Show/Plot pheromone update
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end

%% PLOTDISTRIBUTION
function plotdistribution(model,tau)
%%%%%%%%%%%%%%%%% PLOT THE MAP %%%%%%%%%%%%%%%%%
clf%clearfigure
dt=model.dt; %recupera la triangulación
x=dt.Points(:,1); %coordenadas en x de los puntos
y=dt.Points(:,2); %coordenadas en y de los puntos 
tam=max([x; y])+1; % tamaño máximo del modelo
plot(x,y,'x') %pinta los nodos
triplot(dt,'r');  %# Plot the Delaunay triangulation
for i=1:length(dt.Points) %para cada nodo pinta el número de nodo
    text(dt.Points(i,1),dt.Points(i,2),num2str(i),'FontWeight','bold')
end
% axis([0 1 tam]) %fija los ejes
%%%%%%%%%%%%%%%%%%%%%%
%PLOT THE PHEROMONE%
%%%%%%%%%%%%%%%%%%%%%%
tau=tau./max(max(tau)); %normalizala pheromona
nVar=length(model.x); %numero de nodos
% El codigorecorre cada nodo y para cada vecino pinta una lineacuyo
% grosor depende de cuanta feremonahay en dicha unión
for n=1:nVar %para cada nodo
    nei=model.nei{n}; %vecinos 
    for nv=1:length(nei) % para cada vecinode cada node
        p1x=model.x(n); %coordenada x del nodo origen
        p1y=model.y(n); %coordenada y del nodo origen
        p2x=model.x(nei(nv)); %coordenada x del nodo destino
        p2y=model.y(nei(nv)); %coordenada y del nodo destino
        % lw=round(tau(n,nei(nv))*3); %nivel de feromona 
        lw=round(tau(n,nei(nv))*3); %nivel de feromona 
        if lw>0
            line([p1x p2x], [p1y p2y],'LineWidth',lw)      %pinta una lineacon groso igual a la pheromona
        end
    end
end
%Plot the number of each node
for i=1:length(dt.Points)
    text(dt.Points(i,1),dt.Points(i,2),num2str(i),'FontWeight','bold')
end
box on
xlabel('X coordenate(m)','FontWeight','bold')
ylabel('Y coordenate(m)','FontWeight','bold')
set(gca,'FontSize',12)
set(gca,'FontWeight','bold')
drawnow%refresca la imagen para asegura que se ve
end

%% VARIANTE PLOTDISTRIBUTION
function plotdistribution2(model,tau)
%%%%%%%%%%%%%%%%% PLOT THE MAP %%%%%%%%%%%%%%%%%
clf%clearfigure
dt=model.dt; %recupera la triangulación
x=dt.Points(:,1); %coordenadas en x de los puntos
y=dt.Points(:,2); %coordenadas en y de los puntos 
tam=max([x; y])+1; % tamaño máximo del modelo
plot(x,y,'x') %pinta los nodos
% triplot(dt,'r');  %# Plot the Delaunay triangulation
for i=1:length(dt.Points) %para cada nodo pinta el número de nodo
    text(dt.Points(i,1),dt.Points(i,2),num2str(i),'FontWeight','bold')
end
% axis([0 1 tam]) %fija los ejes
%%%%%%%%%%%%%%%%%%%%%%
%PLOT THE PHEROMONE%
%%%%%%%%%%%%%%%%%%%%%%
tau=tau./max(max(tau)); %normalizala pheromona
nVar=length(model.x); %numero de nodos
% El codigorecorre cada nodo y para cada vecino pinta una lineacuyo
% grosor depende de cuanta feremonahay en dicha unión
for n=1:nVar %para cada nodo
    nei=model.nei{n}; %vecinos 
    for nv=1:length(nei) % para cada vecinode cada node
        p1x=model.x(n); %coordenada x del nodo origen
        p1y=model.y(n); %coordenada y del nodo origen
        p2x=model.x(nei(nv)); %coordenada x del nodo destino
        p2y=model.y(nei(nv)); %coordenada y del nodo destino
        % lw=round(tau(n,nei(nv))*3); %nivel de feromona 
        lw=round(tau(n,nei(nv))*3); %nivel de feromona 
        if lw>0
            line([p1x p2x], [p1y p2y],'LineWidth',lw)      %pinta una lineacon groso igual a la pheromona
        end
    end
end
%Plot the number of each node
for i=1:length(dt.Points)
    text(dt.Points(i,1),dt.Points(i,2),num2str(i),'FontWeight','bold')
end
box on
xlabel('X coordenate(m)','FontWeight','bold')
ylabel('Y coordenate(m)','FontWeight','bold')
set(gca,'FontSize',12)
set(gca,'FontWeight','bold')
drawnow%refresca la imagen para asegura que se ve
end