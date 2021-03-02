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

    
