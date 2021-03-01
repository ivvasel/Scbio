function [BestSol] = aco_TSP_45913106(model,nants,niter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF EMPTY INPUT VARIABLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<1,model=generamodelo_45913106(30,1,0); end %generate a random model
    if nargin<3,niter=100; end
    if nargin<2,nants=40; end      
nest=1;
% food=30;
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
   ant(n).Tabu=nest;%Inicializamos cada hormiga en el hormiguero
   ant(n).Cost=[];%Inicializamos el coste 
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct a solution for each ANT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exitosas=0;
    exito = false;
while(~exito)
for n=1:nants %Recorre todas las hormigas    
    antPath=ant(n).Tabu; %Inicializamos el vector con el camino que lleva la hormiga
    pos=antPath(length(antPath)); %Posición actual obtenida del vector del camino   
    
    %%Comprueba si se han recorrido todos los nodos y salta a la 
    %%siguiente hormiga
    if length(antPath)==length(model.x)
        continue     
    
    %%No ha llegado a la comida y sigue buscando
    else
        auxNei=model.nei{pos};
        checkNei=auxNei;
        allP=[];
        %%calcula tau y eta de la posicion con los posibles vecinos 
        %%que no se hayan visitado ya. 
        %%Calculamos tmb probabilidad de los posibles casos por separado
        for p=1:length(auxNei);  
        
        
        posiblesalto=auxNei(p); %Obtiene un posible salto de model.nei{p}
        comprobar=find(antPath==posiblesalto); %Comprobar si posiblesalto ya se encuentra en el camino de la ant
        if isempty(comprobar)
            %Calculamos p hacia el posiblesalto
            auxTau=tau(pos,posiblesalto);
            auxEta=eta(pos,posiblesalto);            
            allP=[allP (auxTau^alpha)*(auxEta^beta)];        
        
        else
            elimnei=find(checkNei==posiblesalto);
            checkNei(elimnei,:)=[];   
            
        end 
        
        end
    
    
    
        %%Sumamos la probabilidad de los posibles casos
        %calcula tau y eta de la posicion con los posibles vecinos que no se hayan visitado ya
        sumaP = sum(allP);
        P=[];
    
    
        %%Calculamos la probabilidad para cada caso y creamos el vector P
        for p=1:length(allP);
        P=[P allP(p)/sumaP];
        end
    
    
        %%Elegimos de forma euristica el nodo al que iremos
        r=rand;
        C=cumsum(P);
        j=find(r<=C,1,'first'); %selecciona el primero mayor que r
        ant(n).Tabu=[antPath checkNei(j)];%introduce en el path el nodo seleccionado
    
        %%Comprueba si no hay vecinos a los que ir y reinicia el path
        if isempty(checkNei)
            antPath=[nest];
            ant(n).Tabu=antPath;            
        end
        
        %%Comprueba de nuevo si se ha recorrido todos los nodos y se suma a
        %%hormigas exitosas
        antPath=ant(n).Tabu;
        pos=antPath(length(antPath));
        if length(antPath)==length(model.x)
            exitosas=exitosas+1;
            %Añade el coste de la ruta a la matriz
            %ant(n).Cost=length(antPath);
            if exitosas==nants
                exito=true;
                
            end
        end 

    end    
end
end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UPDATE PHEROMONES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sumCoste=[];
        for p=1:nants
            antPath=ant(p).Tabu;
            coste=[ant(p).Cost];%Inicializamos el coste con el valor almacenado
            
            %For para almacenar la distancia entre todos los saltos
            for l=2:length(antPath)
                nodoi=antPath(l-1);
                nodoj=antPath(l);
                coste=[coste model.D(nodoi,nodoj)]; %Añadimos las distancias al vector
                
            end 
            coste=sum(coste); %Sumamos las distancias del vector de una hormiga
            ant(p).Cost=coste; %Almacenamos la distancia de una hormiga
            sumCoste=[sumCoste coste ];%Almacenamos en un vector aux el coste de todas las hormigas
            
            
        end
        
        betterant = find(sumCoste==min(sumCoste));%Apunta a la hormiga con el camino mas corto
        betterpath = ant(betterant).Tabu;
        varFijatau=mrdivide(Q,ant(betterant(1)).Cost);
        %x = mrdivide(B,A) forma alternativa de hacer div
        for l=2:length(betterpath)
                nodoi=betterpath(l-1);
                nodoj=betterpath(l);
                tau(nodoi,nodoj)=tau(nodoi,nodoj)+varFijatau; %Añadimos las distancias al vector
                
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tau=(1-rho)*tau;
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Show/Plot pheromone update
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plotdistribution(model,tau);
end
 
    
    

end

