%% function 1.	tspga_DNI (fmode,model)
% fucntion tspga_DNI(fmode,model)
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
% if nargin<1
     fmode=0;
% end
numberofnodes=8;
tam=10;
% if fmode==3
%     if nargin<2
%          disp('Error: El modo 3 requiere introducir el modelo')
%     end
% end
% if fmode<3
    model=generamodelo_45913106(numberofnodes,tam,fmode);
% end
x=model.x;
y=model.y;
numberofnodes=length(x);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETROS DEL ALGORITMO GENÉTICO %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npar=numberofnodes; % # of optimization variables
Nt=npar; % # of columns in population matrix
maxit=2000; % max number of iterations
popsize=20; % set population size / miembros de la población
mutrate=.05; % set mutation rate
selection=0.5; % fraction of population kept / fracción de miembros sobreviven
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
pop=pop([1:popsize*selection],:);
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
% while iga<maxit
%     iga=iga+1; % increments generation counter
    
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
   
%     for i=1:N %Bucle de emparejamiento e hijos
        
        %Obtenemos quien son los padres
        padre=odds(indicespadres(1)); 
        madre=odds(indicesmadres(1));
        %Cogemos los genes de la matriz pop
        genespadre = pop(padre,:);
        genesmadre = pop(madre,:);
        
       
        hijo1=zeros(1,npar);
        hijo2=zeros(1,npar);
        
        Nn=npar/2;
        at=randi(npar);
        counter1=Nn;
        auxmadre=genesmadre;
        i=at+Nn;
        
        %Hijo 1
        %Genes del padre
        while(counter1~=0)
            i = i +1;
            if i > npar
                i=1;
            end  
            hijo1(i)=genespadre(i);
            counter1=counter1 - 1;
            auxhij=hijo1(i); %Obtiene 
            indi=find(auxmadre==hijo1(i)); %Indice del gen utilizado por el otro padre
            auxmadre(:,indi)=[]; %Contiene los genes restantes que no han sido utilizado por el otro padre
            
        end
        
        j=1;
        counter2=Nn;
     
        %Genes de la madre
        while(counter2~=0)
            i=i+1;
            if i > npar
                i=1;                              
            end
            hijo1(i)=auxmadre(j); 
            counter2=counter2 - 1;
            
            j=j+1;
        end
        
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mutate the population
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTRODUCE EL TU CÓDIGO
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Evalua el coste de la nueva población
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Sort the costs and associated parameters
%     [cost,ind]=sort(cost);
% pop=pop(ind,:);
%     
%     _______________________________________________________
%     Do statistics
%     minc(iga)=min(cost);
%     meanc(iga)=mean(cost);
%         
%     iga
% end %iga
 

%_______________________________________________________
% % Displays the output
% day=clock;
% disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5),day(6)),0))
% %disp(['optimized function is ' ff])
% format short g
% disp(['popsize=' num2str(popsize) 'mutrate=' num2str(mutrate) '# par=' num2str(npar)])
% disp([' best cost=' num2str(cost(1))])
% disp(['best solution']); disp([num2str(pop(1,:))])
%  
% figure(2)
% iters=1:maxit;
% plot(iters,minc,iters,meanc,'--','LineWidth',2);
% hold on
% xlabel('generation','FontWeight','bold','FontSize',12);
% ylabel('cost','FontWeight','bold','FontSize',12);
% legend('Coste Mínimo', 'Coste Medio')
% box on
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')
% plotsolution(x,y,pop)
%  
% function dist=tspfun(pop,model)
% % Cost function for the traveling salesman problem.
% 
% INTRODUCE EL TU CÓDIGO
% 
% function plotsolution(x,y,pop)
%  
% %%%%%%%%%%%%%%%%%%%
% % PINTA LOS NODOS %
% %%%%%%%%%%%%%%%%%%%
% figure(1);
% clf %clear figure
% tam=max([x; y])+1; % tamaño máximo del modelo
% plot(x,y,'x') %pinta los nodos
% axis([0 tam 0 tam]) %fija los ejes
%  
% %%%%%%%%%%%%%%%%%%%%%
% % PINTA LA SOLUCIÓN %
% %%%%%%%%%%%%%%%%%%%%%
% hold on
% sol=pop(1,:);
% for i=1:length(sol)-1
%     p1x=x(sol(i));
%     p2x=x(sol(i+1));
%     p1y=y(sol(i));
%     p2y=y(sol(i+1));
%     line([p1x p2x], [p1y p2y],'LineWidth',2)         
% end
%     p1x=x(sol(1));
%     p2x=x(sol(end));
%     p1y=y(sol(1));
%     p2y=y(sol(end));
%     line([p1x p2x], [p1y p2y],'LineWidth',2)      
% axis square
%  
% hold on
% %Plot the number of each node
% for i=1:length(x)    
%     text(x(i),y(i),num2str(i),'FontWeight','bold','FontSize',12)
% end
% box on
% xlabel('X coordenate (m)','FontWeight','bold')
% ylabel('Y coordenate (m)','FontWeight','bold')
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')
