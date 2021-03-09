% function 1.	tspga_DNI (fmode,model)
% fucntion tspga_DNI(fmode,model)
% Funci�n que utiliza algoritmos gen�ticos para resolver el problema del
% viajante (Travel Sales Problem TSP)
% INPUTS
    % mode: indica el modo de generaci�n los nodos del modelo
    %       est� funci�n llamar� a la funci�n generamodelo_DNI creada en la
    %       primera pr�ctica
    %    modo=0; % los nodos se genera aleatoriamente
    %    modo=1; % los nodos los indica el usuario manualmente
    %    modo=2; % utilizar los nodos almacenados
    %    mode=3; % la estructura del modelo se introduce como una variable
    % model: estructura generada previamente con generamodelo_DNI.m
% OUPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERACION DEL MODELO EN FUNCI�N DE LAS INSTRUCCIONES DEL USUARIO %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin<1
     fmode=3;
% end
numberofnodes=20;
tam=10;
if fmode==3
%     if nargin<2
%          disp('Error: El modo 3 requiere introducir el modelo')
%     end
end
if fmode<3
    model=generamodelo_45913106(numberofnodes,tam,fmode);
end
x=model.x;
y=model.y;
numberofnodes=length(x);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETROS DEL ALGORITMO GEN�TICO %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npar=numberofnodes; % # of optimization variables
Nt=npar; % # of columns in population matrix
maxit=10000; % max number of iterations
popsize=20; % set population size / miembros de la poblaci�n
mutrate=.05; % set mutation rate
selection=0.5; % fraction of population kept / fracci�n de miembros sobreviven
iga=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INICIALIZA LA POBLACI�N Y VECTOR DE COSTES%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pop= ones(popsize,npar);
for i=1:popsize
    pop(i,:)=randperm(npar);
end

cost=tspfun(pop,model); %funci�n que eval�a el fitness de cada miembro de la poblacion
[cost,ind]=sort(cost); %devuelve el vector ordenado y los indices del antiguo vector.
pop=pop(ind,:); %reordena la poblaci�n de acuerdo al coste



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EN CADA ITERACI�N VAMOS A IR GUARDANDO EL VALOR DE COSTE DEL MEJOR INDIVIDUO Y DE LA POBLACI�N MEDIA PARA POSTERIORMENTE PODER VER EL PROGRESO%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 minc(1)=min(cost); % minc contains min of population
 meanc(1)=mean(cost); % meanc contains mean of population
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERO EL VECTOR DE PROBABILIDADES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g. si tengo una poblaci�n de 10 individuos y solo la mitad sobreviven
% para reproducirse entonces el vector de probabilidades ser�a.
% Prob=[5 4 4 3 3 3 2 2 2 2 1 1 1 1 1];
% Seleccionamos a la mejor mitad de la poblaci�n.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PONEMOS EL MUNDO VIRTUAL A FUNCIONAR (MAIN LOOP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while iga<maxit
    iga=iga+1; % increments generation counter
    
    pop=pop([1:popsize*selection],:);
        odds=[];
        k=popsize*selection;
    for i=1:popsize*selection
        for j=1:k
            odds=[i odds];        
        end
            k=k-1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ELIJO A LOS PADRES Y MADRES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % De entre el vector de probabilidades elijo aleatoriamente
    % a los padres y a las madres.
        N=round(popsize*selection/2);
        indicespadres=ceil(length(odds)*rand(1,N));
        indicesmadres=ceil(length(odds)*rand(1,N));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Emparejamientos y generaci�n de los hijos %
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
legend('Coste M�nimo', 'Coste Medio')
box on
set(gca,'FontSize',12)
set(gca,'FontWeight','bold')
plotsolution(x,y,pop)
 
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
        coste=[coste model.D(nodoi,nodoj)]; %A�adimos la distancia entre cada salto del camino seguido por la hormiga
                
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
tam=max([x; y])+1; % tama�o m�ximo del modelo
plot(x,y,'x') %pinta los nodos
axis([0 tam 0 tam]) %fija los ejes
 
%%%%%%%%%%%%%%%%%%%%%
% PINTA LA SOLUCI�N %
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

