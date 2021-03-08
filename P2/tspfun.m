function [cost] = tspfun(pop,model)
%% Calcular el coste de cada individuo y lo almacena en una matrix
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

