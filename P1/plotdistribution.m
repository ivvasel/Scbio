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