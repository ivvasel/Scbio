function myperceptron_AND_53384654
% Funcion principal que realiza las funciones de
 %1) Creación de las variables para el banco de entrenamiento y banco de
    %validación
 %2) Creación de la red neuronal (perceptron)
 %3) Entrenamiento de la red neuronal con el set de valores del banco de
    %entrenamiento
 %4) Validación de la red con el banco de validación.
 %5) Calculo y representación del error cometido.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creamos entradas y salidas de entrenamiento para una funcion AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
muestras_train = 5000;
E1 = round(rand(muestras_train,1));
E2 = round(rand(muestras_train,1));
SE = double(and(E1,E2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muestras_val = 20;
E1V = round(rand(muestras_val,1));
E2V = round(rand(muestras_val,1));
SEV = double(and(E1V,E2V));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos un perceptron para 2 entradas %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myperceptron=initialize_perceptron(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entrenamos el perceptron para un LR, por defecto 0.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LR=0.7;
myperceptronT=train_perceptron(myperceptron,LR,[E1 E2],SE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluamos el perceptron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_est=useperceptron(myperceptronT,[E1V E2V]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculamos y representamos el error cometido
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error=mean(abs(SEV-S_est));
 
close all
figure,
plot(SEV,'ok','LineWidth',2),hold on
plot(S_est,'xr','LineWidth',2)
%set(gcf,'FontWeight','bold')
set(gca,'FontSize',12) %# Fix font size of the text in the current axes 
set(gca,'FontWeight','bold')  %# Fix Bold text in the current axes 
xlabel('Number of test','FontWeight','bold')
ylabel('Output values','FontWeight','bold')
axis([-1 length(SEV)+1 -0.1 1.3])
legend('Correct Values','Perceptron Output')
title('Evaluation of Perceptron Ouput for an AND','FontWeight','bold')
 
function [myperceptron]=initialize_perceptron(n_inputs)
%function [myperceptron]=initialize_perceptron(n_inputs)
%funcion que inicializa un perceptron y le asigna pesos aleatorios Esta función crea e inicializa la estructura myperceptron.weights donde se guardan los pesos del perceptron. Los valores de los pesos iniciales deben ser números aleatorios entre -1 y 1.
% INPUTS:
    % n_inputs: numero de entradas al perceptron
% OUPUT:
      % myperceptron:  estructura con el perceptron 
          % myperceptron.bias:       %Bias del perceptron (e.g. 1)
          % myperceptron.weights:    %Pesos del perceptron
          
%Creamos el struct del perceptron
%n_inputs = 10;
myperceptron = struct;  

%Creamos el array de pesos
%inicializamos los pesos a aleatorios
myperceptron.weights = -1+2*rand(1,n_inputs+1);
myperceptron.bias = [1];
end
 
function [myperceptron]=train_perceptron(myperceptron,LR,input,output)
% function myperceptron=train_perceptron(myperceptron,LR,input,output)
% funcion que modifica los pesos del perceptron para que vaya aprendiendo a
% a partir de los valores de entrada que se le indican
% ESTE PERCEPTRON UTILIZA:
    % Funcion sigma SIGMOIDAL
    % Entrenamiento DELTA RULE
% INPUTS:
   % myperceptron:  estructura con el perceptron
           %myperceptron.bias:       %Bias del perceptron (e.g. 1)
           %myperceptron.weights:    %Pesos del perceptron 
   % LR: learning rate (e.g. 0.7)
   % input: matriz con valores de entrada de entrenamiento (e.g. [E1 E2])
   % output: vector con valores de salida de entrenamiento (e.g. [SE])
% OUPUT:
      % myperceptron:  estructura con el perceptron ya entrenado
          %myperceptron.bias:       %Bias del perceptron (e.g. 1)
          %myperceptron.weights:    %Pesos del perceptron ya entrenado
%and_training = and(E1,E2);
%and_validation = and(E1V,E2V);
%a = 1; %Parámetro "a" de la función sigmoidal

    for i=1:length(input(:,1))
        x = [myperceptron.bias input(i,:)]; %Creamos el vector de inputs
        producto = x.*myperceptron.weights;
        local_field = sum(producto); %Sumatorio del producto entre las entradas y los pesos
        y = 1/(1+exp(-local_field)); %Salida del perceptrón
        out_error = output(i)-y; %Definimos el error como la resta de la salida actual y la salida deseada
        myperceptron.weights = myperceptron.weights + (LR.*out_error.*x); %Actualizamos los pesos
    end
  
end
 
function [out]=useperceptron(myperceptron,input)
% function out=useperceptron(myperceptron,input)
% funcion que utiliza el perceptron para calcular las salidas a partir de
% las entradas de acuerdo con lo que haya aprendido el perceptron en la
% fase de entrenamiento
%a=1;
    for i=1:length(input(:,1))
        x = [myperceptron.bias input(i,:)]; %Creamos el vector de inputs
        local_field = sum(x.*myperceptron.weights); %Sumatorio del producto entre las entradas y los pesos
        out(i) = 1./(1+exp(-local_field)); %Salida del perceptrón 
    end


end


end