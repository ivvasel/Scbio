%%PRACTICA 4: NEURAL NETWORKS 
% Realizada por: Iv醤 Valencia y Ra鷏 L髉ez.
% Este archivo contiene las 4 subfunciones:
% 1. myperceptron_AND
% 2. myperceptron_XOR
% 3. mynn_XOR
% 4. matlabperceptron

%% Interfaz de usuario
disp('Elige una funci髇: ')
disp('1: Funcion AND')
disp('2: Funcion XOR')
disp('3: XOR multicapa')
disp('4: Funcion XOR con toolbox de matlab')

funcion=input('Elecci髇: ');

switch funcion
    case 1
        myperceptron_AND_45913106;
        
    case 2
        myperceptron_XOR_45913106;
        
    case 3
        mynn_XOR_45913106();
    
    case 4
        matlabperceptron_45913106;
        
end


function myperceptron_AND_45913106
% Funcion principal que realiza las funciones de
 %1) Creaci贸n de las variables para el banco de entrenamiento y banco de
    %validaci贸n
 %2) Creaci贸n de la red neuronal (perceptron)
 %3) Entrenamiento de la red neuronal con el set de valores del banco de
    %entrenamiento
 %4) Validaci贸n de la red con el banco de validaci贸n.
 %5) Calculo y representaci贸n del error cometido.
    
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
%funcion que inicializa un perceptron y le asigna pesos aleatorios Esta funci贸n crea e inicializa la estructura myperceptron.weights donde se guardan los pesos del perceptron. Los valores de los pesos iniciales deben ser n煤meros aleatorios entre -1 y 1.
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
%a = 1; %Par谩metro "a" de la funci贸n sigmoidal

    for i=1:length(input(:,1))
        x = [myperceptron.bias input(i,:)]; %Creamos el vector de inputs
        producto = x.*myperceptron.weights;
        local_field = sum(producto); %Sumatorio del producto entre las entradas y los pesos
        y = 1/(1+exp(-local_field)); %Salida del perceptr贸n
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
        out(i) = 1./(1+exp(-local_field)); %Salida del perceptr贸n 
    end


end


end

function myperceptron_XOR_45913106
% Funcion principal que realiza las funciones de
 %1) Creaci贸n de las variables para el banco de entrenamiento y banco de
    %validaci贸n
 %2) Creaci贸n de la red neuronal (perceptron)
 %3) Entrenamiento de la red neuronal con el set de valores del banco de
    %entrenamiento
 %4) Validaci贸n de la red con el banco de validaci贸n.
 %5) Calculo y representaci贸n del error cometido.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creamos entradas y salidas de entrenamiento para una funcion AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
muestras_train = 5000;
E1 = round(rand(muestras_train,1));
E2 = round(rand(muestras_train,1));
SE = double(xor(E1,E2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muestras_val = 20;
E1V = round(rand(muestras_val,1));
E2V = round(rand(muestras_val,1));
SEV = double(xor(E1V,E2V));
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
title('Evaluation of Perceptron Ouput for an XOR','FontWeight','bold')
 
function [myperceptron]=initialize_perceptron(n_inputs)
%function [myperceptron]=initialize_perceptron(n_inputs)
%funcion que inicializa un perceptron y le asigna pesos aleatorios Esta funci贸n crea e inicializa la estructura myperceptron.weights donde se guardan los pesos del perceptron. Los valores de los pesos iniciales deben ser n煤meros aleatorios entre -1 y 1.
% INPUTS:
    % n_inputs: numero de entradas al perceptron
% OUPUT:
      % myperceptron:  estructura con el perceptron 
          % myperceptron.bias:       %Bias del perceptron (e.g. 1)
          % myperceptron.weights:    %Pesos del perceptron
          
%Creamos el struct del perceptron
%n_inputs = 10;
myperceptron = struct;  

%Creamos un array de pesos con el mismo n煤mero de elementos que n_inputs

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
a = 1; %Par谩metro "a" de la funci贸n sigmoidal

    for i=1:length(input(:,1))
        x = [myperceptron.bias input(i,:)]; %Creamos el vector de inputs
        producto = x.*myperceptron.weights;
        local_field = sum(producto); %Sumatorio del producto entre las entradas y los pesos
        y = 1/(1+exp(-a*local_field)); %Salida del perceptr贸n
        out_error = output(i)-y; %Definimos el error como la resta de la salida actual y la salida deseada
        myperceptron.weights = myperceptron.weights + LR.*out_error.*x; %Actualizamos los pesos
    end
  
end
 
function [out]=useperceptron(myperceptron,input)
% function out=useperceptron(myperceptron,input)
% funcion que utiliza el perceptron para calcular las salidas a partir de
% las entradas de acuerdo con lo que haya aprendido el perceptron en la
% fase de entrenamiento
a=1;
    for i=1:length(input(:,1))
        x = [myperceptron.bias input(i,:)]; %Creamos el vector de inputs
        local_field = sum(x.*myperceptron.weights); %Sumatorio del producto entre las entradas y los pesos
        out(i) = 1/(1+exp(-a*local_field)); %Salida del perceptr贸n 
    end


end


end

function mynn_XOR_45913106(VT)
%function mynn_XOR_profesor_AND_profesor
% Funcion principal que realiza las funciones de
 %1) Creaci贸n de las variables para el banco de entrenamiento y banco de
    %validaci贸n
 %2) Creaci贸n de la red neuronal de DOS CAPAS
 %3) Entrenamiento de la red neuronal con el set de valores del banco de
    %entrenamiento
 %4) Validaci贸n de la red con el banco de validaci贸n.
 %5) Calculo y representaci贸n del error cometido.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creamos entradas y salidas de entrenamiento para una funcion AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
if nargin < 1   
    VT=50000;                %numero de muestras de entrada
end
NV=round(VT);           %aseguramos que VT sea un numero entero
E1=round(rand(NV,1));   %vector de NV valores de entrada en pin 1
E2=round(rand(NV,1));   %vector de NV valores de entrada en pin 2
SE=double(xor(E1,E2));   %vector salida ideal para las entradas E1 y E

%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
NVV=20;
E1V=round(rand(NVV,1));
E2V=round(rand(NVV,1));
SEV=double(xor(E1V,E2V));
 
% Inicializamos un perceptron para 2 entradas %%%
mynn=initialize_nn;
% Entrenamos el perceptron para un LR, por defecto 0.7
LR=0.7;
mynnT=train_nn(mynn,LR,[E1 E2],SE);
%evaluamos el perceptron
S_est=usenn(mynnT,[E1V E2V]);

%Calculamos y representamos el error cometido
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
title('Evaluation of MULTILAYER Ouput for an XOR','FontWeight','bold')
 
function [mynn]=initialize_nn
%function [myperceptron]=initialize_perceptron(n_inputs)
%funcion que inicializa un perceptron y le asigna pesos aleatorios
%funcion que inicializa una estructura de perceptron
% INPUTS:
    % n_inputs: numero de entradas al perceptron
% OUTPUTS:
    % myperceptron:  estructura con el perceptron
           %myperceptron.bias:       %Bias del perceptron (e.g. 1)
           %myperceptron.weights:    %Pesos del perceptron (habr谩 tantos
                                     %como indice n_inputs +1 )
rand('state',sum(100*clock));  %inicializa random en funci贸n del reloj
mynn.bias = [1 1 1]; %vector con los 3 bias
mynn.weights = -1+2.*(rand(3,3)); %Pesos aleatorios
end
 
function [mynnT]=train_nn(mynn,LR,input,output)
% funcion que modifica los pesos la red para que vaya aprendiendo 
% ESTE PERCEPTRON UTILIZA:
    % Funcion sigma SIGMOIDAL
    % Entrenamiento BACKPROPAGATION
% INPUTS:
   % mynn:  estructura con el perceptron
           %mynnT.bias:       %Bias del perceptron (e.g. 1)
           %mynnT.weights:    %Pesos del perceptron 
   % LR: learning rate (e.g. 0.7)
   % input: matriz con valores de entrada de entrenamiento (e.g. [E1 E2])
   % output: vector con valores de salida de entrenamiento (e.g. [SE])
% OUPUT:
      % mynnT:  estructura con el perceptron ya entrenado
          %mynnT.bias:       %Bias del perceptron (e.g. 1)
          %mynnT.weights:    %Pesos del perceptron ya entrenado  
[nsamples,ninputs]=size(input);
weights=mynn.weights;
bias=mynn.bias;
coeff = LR;
out=zeros(nsamples,1)';

for i=1:nsamples
    x = input(i,:);
    
    h1 = x(1).*weights(1,2)+x(2).*weights(1,3)+weights(1,1).*bias(1,1);
    h2 = x(2).*weights(2,3)+x(1).*weights(2,2)+weights(2,1).*bias(1,2);
    
    y1 = 1./(1+exp(-h1));
    y2 = 1./(1+exp(-h2));
    
    h3 = y1.*weights(3,2)+ y2.*weights(3,3)+weights(3,1).*bias(1,3);
    
    y3 = 1./(1+exp(-h3));
    
    delta3 = y3.*(1-y3).*(output(i)-y3); %Delta para la 3a neurona
    delta1 = y1.*(1-y1).*weights(3,2).*delta3; %Delta para las neuronas 1 y 2
    delta2 = y2.*(1-y2).*weights(3,3).*delta3;
    
    weights(1,:) = weights(1,:)+(coeff*delta1.*[bias(1,1) x]); %Actualizamos pesos en la 1a neurona
    weights(2,:) = weights(2,:)+(coeff*delta2.*[bias(1,2) x]); %Actualizamos pesos en la 2a neurona
    weights(3,:) = weights(3,:)+(coeff*delta3.*[bias(1,3) y1 y2]); %Actualizamos pesos en la 3a neurona
end

mynnT.bias = bias;
mynnT.weights = weights; %Par谩metros de la NN entrenada

%INCLUIR C脫DIGO
end
function [out]=usenn(mynn,input)
% function out=useperceptron(myperceptron,input)
% funcion que utiliza el perceptron para calcular las salidas a partir de
% las entradas de acuerdo con lo que haya aprendido el perceptron en la
% fase de entrenamiento
[nsamples,ninputs]=size(input);
weights=mynn.weights;
bias=mynn.bias;
out=zeros(nsamples,1)'; 
for i=1:nsamples
    x = input(i,:);
    
    h1 = x(1).*weights(1,2)+x(2).*weights(1,3)+weights(1,1).*bias(1,1);
    h2 = x(2).*weights(2,3)+x(1).*weights(2,2)+weights(2,1).*bias(1,2);
    
    y1 = 1./(1+exp(-h1));
    y2 = 1./(1+exp(-h2));
    
    h3 = y1.*weights(3,2)+ y2.*weights(3,3)+weights(3,1).*bias(1,3);
    
    out(i) = 1./(1+exp(-h3));
    
end

end
end

function matlabperceptron_45913106
% Funcion principal que realiza las funciones de
 %1) Creaci贸n de las variables para el banco de entrenamiento y banco de
    %validaci贸n
 %2) Creaci贸n de la red neuronal (perceptron)
 %3) Entrenamiento de la red neuronal con el set de valores del banco de
    %entrenamiento
 %4) Validaci贸n de la red con el banco de validaci贸n.
 %5) Calculo y representaci贸n del error cometido.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creamos entradas y salidas de entrenamiento para una funcion AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
muestras_train = 1000;
E1 = round(rand(muestras_train,1));
E2 = round(rand(muestras_train,1));
SE = double(xor(E1,E2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muestras_val = 10;
E1V = round(rand(muestras_val,1));
E2V = round(rand(muestras_val,1));
SEV = double(xor(E1V,E2V));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos un perceptron para 2 entradas %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
net=patternnet(2,'trainrp');
net=configure(net,[E1 E2]',SE');
net=init(net);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entrenamos el perceptron para un LR, por defecto 0.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
net.divideParam.trainRatio = 0.7;
netT=train(net,[E1 E2]',SE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluamos el perceptron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_est=netT([E1V E2V]');

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
title('Evaluation of Perceptron Ouput for an XOR','FontWeight','bold')
end

