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

