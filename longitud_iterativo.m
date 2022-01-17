
function [L]=longitud_iterativo(g,T,d)
%e=error=(g*(T^2)/(2*pi))*tanh(2*pi*d/L)-L

e=zeros(40,2); %Inicializamos la matriz donde almacenaremos los errores
%para poder despues buscar en el el valor mínimo
i=1;
for L=1:10:400 %olas de 1 a 400 m de eslora en saltos de 10 m
    error=(g*(T^2)/(2*pi))*tanh(2*pi*d/L)-L;
    e(i,1)=L;%Almacenamos la longitud L del error e
    e(i,2)=error; %Almacenamos el error para la longitud L de onda
    i=i+1;
end
 emin=min(abs(e(:,2))); %Hallamos el valor más pequeño
%Ahora recorremos la segunda columna de errores de la matriz e donde  
%buscamos la longitud L que le corresponde al error emin
for i=1:30
    if emin == abs(e(i,2))
       Lmin=e(i,1);
    end
end
%Ahora realizamos exactamente lo mismo que ante pero con mayor precisión
i=1;
Linf=Lmin-10;
Lsup=Lmin+10;
n=(Lsup-Linf)/40;
for L=Linf:n:Lsup %olas de 1 a 400 m de eslora en saltos de (Lsup-Linf)/40
    error=(g*(T^2)/(2*pi))*tanh(2*pi*d/L)-L;
    e(i,1)=L;
    e(i,2)=error;
    i=i+1;
end
 emin=min(abs(e(:,2)));
for i=1:30
    if emin == abs(e(i,2))
       Lmin=e(i,1);
    end
end
L=Lmin; %Renonbramos la longitud de onda L al valor minimo Lmin
end