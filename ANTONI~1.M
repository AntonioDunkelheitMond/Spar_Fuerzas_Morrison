clc
clear
close all

load nodos.list
DATOS = importdata('datos.inp'); %cargamos el archivo datos.inp
DATOS = table(DATOS.textdata, DATOS.data); %creamos una tabla con nuestros datos
DATOS1 = table2array(DATOS(:,2)); %extraemos la columna de la tabla y creamos una matriz


%% Extraemos los datos del documento datos

densidad=DATOS1(1,1); %Densidad, kg*10^3/m^3
d=DATOS1(2,1);        %Profundidad, m
g=DATOS1(3,1);        %Constante gravedad, m/s^2
T=DATOS1(4,1);        %Perido, s
w=DATOS1(5,1);        %Frecuencia, Hz
if T~=0&&w==0
    w=2*pi/T;
end
L=DATOS1(6,1);        %Longitud de onda, m
if L~=0
    k=2*pi/L;         %Número de onda rad/m
end
if w==0
    w=sqrt(g*k*tanh(k*d)); 
end
if T==0
    T=2*pi/w; 
end
if L==0 
    [L]=longitud_iterativo(g,T,d);
end
k=2*pi/L;             %Número de onda rad/m

H=DATOS1(7,1);         %altura de ola, m
Cd=DATOS1(8,1);
Cm=DATOS1(9,1);

%% Inicializamos los vectores donde alamacenaremos los resultados

m=size(nodos,1); %Obtenemos el número de elementos que tenemos a partir del 
                 %tamaño de la matriz nodos

Txyz_refuerzos=zeros(m,6); %En este vector extraeremos del archivo nodos las coordenas
                 %finales e iniciales de cada elemento
Diam=zeros(m,1); %En este vector alamacenaremos los valores del diámetro de 
                 %cada elemento
for i=1:m
    
    Txyz_refuerzos(i,:)=nodos(i,2:7);%Extraemos las coordenadas alamacenas en las columanas
    %2 a la 7, en donde de 2-4 se encuentran las coordenadas del
    %punto inicial, y de 5-7 se encuentran las coordenadas del último punto
    
    Diam(i,1)=nodos(i,8);  %y se utilizara para introducir la z, este vector 
             %son las coordenadas de los vectores de posicion final e inicial
end

%F=C*abs(Vn)*Vnx+K*anx, (kg/m^3)*(m^2)*(m/s^2)
Vnormal_total=zeros(3*m,6); %Es 3*m filas ya que cada elemento tendra 1 fila para el
                            %punto inicial, 1 fila para el punto
                            %intermedio, y una tercera fila para el punto
                            %final. Las columnas son vx,vy,z
anormal_total=zeros(3*m,6); %Igual que para las velocidades

vv=0;%variable auxiliar que me permitirá alamacenar las velocidades y aceleraciones
%normales en sus vectores

%Ahora inicializaremos tanto el vector donde alamacenaremos las fuerzas
%distribuidas, FNxm, en donde la primera columna corresponde a la fuerza
%distribuida en el primer punto del elemento, que corresponde al punto más
%"bajo" o mas "a popa", la segunda columna corresponderá al punto intermdio
%y la tercera columna corresponderá al punto final del elemento que se
%encuentra mas a proa o mas arriba del elmento en cuestion.
FNxm_brazing=zeros(m,3);%Vector que almacena las fuerzas distribuidas en x, 1ºcolumna
                %corresponde al primer punto, la 2º para el intermedio y
                %la 3º para el ultimo punto
FNzm_brazing=zeros(m,3); %Igual que para FNxm, pero para la dirección z
FNx_brazing=zeros(m,1); %Vector donde almacenaremos los valores totales de las fuerzss en x
FNz_brazing=zeros(m,1); %Vector donde almacenaremos los valores totales de las fuerzss en z
n_axial_refuerzos=zeros(m,4);     %Vector donde almacenaremos los vectores directores de cada elemento

%% Intervalo de tiempo de estudio

t_inicial=DATOS1(10,1); %Punto de inicio del calculo
t_final=DATOS1(11,1); %Punto final del calculo
nt=1;  %Numero de pasos donde haremos los calculos que se van a evaluar, es decir
       %hallaremos las fuerzas cada 2 segundos si es de 0-150s el rango.
nwt=((t_final-t_inicial)/nt)+1;%Se le añade el 1, para que me tenga en cuenta los extremos
F_refuerzosx_rangot=zeros(nwt,m);%En este vector se alamcenará los valores de las fuerzas
                 %en los distintos puntos de tiempo t. en donde las
                 %columnas corresponden a distintos elementos brazing y las
                 %filas a los distinto valoroes de tiempo t
F_refuerzosz_rangot=zeros(nwt,m);
                 
%% Fuerzas en los refuerzos, aplicamos la expresión de Morrison
t=t_inicial;
while t <= t_final
    
    wt=w*t; %Hallamos el valor del producto w*t
    
    for i=1:m %Recorremos todos los elementos desde el primer brazin i al ultimo m

        %Valores constantes 
        C=Cd*densidad*(Diam(i,1)/2);%constante que va delante de la fuerza de arrastre
        K=Cm*densidad*pi*((Diam(i,1)/2)^2); %cconstante que va delante de la fuerza de inercia

        %Vector direccional del elemento i. Lo obtendremos haciendo uso de las
        %coordenas iniciales y finales del elemento en cuestión
        ci=zeros(1,3);
        e1=Txyz_refuerzos(i,1:3); %punto inicial del elemento, (punto más bajo)
        e2=Txyz_refuerzos(i,4:6); %punto final del elemento, (punto más alto)
        ci(1,:)=(e2-e1)/norm(e2-e1); %vector unitario axial a la barra
        n_axial_refuerzos(i,1)=i; %Almacenamos a que elemento le pertenece este vector axial
        n_axial_refuerzos(i,2:4)=ci(1,:); %alamacenamos loa vectores directores
        long=norm(e2-e1);  %Tambien podemos hallar la longitud del elemento  

        %Obtenemos las coordenadas del punto intermedio
        eintermedio=(e1+e2)/2; %Nos da el punto intermedio del brazing
        %Vector coordenadas del elemento i
        xyz=zeros(1,9);
        xyz(1,1:3)=Txyz_refuerzos(i,1:3);%almacenamos las coordenadas del punto inicial
        xyz(1,4:6)=eintermedio(1,:);%almacenamos las coordenas del punto intermedio
        xyz(1,7:9)=Txyz_refuerzos(i,4:6);%almacenamos las coordenadas del punto final


        v=0;%variable intermedia para el bucle, que nos permitirá extraer las coordenas
        %x e y del archivo nodos, y que se reseteará de nuevo a 1 en cada nueva
        %iteración de i
        for j=1:3 %3 puntos de integración para la regla de simpson

            x=xyz(1,1+v);   %Extraemos la coordenada de x, del punto inicial, intermedio, final
            z=xyz(1,3+v); %Extraemos la coordenada de z, del punto inicial, intermedio, final

            %Dependiendo de si tenemos aguas profundas o aguas limitas se
            %utilizará una u otra expresión para la reducción de la amplitud
            %con la profundidad
            if d/L < 0.5  %Aguas transitorias,general water depth, poco profundas
                %La velocidad en la dirección x es:
                vx=(pi*H/T)*(cosh(k*(d+z))/sinh(k*d))*cos(k*x-wt);
                %La velocidad en la dirección z es:
                vz=(pi*H/T)*(sinh(k*(d+z))/sinh(k*d))*sin(k*x-wt);
            else %d/L > 0.5 1/2 aguas profundas
                %La velocidad en la dirección x es:
                vx=(pi*H/T)*(exp(k*z))*cos(k*x-wt);
                %La velocidad en la dirección z es:
                vz=(pi*H/T)*(exp(k*z))*sin(k*x-wt);
            end

            %Almacenamos el vector de velocidades en un vector columna
            vi_sin_normalizar=zeros(3,1);
            vi_sin_normalizar(1,1)=vx;
            vi_sin_normalizar(3,1)=vz;

            %Obtenemos el triple producto para la fuerza arrastre c1xVixc1, de
            %forma que obtenemos la velocidad normal al elemento
            Vnormal=cross(cross(ci,vi_sin_normalizar),ci);
            Vnormal=Vnormal';%Pasamos de un vector fila a un vector columna
            %Calculamos la magnitud del vector normal
            Vnormal_abs=norm(Vnormal); %es el vector normal de la velocidad nomal no la Vi

            %Dependiendo de si tenemos aguas profundas o aguas limitas se
            %utilizará una u otra expresión para la reducción de la amplitud
            %con la profundidad
            if d/L < 0.5  %Aguas transitorias,general water depth, poco profundas
                %La aceleración en la dirección x es:
                ax=(2*(pi^2)*H/(T^2))*(cosh(k*(d+z))/sinh(k*d))*sin(k*x-wt);
                %La aceleración en la dirección z es:
                az=(-2*(pi^2)*H/(T^2))*(sinh(k*(d+z))/sinh(k*d))*cos(k*x-wt);
            else %d/L > 0.5 1/2 aguas profundas
                %La aceleración en la dirección x es:
                ax=(2*(pi^2)*H/(T^2))*(exp(k*z))*sin(k*x-wt);
                %La aceleración en la dirección z es:
                az=(-2*(pi^2)*H/(T^2))*(exp(k*z))*cos(k*x-wt);
            end

            %Almacenamos el vector de aceleraciones en un vector columna
            ai_sin_normalizar=zeros(3,1);
            ai_sin_normalizar(1,1)=ax;
            ai_sin_normalizar(3,1)=az;

            %Obtenemos el triple producto para la fuerza inercia c1xaixc1, de
            %forma que obtenemos la velocidad normal al elemento
            anormal=cross(cross(ci,ai_sin_normalizar),ci); 
            anormal=anormal'; %como esta en fila la ponemos en columna

            %Calculamos las fuerzas de arrastre e inercia
            FD=C*Vnormal_abs*Vnormal; %Valor de la fuerza de arrastre en N/m
            FI=K*anormal; %Valor de la fuerza de inercia en N/m
            FNmi=FD+FI; %Valor de la fuerza total en N/m

            %Almacenamos la fuerza distribudia, la velocidad y la aceleracion normal
            FNmi=FNmi';
            FNxm_brazing(i,j)=FNmi(1,1); %Almacenamos las fuerzas en la dirección x
            FNzm_brazing(i,j)=FNmi(1,3); %Almacenamos las fuerzas en la dirección z
            Vnormal=Vnormal'; 
            Vnormal_total(j+vv,1)=i; %Almacenamos para que elemento es la fuerza
            Vnormal_total(j+vv,2)=t; %Almacenamos para que tiempo t se ha calculado
            Vnormal_total(j+vv,3)=j; %Almacenamos para que punto de integracion 1,2 o3
            Vnormal_total(j+vv,4:6)=Vnormal(1,:); %Almacenamos las velocidades de los
                                     %tres puntos de integracion en filas distintas
            anormal=anormal';
            anormal_total(j+vv,1)=i; %Almacenamos para que elemento es la fuerza
            anormal_total(j+vv,2)=t;%Almacenamos para que tiempo t se ha calculado
            anormal_total(j+vv,3)=j;%Almacenamos para que punto de integracion 1,2 o3
            anormal_total(j+vv,4:6)=anormal(1,:);%Almacenamos las aceleraciones de los
                                     %tres puntos de integracion en filas distintas
            
            v=v+3;%Nos movemos al siguiente set de coordenadas x,y,z   
        end
        
        vv=vv+3; %Nos movemos al siguiente set en el vector Vnormal_total y anormal_total
        
        %El espaciado de la integración es, para el caso de los brazing la
        %mitad de la longitud
        espaciado=(long/2);
        %Obtenemos la fuerza en x haciendo uso de la fuerza distribuida FNm,
        %para ello se aplicará la primera regla de simpson
        FNxi_brazing=espaciado*(1/3)*(1*FNxm_brazing(i,1)+4*FNxm_brazing(i,2)...
            +1*FNxm_brazing(i,3));
        FNx_brazing(i,1)=FNxi_brazing;
        
        FNzi_brazing=espaciado*(1/3)*(1*FNzm_brazing(i,1)+4*FNzm_brazing(i,2)...
            +1*FNzm_brazing(i,3));
        FNz_brazing(i,1)=FNzi_brazing;
        
        %Almacenamos las fuerzas de todos los elementos en las distintas
        %columna, i, para el mismo valor de tiempo. Ponemos t-(t_inicial-1)
        %de forma que si t_inicial no se parte de cero, no nos ponga las
        %demas casillas como cero y al final los resultados
        F_refuerzosx_rangot(t-(t_inicial-1),i)=FNx_brazing(i,1); 
        F_refuerzosz_rangot(t-(t_inicial-1),i)=FNz_brazing(i,1);
    end
    
    t=t+nt;%Pasamos al siguiente paso de tiempo
end

% Hallamos las fuerzas total de todos los refuerzos, por lo que hacemos el
% sumatorio de todas las fuerzas para el mismo valor de tiempo t
Ftotalx_refuerzos_rangot=zeros(nwt,1);
t=t_inicial;
while t <= t_final %Recorremos todos las filas temporales
    for i=1:m %Recorremos todos los elementos
        Ftotalx_refuerzos_rangot(t-(t_inicial-1),1)=...
            Ftotalx_refuerzos_rangot(t-(t_inicial-1),1)...
            +F_refuerzosx_rangot(t-(t_inicial-1),i);
    end
    t=t+nt;
end

plot(linspace(t_inicial,t_final,nwt),Ftotalx_refuerzos_rangot,'r-');
xlabel('Rango de tiempo, t (s)');
ylabel('Fuerzas total x, z pilares y refuerzos (N)');
title('Fuerzas total x, pilares y refuerzos (N) vs t (s)');
hold on

Ftotalz_refuerzos_rangot=zeros(nwt,1);
t=t_inicial;
while t <= t_final
    for i=1:m
        Ftotalz_refuerzos_rangot(t-(t_inicial-1),1)=...
            Ftotalz_refuerzos_rangot(t-(t_inicial-1),1)...
            +F_refuerzosz_rangot(t-(t_inicial-1),i);
    end
    t=t+nt;
end

plot(linspace(t_inicial,t_final,nwt),Ftotalz_refuerzos_rangot,'b-');
legend('Fuerza total x pilares y refuerzos','Fuerza total z pilares y refuerzos');
hold off

figure
%% Ahora vamos a calcular las fuerzas de froude krylov para las caras verticales del ponton o soft tank 

load ponton_vertical.list

m=size(ponton_vertical,1)/4; %Números de planchas, se divide por 4 ya que cada plancha tienen 4 puntos
                    %que definen sus vertices y que estan almacenados
                    %en la matriz ponton
n_directores_ponton=zeros(m,4); %vector donde almacenaremos los vectores directores
                                %de las planchas
Txyz_ponton=zeros(4*m,3);
Txyz_ponton(:,:)=ponton_vertical(:,2:4);
FNxm_ponton=zeros(m,3);
FNx_ponton=zeros(m,1);
L_plancha_ponton=zeros(m,1);
   
F_ponton_verticalx_rangot=zeros(nwt,m);%En este vector se alamcenará los valores de las fuerzas
                 %en los distintos puntos de tiempo t. en donde las
                 %columnas corresponden a distintos elementos plancha y las
                 %filas a los distinto valoroes de tiempo t         
t=t_inicial;
while t <= t_final
   
  wt=w*t; %Hallamos el valor del producto w*t
  v=0;    %variable auxiliar que se utilizará en el bucle para coger y saltar de una 
          %plancha a otra
    
    for i=1:m 
        
        %Extraemos los vertices de la planchas
        e1(1,:)=Txyz_ponton(1+v,:);
        e2(1,:)=Txyz_ponton(2+v,:);
        e3(1,:)=Txyz_ponton(3+v,:);
        e4(1,:)=Txyz_ponton(4+v,:);
        
        %Primero, haciendo uso de las vertices del ponton calcularemos los vectores directores                                
        %c1 y c2 son vectores unitarios contenidos en el plano de las planchas
        c1=(e1-e2)/norm(e1-e2); %vector unitario axial al lado corto de la plancha
        c2=(e4-e2)/norm(e4-e2); %vector unitario axial al lado largo de la plancha
        n_director=cross(c1,c2);%vector normal unitario a la superficie
        
        n_directores_ponton(i,1)=i; %Almacenamos el elemento 
        n_directores_ponton(i,2:4)=n_director(1,:); %Almacenamos sus vecto n correspondiente
        
        %Hallamos la longitdud de la plancha, del lado más largo para la
        %integral de froude krylov
        L_plancha=norm(e4-e2);
        L_plancha_ponton(i,1)=L_plancha; 
        
        %Obtenemos las coordenadas del punto intermedio de la plancha,
        %punto medio vertical para la integración numérica. Para ello
        %haremos uso de los puntos up and down en popa, que tienen iguales
        %coordenadas x e y, pero distinta altura z.
        eintermedio=(e1+e2)/2; %Nos da el punto intermedio 
        
        %Para poder hallar el angulo de desfase, kx-wt necesitamos hallar el valor
        %de x en los puntos de integración, que al ser una integración vertical,
        %se aproximará por la posción intermedia de x en la plancha, ya que al estar
        %inclinada es cierto que la ola no es igual en todos sus puntos.
        %Se hallara en el punto intermedio de la plancha longitudinalemente. 
        x=(e1(1,1)+e3(1,1))/2;
        %Creamos un vector zzz donde almacenaremos los puntos verticales de
        %integración, y que se utilizará en el siguiente loop
        zzz=zeros(1,3);
        zzz(1,1)=e1(1,3); %Punto inicial, que es el punto mas bajo de la plancha, integraremos
        zzz(1,2)=eintermedio(1,3);  %hacia arriba la fuerza FNx en la plancha
        zzz(1,3)=e2(1,3); %Punto final, punto mas alto de la plancha               
               
        for j=1:3 %Tres puntos de integración, j=1,2,3
            
            z=zzz(1,j);
            if d/L < 0.5 %Aguas transitorias,general water depth
                FNmi=(-densidad)*(-w*pi*H/(k*T))*L_plancha*(cosh(k*(d+z))/sinh(k*d))*sin(k*x-wt);
            else %d/L > 0.5 1/2 aguas profundas
                FNmi=(-densidad)*(-w*pi*H/(k*T))*L_plancha*exp(k*z)*sin(k*x-wt);
            end
            FNxm_ponton(i,j)=FNmi;%Almacenamos los valores de las fuerzas 
                                   %en los puntos de integración
        end
        
        %Para poder hacer la integración numérica, hay que hayar el
        %espaciado
        espaciado=abs(zzz(1,1)-zzz(1,2));
        %Obtenemos la fuerza en x haciendo uso de la fuerza distribuida FNm,
        %para ello se aplicará la primera regla de simpson
        FNxi_ponton=espaciado*(1/3)*(1*FNxm_ponton(i,1)+4*FNxm_ponton(i,2)+1*FNxm_ponton(i,3));
        FNx_ponton(i,1)=FNxi_ponton;
        
        %Almacenamos las fuerzas de todos los elementos en las distintas
        %columna, i, para el mismo valor de tiempo. Ponemos t-(t_inicial-1)
        %de forma que si t_inicial no se parte de cero, no nos ponga las
        %demas casillas como cero y al final los resultados
        F_ponton_verticalx_rangot(t-(t_inicial-1),i)=FNx_ponton(i,1); 
        
        v=v+4;%Pasamos al siguiente set de vertices de la siguiente plancha
    end
    
     t=t+nt;%Pasamos al siguiente paso de tiempo
end

%% Se van a calcular las fuerzas verticales para la cara superior e inferior horizontal del ponton

load ponton_horizontal.list


m=size(ponton_horizontal,1)/4; %Números de planchas, se divide por 4 ya que cada plancha tienen 4 puntos
                    %que definen sus vertices y que estan almacenados
                    %en la matriz ponton
Valores_logicos_ponton_horizontal=zeros(m,1); %PARA saber si es una plancha horizontal superior o inferior
%necesitamos decirselo al programa para que alterne el vector director de
%la superficie y esto se hara utilizando los valores lógicos introducidos
%en mismo archivo de los vertices, en la última columna, si tenemos 1,
%entonce se trata de una plancha superior y su vector normal apunta hacia
%arriba. si su valor es 0, es una plancha inferior y su vector normal
%apunta hacia abajo

Valores_logicos_ponton_horizontal=ponton_horizontal(:,5); %Estan alamcenados en al columna 5

n_directores_ponton_horizontal=zeros(m,4); %vector donde almacenaremos los vectores 
                                %directores de las planchas
Txyz_ponton_horizontal=zeros(4*m,3); %Matriz donde alamcenaremos los vertices
Txyz_ponton_horizontal(:,:)=ponton_horizontal(:,2:4); %Extraemos los vertices de la matriz
FNzm_ponton_horizontal=zeros(m,5); %Vector donde se alamcenan las fuerzas distribuidas,
                          %Por lo que en cada columna aparecerá el valor de
                          %la fuerza a integrar en un punto de integración,
                          %que tendremos 5 puntos.
FNz_ponton_horizontal=zeros(m,1);  %Vector donde se alamcena las fuerzas ya integradas para
                                   %un punto de tiempo t, y en donde cada
                                   %fila le corresponde a un elemento
   
F_ponton_horizontalz_rangot=zeros(nwt,m);%En este vector se alamcenará los valores de las fuerzas
                 %en los distintos puntos de tiempo t. en donde las
                 %columnas corresponden a distintos elementos plancha y las
                 %filas a los distinto valoroes de tiempo t         
t=t_inicial;
while t <= t_final
   
  wt=w*t; %Hallamos el valor del producto w*t
  v=0;    %variable auxiliar que se utilizará en el bucle para coger y saltar de una 
          %plancha a otra
  vv=1;   %Con esta variable saltaremos de 4 en 4 los coeficientes logicos
  for i=1:m 
        
        %Extraemos los vertices de la planchas
        e1(1,:)=Txyz_ponton_horizontal(1+v,:);
        e2(1,:)=Txyz_ponton_horizontal(2+v,:);
        e3(1,:)=Txyz_ponton_horizontal(3+v,:);
        e4(1,:)=Txyz_ponton_horizontal(4+v,:);
        
        %Primero, haciendo uso de las vertices del ponton calcularemos los vectores directores                                
        %c1 y c2 son vectores unitarios contenidos en el plano de las planchas
        c1=(e2-e1)/norm(e2-e1); %vector unitario axial al lado corto de la plancha
        c2=(e3-e1)/norm(e3-e1); %vector unitario axial al lado largo de la plancha
        n_director=cross(c1,c2);%vector normal unitario a la superficie
        
        %Ahora comprobamos si la plancha es la superior o la inferior
        if Valores_logicos_ponton_horizontal(vv,1)==0
            n_director=n_director*(-1);
        end
        
        n_directores_ponton_horizontal(i,1)=i; %Almacenamos elemento al que correponde
        n_directores_ponton_horizontal(i,2:4)=n_director(1,:); %Almacenamos el vector n
        
        %Obtenemos las coordenadas del punto intermedio de la plancha, para
        %ello haremos uso de las coordenas de los puntos extremos que se
        %encuentran en popa y en proa 
        eintermedio=(e2+e3)/2;
        %AUNQUE para que nos de resultados fiables hay que aumentar el
        %número de puntos de integración. Esto se ve claramente en el
        %ejemplo añadido para el cáculo del área para un rombo.
        
        %Para poder hallar el angulo de desfase, k*x-w*t, necesitamos hallar el valor
        %de x, que al ser ahora la planchahorizontal, no nos vale considerar una 
        %aproximación como si hicimos para la planchas verticales
        z=(e1(1,3)+e4(1,3))/2;
        %Extraemos las coordenas de x, ya que ahora vamos integrar en la
        %dirección de x y no en la dirección vertical de z
        xxx=zeros(1,5);
        xxx(1,1)=e3(1,1); %Punto inicial,corresponde al punto x de babor, punto inicial
        xintermedio1=(e1(1,1)+e3(1,1))/2;
        xxx(1,2)=xintermedio1; %coordenada entre el punto de popa x y el de babor x 
        xxx(1,3)=eintermedio(1,1);  %Coordenada del punto intermedio
        xintermedio2=(e1(1,1)+e2(1,1))/2;%coordenada x entre el punto x de estribor y el de popa x
        xxx(1,4)=xintermedio2;
        xxx(1,5)=e2(1,1); %Punto final, corresponde al punto x de estribor, punto final  
        
        
        %Ya que estamos integrando un rombo, la doble integral dydx, la vamos a
        %reducir a una sola integral dx, ya que ninguno de sus terminos
        %depende de y, por lo tenemos que obtener los valores de la
        %"altura" de la barra de integración, que va desde la "popa a proa"
        %de la estructura offshore, y nosotros nos movemos integrando de
        %"babor a estribor" en la direccion x. Ya que la ola viene de
        %"costado" segun nuestra nomenclatura.
        yyy=zeros(1,5);
        yyy(1,1)=e3(1,2); %Punto inicial,corresponde al punto y de babor, punto inicial
        y_para_xintermedio1=abs((e1(1,2)+e3(1,2))/2);
        yyy(1,2)=y_para_xintermedio1; %coordenada entre el punto de popa y y el de babor y
        yyy(1,3)=abs(e1(1,2)); %Coordenada y, del el punto intermedio o central del rombo
        y_para_xintermedio2=abs((e1(1,2)+e2(1,2))/2); 
        yyy(1,4)=y_para_xintermedio2;%coordenada y entre el punto y de estribor y el de popa y
        yyy(1,5)=e2(1,2); %Valor de y en la coordenada de estribor
               
        for j=1:5 %Cinco puntos de integración, j=1,2,3,4,5 para la integración de Simpson
            
            x=xxx(1,j);
            y=yyy(1,j);
            if d/L < 0.5 %Aguas transitorias,general water depth. 
                %Se multiplica por 2 el valor de las y porque, este valor
                %va desde el punto de popa al eje de referencia, y nosotros
                %queremos que vaya hasta el otro extremo de proa del rombo,
                %de forma que estamos haciendo de la integral de golpe
                FNmi=(-densidad)*(-w*pi*H/(k*T))*(2*y)*(cosh(k*(d+z))/sinh(k*d))*sin(k*x-wt);
            else %d/L > 0.5 1/2 aguas profundas
                FNmi=(-densidad)*(-w*pi*H/(k*T))*(2*y)*exp(k*z)*sin(k*x-wt);
            end
            FNzm_ponton_horizontal(i,j)=FNmi;%Almacenamos los valores de las fuerzas 
                                   %en los puntos de integración
        end
        
        %Para poder hacer la integración numérica, hay que hayar el
        %espaciado. Para ello hallamos la distancia entre el punto inicial
        %más a popa y el punto intermedio
        espaciado=abs(xxx(1,1)-xxx(1,2));
        %Obtenemos la fuerza en x haciendo uso de la fuerza distribuida FNm,
        %para ello se aplicará la primera regla de simpson
        FNzi_ponton_horizontal=espaciado*(1/3)*(1*FNzm_ponton_horizontal(i,1)+...
            4*FNzm_ponton_horizontal(i,2)+2*FNzm_ponton_horizontal(i,3)+...
            4*FNzm_ponton_horizontal(i,4)+1*FNzm_ponton_horizontal(i,5));
        FNz_ponton_horizontal(i,1)=FNzi_ponton_horizontal;
        
        %Almacenamos las fuerzas de todos los elementos en las distintas
        %columna, i, para el mismo valor de tiempo. Ponemos t-(t_inicial-1)
        %de forma que si t_inicial no se parte de cero, no nos ponga las
        %demas casillas como cero y al final los resultados
        F_ponton_horizontalz_rangot(t-(t_inicial-1),i)=FNz_ponton_horizontal(i,1); 
        
        v=v+4;%Pasamos al siguiente set de vertices de la siguiente plancha
        vv=vv+4; %Pasamos al siguiente set de valores lógicos
   end
    
     t=t+nt;%Pasamos al siguiente paso de tiempo
end

% Hallamos las fuerzas total de todos el ponton, por lo que hacemos el
% sumatorio de todas las fuerzas para el mismo valor de tiempo t
Ftotalx_ponton_verticalx_rangot=zeros(nwt,1);
t=t_inicial;
while t <= t_final
    for i=1:m
        Ftotalx_ponton_verticalx_rangot(t-(t_inicial-1),1)=...
            Ftotalx_ponton_verticalx_rangot(t-(t_inicial-1),1)...
            +F_ponton_verticalx_rangot(t-(t_inicial-1),i);
    end
    t=t+nt;
end

plot(linspace(t_inicial,t_final,nwt),Ftotalx_ponton_verticalx_rangot,'r-');
xlabel('Rango de tiempo, t (s)');
ylabel('Fuerzas total x,zde soft tank (N)');
title('Fuerzas total x, z de soft tank (N) vs t (s)');

hold on

Ftotalz_ponton_horizontalz_rangot=zeros(nwt,1);
t=t_inicial;
while t <= t_final
    for i=1:m
        Ftotalz_ponton_horizontalz_rangot(t-(t_inicial-1),1)=...
            Ftotalz_ponton_horizontalz_rangot(t-(t_inicial-1),1)...
            +F_ponton_horizontalz_rangot(t-(t_inicial-1),i);
    end
    t=t+nt;
end

plot(linspace(t_inicial,t_final,nwt),Ftotalz_ponton_horizontalz_rangot,'b-');
legend('Fuerza total x Soft tank','Fuerza total z Soft tank');

hold off

figure
%% Ahora vamos a calcular las fuerzas de froude krylov para las planchas verticales de arfada

load plancha_arfada_vertical.list

m=size(plancha_arfada_vertical,1)/4; %Números de planchas, se divide por 4 ya 
                %que cada plancha tienen 4 puntos que definen sus vertices
                %y que estan almacenados en la matriz ponton
n_directores_plancha_arfada_vertical=zeros(m,4); %vector donde almacenaremos los vectores directores
                                %de las planchas
Txyz_arfada=zeros(4*m,3); %Matriz donde almacenaremos las coordenas de los vertices de las
                          % planchas de arfada
Txyz_arfada(:,:)=plancha_arfada_vertical(:,2:4); %Extraemos los vertices

FNxm_arfada=zeros(m,3);%Aquí almacenaremos los valores de las fuerzas en los 
                       %puntos de integracion de Simpson, 1º regla
FNx_arfada=zeros(m,1); %Aquí alamacenaremos los valores de las fuerzas finales por plancha
L_plancha_arfada=zeros(m,1); %Aquí almacenaremos las longitudes de integración

F_plancha_arfada_verticalx_rangot=zeros(nwt,m);%En este vector se alamcenará los valores de las fuerzas
                 %en los distintos puntos de tiempo t. en donde las
                 %columnas corresponden a distintos elementos plancha y las
                 %filas a los distinto valoroes de tiempo t
t=t_inicial;
while t <= t_final
   
  wt=w*t; %Hallamos el valor del producto w*t
  v=0;    %variable auxiliar que se utilizará en el bucle para coger y saltar de una 
          %plancha a otra
    
    for i=1:m 
        
        %Extraemos los vertices de la planchas
        e1(1,:)=Txyz_arfada(1+v,:);
        e2(1,:)=Txyz_arfada(2+v,:);
        e3(1,:)=Txyz_arfada(3+v,:);
        e4(1,:)=Txyz_arfada(4+v,:);
        
        %Primero, haciendo uso de las vertices del ponton calcularemos los vectores directores                                
        %c1 y c2 son vectores unitarios contenidos en el plano de las planchas
        c1=(e1-e2)/norm(e1-e2); %vector unitario axial al lado corto de la plancha
        c2=(e4-e2)/norm(e4-e2); %vector unitario axial al lado largo de la plancha
        n_director=cross(c1,c2);%vector normal unitario a la superficie
        
        n_directores_plancha_arfada_vertical(i,1)=i; %Almacenamos el elemento
        n_directores_plancha_arfada_vertical(i,2:4)=n_director(1,:); %Almacenamos su vector director
        
        %Hallamos la longitdud de la plancha, del lado más largo para la
        %integral de froude krylov
        L_plancha=norm(e4-e2);
        L_plancha_arfada(i,1)=L_plancha; 
        
        %Obtenemos las coordenadas del punto intermedio de la plancha
        eintermedio=(e1+e2)/2; %Nos da el punto intermedio 
        
        %Para poder hallar el angulo de desfase necesitamos hallar el valor
        %de x, que se tomará el valor medio entre e1 y e2, que se encontrara en
        %el punto medio de la plancha. 
        x=(e1(1,1)+e3(1,1))/2;
        zzz=zeros(1,3);
        zzz(1,1)=e1(1,3); %Punto inicial, que es el punto mas bajo de la plancha, integraremos
        zzz(1,2)=eintermedio(1,3);  %hacia arriba la fuerza FNx en la plancha
        zzz(1,3)=e2(1,3); %Punto final, punto mas alto de la plancha               
               
        for j=1:3 %Tres puntos de integración, j=1,2,3
            
            z=zzz(1,j);
            if d/L < 0.5 %Aguas transitorias,general water depth
                FNmi=(-densidad)*(-w*pi*H/(k*T))*L_plancha*(cosh(k*(d+z))/sinh(k*d))*sin(k*x-wt);
            else %d/L > 0.5 1/2 aguas profundas
                FNmi=(-densidad)*(-w*pi*H/(k*T))*L_plancha*exp(k*z)*sin(k*x-wt);
            end
            FNxm_arfada(i,j)=FNmi;%Almacenamos los valores de las fuerzas 
                                  %en los puntos de integración
        end
        
        %Para poder hacer la integración numérica, hay que hayar el
        %espaciado
        espaciado=abs(zzz(1,1)-zzz(1,2));
        %Obtenemos la fuerza en x haciendo uso de la fuerza distribuida FNm,
        %para ello se aplicará la primera regla de simpson
        FNxi_arfada=espaciado*(1/3)*(1*FNxm_arfada(i,1)+4*FNxm_arfada(i,2)+1*FNxm_arfada(i,3));
        FNx_arfada(i,1)=FNxi_arfada;
        
        F_plancha_arfada_verticalx_rangot(t-(t_inicial-1),i)=FNx_arfada(i,1);
        
        v=v+4;%Pasamos al siguiente set de vertices de la siguiente plancha
    end
    
    t=t+nt;%Pasamos al siguiente paso de tiempo
end

%% Se van a calcular las fuerzas verticales para la cara superior horizontal del ponton

load plancha_arfada_horizontal.list

m=size(plancha_arfada_horizontal,1)/4; %Números de planchas, se divide por 4 ya que cada plancha tienen 4 puntos
                    %que definen sus vertices y que estan almacenados
                    %en la matriz ponton
Valores_logicos_plancha_arfada_horizontal=zeros(m,1); %PARA saber si es una plancha horizontal superior o inferior
%necesitamos decirselo al programa para que alterne el vector director de
%la superficie y esto se hara utilizando los valores lógicos introducidos
%en mismo archivo de los vertices, en la última columna, si tenemos 1,
%entonce se trata de una plancha superior y su vector normal apunta hacia
%arriba. si su valor es 0, es una plancha inferior y su vector normal
%apunta hacia abajo

Valores_logicos_plancha_arfada_horizontal=plancha_arfada_horizontal(:,5); %Estan alamcenados en al columna 5

n_directores_plancha_arfada_horizontal=zeros(m,4); %vector donde almacenaremos los vectores directores
                                %de las planchas
Txyz_plancha_arfada_horizontal=zeros(4*m,3);
Txyz_plancha_arfada_horizontal(:,:)=plancha_arfada_horizontal(:,2:4);

FNzm_plancha_arfada_horizontal=zeros(m,5); %matriz donde almacenaremos las fuerzas
%distribuidas, para cada uno de los 5 puntos de integración
FNz_plancha_arfada_horizontal=zeros(m,1);  %Vector columna donde almacenaremos las
%fuerzas una vez integradas
   
F_plancha_arfada_horizontalz_rangot=zeros(nwt,m);%En este vector se alamcenará los valores de las fuerzas
                 %en los distintos puntos de tiempo t. en donde las
                 %columnas corresponden a distintos elementos plancha y las
                 %filas a los distinto valoroes de tiempo t         
t=t_inicial;
while t <= t_final
   
  wt=w*t; %Hallamos el valor del producto w*t
  v=0;    %variable auxiliar que se utilizará en el bucle para coger y saltar de una 
          %plancha a otra
  vv=1;   %Con esta variable saltaremos de 4 en 4 los coeficientes logicos
  for i=1:m
        %Extraemos los vertices de la planchas
        e1(1,:)=Txyz_plancha_arfada_horizontal(1+v,:);
        e2(1,:)=Txyz_plancha_arfada_horizontal(2+v,:);
        e3(1,:)=Txyz_plancha_arfada_horizontal(3+v,:);
        e4(1,:)=Txyz_plancha_arfada_horizontal(4+v,:);
        
        %Primero, haciendo uso de las vertices del ponton calcularemos los vectores directores                                
        %c1 y c2 son vectores unitarios contenidos en el plano de las planchas
        c1=(e2-e1)/norm(e2-e1); %vector unitario axial al lado corto de la plancha
        c2=(e3-e1)/norm(e3-e1); %vector unitario axial al lado largo de la plancha
        n_director=cross(c1,c2);%vector normal unitario a la superficie
        
        %Ahora comprobamos si la plancha es la superior o la inferio
        if Valores_logicos_plancha_arfada_horizontal(vv,1)==0
            n_director=n_director*(-1);
        end
        
        n_directores_plancha_arfada_horizontal(i,1)=t; %Almacenamos el elemento
        n_directores_plancha_arfada_horizontal(i,2:4)=n_director(1,:); %Almacenamos n
        
        %Obtenemos las coordenadas del punto intermedio de la plancha, para
        %ello haremos uso de las coordenas de los puntos extremos que se
        %encuentran en popa y en proa 
        eintermedio=(e2+e3)/2;
        %AUNQUE para que nos de resultados fiables hay que aumentar el
        %número de puntos de integración. Esto se ve claramente en el
        %ejemplo añadido para el cáculo del área para un rombo.
        
        %Para poder hallar el angulo de desfase, k*x-w*t, necesitamos hallar el valor
        %de x, que al ser ahora la planchahorizontal, no nos vale considerar una 
        %aproximación como si hicimos para la planchas verticales
        z=(e1(1,3)+e4(1,3))/2;
        %Extraemos las coordenas de x, ya que ahora vamos integrar en la
        %dirección de x y no en la dirección vertical de z
        xxx=zeros(1,5);
        xxx(1,1)=e3(1,1); %Punto inicial,corresponde al punto x de babor, punto inicial
        xintermedio1=(e1(1,1)+e3(1,1))/2;
        xxx(1,2)=xintermedio1; %coordenada entre el punto de popa x y el de babor x 
        xxx(1,3)=eintermedio(1,1);  %Coordenada del punto intermedio
        xintermedio2=(e1(1,1)+e2(1,1))/2;%coordenada x entre el punto x de estribor y el de popa x
        xxx(1,4)=xintermedio2;
        xxx(1,5)=e2(1,1); %Punto final, corresponde al punto x de estribor, punto final  
        
        
        %Ya que estamos integrando un rombo, la doble integral dydx, la vamos a
        %reducir a una sola integral dx, ya que ninguno de sus terminos
        %depende de y, por lo tenemos que obtener los valores de la
        %"altura" de la barra de integración, que va desde la "popa a proa"
        %de la estructura offshore, y nosotros nos movemos integrando de
        %"babor a estribor" en la direccion x. Ya que la ola viene de
        %"costado" segun nuestra nomenclatura.
        yyy=zeros(1,5);
        yyy(1,1)=e3(1,2); %Punto inicial,corresponde al punto y de babor, punto inicial
        y_para_xintermedio1=abs((e1(1,2)+e3(1,2))/2);
        yyy(1,2)=y_para_xintermedio1; %coordenada entre el punto de popa y y el de babor y
        yyy(1,3)=abs(e1(1,2)); %Coordenada y, del el punto intermedio o central del rombo
        y_para_xintermedio2=abs((e1(1,2)+e2(1,2))/2); 
        yyy(1,4)=y_para_xintermedio2;%coordenada y entre el punto y de estribor y el de popa y
        yyy(1,5)=e2(1,2); %Valor de y en la coordenada de estribor
               
        for j=1:5 %Cinco puntos de integración, j=1,2,3,4,5 para la integración de Simpson
            
            x=xxx(1,j);
            y=yyy(1,j);
            if d/L < 0.5 %Aguas transitorias,general water depth. 
                %Se multiplica por 2 el valor de las y porque, este valor
                %va desde el punto de popa al eje de referencia, y nosotros
                %queremos que vaya hasta el otro extremo de proa del rombo,
                %de forma que estamos haciendo de la integral de golpe
                FNmi=(-densidad)*(-w*pi*H/(k*T))*(2*y)*(cosh(k*(d+z))/sinh(k*d))*sin(k*x-wt);
            else %d/L > 0.5 1/2 aguas profundas
                FNmi=(-densidad)*(-w*pi*H/(k*T))*(2*y)*exp(k*z)*sin(k*x-wt);
            end
            FNzm_plancha_arfada_horizontal(i,j)=FNmi;%Almacenamos los valores de las fuerzas 
                                   %en los puntos de integración
        end
        
        %Para poder hacer la integración numérica, hay que hayar el
        %espaciado. Para ello hallamos la distancia entre el punto inicial
        %más a popa y el punto intermedio
        espaciado=abs(xxx(1,1)-xxx(1,2));
        %Obtenemos la fuerza en x haciendo uso de la fuerza distribuida FNm,
        %para ello se aplicará la primera regla de simpson
        FNzi_plancha_arfada_horizontal=espaciado*(1/3)*(1*FNzm_plancha_arfada_horizontal(i,1)+...
            4*FNzm_plancha_arfada_horizontal(i,2)+2*FNzm_plancha_arfada_horizontal(i,3)+...
            4*FNzm_plancha_arfada_horizontal(i,4)+1*FNzm_plancha_arfada_horizontal(i,5));
        FNz_plancha_arfada_horizontal(i,1)=FNzi_plancha_arfada_horizontal;
        
        %Almacenamos las fuerzas de todos los elementos en las distintas
        %columna, i, para el mismo valor de tiempo. Ponemos t-(t_inicial-1)
        %de forma que si t_inicial no se parte de cero, no nos ponga las
        %demas casillas como cero y al final los resultados
        F_plancha_arfada_horizontalz_rangot(t-(t_inicial-1),i)=FNz_plancha_arfada_horizontal(i,1); 
        
        v=v+4;%Pasamos al siguiente set de vertices de la siguiente plancha
        vv=vv+4; %Pasamos al siguiente set de valores lógicos
   end
    
     t=t+nt;%Pasamos al siguiente paso de tiempo
end

% Hallamos las fuerzas total de todas las planchas de arfada por lo que hacemos el
% sumatorio de todas las fuerzas para el mismo valor de tiempo t
Ftotalx_plancha_arfada_verticalx_rangot=zeros(nwt,1);
t=t_inicial;
while t <= t_final
    for i=1:m
        Ftotalx_plancha_arfada_verticalx_rangot(t-(t_inicial-1),1)=...
            Ftotalx_plancha_arfada_verticalx_rangot(t-(t_inicial-1),1)...
            +F_plancha_arfada_verticalx_rangot(t-(t_inicial-1),i);
    end
    t=t+nt;
end

plot(linspace(t_inicial,t_final,nwt),Ftotalx_plancha_arfada_verticalx_rangot,'r-');
xlabel('Rango de tiempo, t (s)');
ylabel('Fuerzas total x de planchas arfada (N)');
title('Fuerzas total x,z de planchas de arfada (N) vs t (s)');

hold on

Ftotalz_plancha_arfada_horizontalz_rangot=zeros(nwt,1);
t=t_inicial;
while t <= t_final
    for i=1:m
        Ftotalz_plancha_arfada_horizontalz_rangot(t-(t_inicial-1),1)=...
            Ftotalz_plancha_arfada_horizontalz_rangot(t-(t_inicial-1),1)...
            +F_plancha_arfada_horizontalz_rangot(t-(t_inicial-1),i);
    end
    t=t+nt;
end

plot(linspace(t_inicial,t_final,nwt),Ftotalz_plancha_arfada_horizontalz_rangot,'b-');
legend('Fuerza total x Planchas arfada','Fuerza total z Planchas arfada');

hold off

figure
%% Calculamos ahora la fuerza ejercida en el cilindro, hard tank sumergido segun la teoria de la difracción

%La fuerza en y, y en z son cero, ya que la ola se propaga en x, y en la
%direcion transversal y, se compensa.
m=21; %Numero de puntos de integracion para la primera regla de simpson
zsup=DATOS1(12,1);            %Coordenada superior del cilindro
zinf=DATOS1(13,1);            %Coordenada inferior del cilindro
espaciado=abs((zsup-zinf)/m); %Hallamos el espaciado de la integración
D_cilindro=DATOS1(14,1);      %Extraemos el diamtro del cilindro

kr=k*D_cilindro/2; %Hallamos el producto kr

FNx_cilindro=0;         %Vector final de la fuerza en la dirección x

%Calculamos los siguientes valores constante de las funciones de Bessel
J0=besselj(0,kr);
J1=besselj(1,kr);
Y0=bessely(0,kr);
Y1=bessely(1,kr);
J_prima=J0-(1/kr)*J1;
Y_prima=Y0-(1/kr)*Y1;
A=1/sqrt(J_prima^2 + Y_prima^2);
alpha=atan(J_prima/Y_prima);

F_cilindrox_rangot=zeros(nwt,1);%En este vector se alamcenará los valores de las fuerzas
                 %en los distintos puntos de tiempo t. en donde las
                 %columnas corresponden a distintos elementos cilindro y las
                 %filas a los distinto valoroes de tiempo t
t=t_inicial;
while t <= t_final
 
    wt=w*t; %Hallamos el valor del producto w*t
    
    z=0; %Reseteamos el valor de z para cada paso de tiempo
    Sumatorio=0; %Reseteamos tamvien el sumatorio
    for i=1:m  %Recorremos todos los puntos de integración
        FNxi=(2*densidad*g*H/k)*(cosh(k*(d+z))/cosh(k*d))*A*cos(wt-alpha);
        %Multiplicamos la fuerza por el integrando correspondiente
        if i==1 || i==m
            FNxi=FNxi*1;        
        elseif rem(i,2)==0 %Numero par
            FNxi=FNxi*4;    
        elseif rem(i,2)==1 %Numero impar
            FNxi=FNxi*2;
        end
        Sumatorio=Sumatorio+FNxi;
        z=z-espaciado; %Pasamos al siguiente punto de integración

    end
    
      FNx_cilindro=(1/3)*espaciado*Sumatorio;
      
      F_cilindrox_rangot(t-(t_inicial-1),1)=FNx_cilindro(1,1); %Solo tenemos 1 cilindro por lo que (1,1)
            
      t=t+nt;%Pasamos al siguiente paso de tiempo
end

Ftotalx_cilindro_rangot=zeros(nwt,1);
t=t_inicial;
while t <= t_final
    for i=1:1
        Ftotalx_cilindro_rangot(t-(t_inicial-1),1)=...
            Ftotalx_cilindro_rangot(t-(t_inicial-1),1)...
            +F_cilindrox_rangot(t-(t_inicial-1),i);
    end
    t=t+nt;
end

%Representamos la fuerza que sufre el cilindro
plot(linspace(t_inicial,t_final,nwt),Ftotalx_cilindro_rangot,'g-');
xlabel('Rango de tiempo, t (s)');
ylabel('Fuerzas en el Hard tank total x (N)');
title('Fuerzas en el Hard tank total x (N) vs t (s)');

figure

%% Ahora vamos a hallar el valor total de las fuerzas en el eje x y eje z

%Calculamos la fuerza total en el eje x
Ftotalx=Ftotalx_refuerzos_rangot+Ftotalx_ponton_verticalx_rangot+...
        Ftotalx_plancha_arfada_verticalx_rangot+Ftotalx_cilindro_rangot;
   
plot(linspace(t_inicial,t_final,nwt),Ftotalx,'r-');
xlabel('Rango de tiempo, t (s)');
ylabel('Fuerzas total x (N)');
title('Fuerzas total x (N) vs t (s)');

figure

%Calculamos la fuerza total en el eje z
Ftotalz=Ftotalz_refuerzos_rangot+Ftotalz_ponton_horizontalz_rangot +...
        Ftotalz_plancha_arfada_horizontalz_rangot;
    
plot(linspace(t_inicial,t_final,nwt),Ftotalz,'b-');
xlabel('Rango de tiempo, t (s)');
ylabel('Fuerzas total z (N)');
title('Fuerzas total z (N) vs t (s)');

figure

%A parte de representar las fuerzas x,z individualmente, se va a realizar
%la representación de ambas fuerzas de forma que se pueda ver más
%claramente su diferencia en magnitud

plot(linspace(t_inicial,t_final,nwt),Ftotalx,'r-');
xlabel('Rango de tiempo, t (s)');
ylabel('Fuerzas total en (N)');
title('Fuerzas total x,z (N) vs t (s)');

hold on
plot(linspace(t_inicial,t_final,nwt),Ftotalz,'b-');
legend('Fuerza total x','Fuerza total z');
hold off


clear FNmi ci g i j K C FD FI long v vx vz ax az nwt...
    ai_sin_normalizar vi_sin_normalizar L_plancha n_director x z xyz m Vnormal_abs...
    Cm Cd c1 c2 anormal Vnormal espaciado zzz
clear e1 e2 e3 e4 eintermedio Txyz_arfada Txyz_ponton Txyz_refuerzos J0 J1...
    Y0 Y1 Y_prima J_prima A alpha kr DATOS1 Diam vv Sumatorio  ...
    D_cilindro 
clear wt t xintermedio1 yintermedio1 xxx y y_para_xintermedio1 y_para_xintermedio2...
    yyy zinf zsup xintermedio2 yintermedio2
clear FNxi_brazing FNxi_ponton FNxi_arfada FNzi_brazing FNzi_ponton_horizontal...
    FNzi_plancha_arfada_horizontal FNx_cilindro ...
    FNz_plancha_arfada_horizontal FNzm_plancha_arfada_horizontal...
    FNzm_ponton_horizontal FNz_ponton_horizontal...
    FNx_arfada FNxm_arfada FNx_brazing FNxm_brazing FNz_brazing FNzm_brazing...
    FNx_ponton FNxm_ponton FNxi...
    Valores_logicos_plancha_arfada_horizontal...
    Valores_logicos_ponton_horizontal
clear Txyz_plancha_arfada_horizontal Txyz_ponton_horizontal L_plancha_arfada...
    L_plancha_ponton

%close all
    
    
    
    
    
    
    
    
    
