% LAB 3 - MSIM
% Autor: Bernardo Rocha & Vasco Ara�jo
% N�mero de Aluno: 89867 & 90817
% Maio 2020; �ltima Revis�o: 13/05/2020   

%% Exerc�cio 5

clear
close all
clc


%inicializa��o de variaveis
g = 9.8;
k = 3;
L = 0.5;
M = 0.15;
l = 0.4;
m = 0.2;
beta = 0.1;
J = m*l^2+(M*L^2)/3;

%condi�oes iniciais
teta0 = 0;
d_teta0 = pi/4;

stop_time = 10;

sim5 = sim('system5');

figure(1)
plot(sim5.tout, sim5.teta);
grid on;
xlabel('Tempo [s]');
ylabel('teta(t)');
title('Evolu��o do �ngulo de teta(t)');

figure(2)
plot(sim5.tout, sim5.d_teta);
grid on;
xlabel('Tempo [s]');
ylabel('d teta(t)');
title('Evolu��o da velocidade de teta(t)');

figure(3)
plot(sim5.teta, sim5.d_teta);
xlabel('teta(t)');
ylabel('d teta(t)');
title('Tra�ado de fase');

%% Modelo SIMULINK 5

system5


%% Exerc�cio 6

clear
close all
clc


%inicializa��o de variaveis
g = 9.8;
k = 3;
L = 0.5;
M = 0.15;
l = 0.4;
m = 0.2;
beta = 0.1;
J = m*l^2+(M*L^2)/3;

stop_time = 10;


%definir as matrizes de estado
A =  [          0                    1      ;
     ((-k+g*(m*l+(M*L/2)))/J)    (-beta/J) ];

B =  [        0          ;         1/J     ];

C =  [        1                     0       ;
              0                     1      ];

D =  [        0          ;          0      ];
  
x0 = [        0          ;        pi/4     ];  % condi�oes iniciais

T =  [        0                     0      ];

sim6 = sim('system6');

figure(1)
plot(sim6.tout, sim6.y(:,1));
grid on;
xlabel('Tempo [s]');
ylabel('teta(t)');
title('Evolu��o do �ngulo de teta(t)');

figure(2)
plot(sim6.tout, sim6.y(:,2));
grid on;
xlabel('Tempo [s]');
ylabel('d teta(t)');
title('Evolu��o da velocidade de teta(t)');

figure(3)
plot(sim6.y(:,1), sim6.y(:,2));
xlabel('teta(t)');
ylabel('d teta(t)');
title('Tra�ado de fase');

%%
% Teoricamente obtem-se que a matriz C � dada por C = [ 1 0 ], uma vez que
% se tem y = C $\cdot$ x. No entanto, para efeitos de simula��o a matriz C
% tem de ser 2x2 devido ao bloco State-Space do SIMULINK, o que at� acaba
% por ser conveniente pois fazendo C = [ 1 0 ; 0 1 ] podemos observar o
% comportamento de $x_2$, que neste caso representa a velocidade de
% $\theta$(t) de forma a fazermos o tra�ado de fase do sistema.

%% Modelo SIMULINK 6

system6

%% Exerc�cio 7

clear
close all
clc


%inicializa��o de variaveis
g = 9.8;
k = 3;
L = 0.5;
M = 0.15;
l = 0.4;
m = 0.2;

J = m*l^2+(M*L^2)/3;

stop_time = 10;


%definir as matrizes de estado
B =  [        0          ;         1/J     ];

C =  [        1                     0       ;
              0                     1      ];

D =  [        0          ;          0      ];
  
x0 = [        0          ;        pi/4     ];  % condi�oes iniciais

T =  [        0                     0      ];

for i=0:1:1
    beta = i;
    
    A =  [          0                    1      ;
     ((-k+g*(m*l+(M*L/2)))/J)    (-beta/J) ];
 
    sim7 = sim('system6');
    
    figure(1)
    plot(sim7.tout, sim7.y(:,1));
    grid on
    hold on
    
    figure(2)
    plot(sim7.tout, sim7.y(:,2));
    grid on
    hold on
    
    figure(3)
    plot(sim7.y(:,1), sim7.y(:,2));
    hold on
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%respostas no plano de estado de diferentes condi�oes inicias
 
for beta = 0:1:1
    %recalcular matriz A
    A =  [          0                    1      ;
        ((-k+g*(m*l+(M*L/2)))/J)    (-beta/J) ];
    for p = -5:2:5
        %mudar condi��o inicial
        x0 = [p   p]; 
        sim7_p = sim('system6');
        if beta == 0
            figure(4)
        else figure(5)
        end
        plot(sim7_p.y(:,1), sim7_p.y(:,2));
        hold on
        grid on
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %faz o campo dos vectores (quiver)
    X = linspace(-max(abs(sim7_p.y(:,1)))-1, max(abs(sim7_p.y(:,1)))+1, 15);
    Y = linspace(-max(abs(sim7_p.y(:,2)))-1, max(abs(sim7_p.y(:,2)))+1, 15);

    for i = 1:size(X,2)
        for j = 1:size(Y,2) 
            q = A*[X(i);Y(j)];
            u(j,i) = q(1);
            v(j,i) = q(2);
        end
    end
    quiver(X, Y, u, v,0.7);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculo dos valores e vetores proprios
    [vectors, values] = eig(A);
    fprintf('\n\nPara beta = %d:\n Vectores pr�prios:\n', beta);
    disp(vectors);
    fprintf('\n Valores pr�prios:\n');
    disp(values);
end

%calcular valores e vetores proprios de A para o valor de beta da questao 5
beta = 0.1;
A =  [          0                    1      ;
      ((-k+g*(m*l+(M*L/2)))/J)    (-beta/J) ]; 
[vectors, values] = eig(A);
fprintf('\n\nPara beta = %.1f:\n Vectores pr�prios:\n', beta);
disp(vectors);
fprintf('\n Valores pr�prios:\n');
disp(values);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%graficos

figure(1)
xlabel('Tempo [s]');
ylabel('teta(t)');
title('Evolu��o do �ngulo de teta(t)');
legend('beta = 0', 'beta = 1');

figure(2)
xlabel('Tempo [s]');
ylabel('d teta(t)');
title('Evolu��o da velocidade de teta(t)');
legend('beta = 0', 'beta = 1');

figure(3)
xlabel('teta(t)');
ylabel('d teta(t)');
title('Tra�ado de fase');
legend('beta = 0', 'beta = 1');

for r = 4:1:5
    figure(r)
    xlabel('teta(t)');
    ylabel('d teta(t)');
    if r == 4
        title('Tra�ado de fase com beta = 0');
    elseif r == 5
        title('Tra�ado de fase com beta = 1');
    end
    leg = legend('(-5, -5)', '(-3, -3)', '(-1, -1)', '(1, 1)', '(3, 3)','(5, 5)');
    title(leg, 'Pontos iniciais (teta,dteta)');
end
    
%%
% Para $\beta$ = 0, o sistema n�o tem amortecimento, logo como se pode confirmar
% pela Figura 1, o �ngulo $\theta$ oscila indefinidamente de forma 
% constante sem perturba��es, enquanto que para $\beta$ = 1 o sistema tem
% amortecimento total e portanto converge bastante r�pido para 0.
% Para o caso de $\beta$ = 0 podemos verificar que obtemos no nosso plano
% de fase uma din�mica do tipo "Centro", isto deve-se ao facto de os
% valores pr�prios serem imagin�rios puros. Intuitivamente, percebemos que
% o sistema sem qualquer tipo de atr�to, permanecer� em oscila��o
% indeterminadamente. J� no caso de $\beta$ = 1 os nossos valores pr�rpios
% s�o negativos, logo as nossas solu��es ir�o tender para a origem,
% orientando-se e acompanhando os vectores pr�prios existentes. Fazendo
% novamente uma an�lise intuitiva, visto que existe uma for�a de atr�to
% presente no sistema, a oscila��o ir� ser amortecida at� se anular no
% limite do tempo.


%% Exerc�cio 8

clear
close all
clc


%inicializa��o de variaveis
g = 9.8;
k = 3;
L = 0.5;
M = 0.15;
l = 0.4;
m = 0.2;
J = m*l^2+(M*L^2)/3;

beta = 1;

stop_time = 10;


%definir as matrizes de estado
A =  [          0                    1      ;
     ((-k+g*(m*l+(M*L/2)))/J)    (-beta/J) ];

B =  [        0          ;         1/J     ];

C =  [        1                     0       ;
              0                     1      ];

D =  [        0          ;          0      ];
  
T =  [        0                     0      ];

%calculo das condi�oes iniciais
[vectors, values] = eig(A);

for i = 1:2
    x0 = [vectors(1,i) vectors(2,i)];

    sim8 = sim('system6');
    figure(i)
    plot(sim8.y(:,1), sim8.y(:,2));
    title("Tra�ado de fase com condi��es iniciais: x0 = [" + vectors(1,i) + " ; " + vectors(2,i) + "]");
    xlabel('teta(t)');
    ylabel('d teta(t)');
end


%% Exerc�cio 9

clear
close all
clc

%inicializa��o de variaveis
g = 9.8;
k = 0.35;
L = 0.25;
M = 0.1;
beta = 0.001;

m_lento = 0;
l_lento = 0;
m_allegro = 0;
l_allegro = 0;


%definir BPMs que queremos (grupo 36 -> X = 3; Y = 6)
X = 3;
Y = 6;

BPM_lento = 50 + 2*(X+1);

BPM_allegro = 150 - 2*(Y+2);

%para fazermos o dimensionamento iremos correr o modelo para varios valores
%de m e l at� encontrarmos o q estiver mais proximo do pretendido
m_array = linspace(0, M, 200);
l_array = linspace(0.05, L, 200);

%calcula o BPM para cada valor de m e l
for i = 1:size(m_array,2)
    for j = 1:size(l_array,2)
        m = m_array(i);
        l = l_array(j);
        
        J = m*l^2+(M*L^2)/3;
        
        wn(i,j) = sqrt((k-g*(m*l+(M*L/2)))/J);
        
        epsilon(i,j) = (beta/(2*J*wn(i,j)));
        
        %frequencia do oscilador amortecido -> mais fiel q frequencia
        %natural
        wa(i,j) = wn(i,j)*sqrt(1-epsilon(i,j)^2);
        
        %para eliminar os valores complexos
        if ~isreal(wa(i,j))
            wa(i,j) = NaN;
        end
        
         BPM_matrix(i,j) = wa(i,j)*60/pi;          
   end
end

%encontra o primeiro valor de m que permite um BPM maior que o BPM_allegro e um BPM
%menor que o BPM_lento
for i = 1:size(m_array,2)
    m_dimensionamento_lento = 0;
    m_dimensionamento_allegro = 0;
    for j = 1:size(l_array,2)
        if BPM_matrix(i,j) >= BPM_allegro
            m_dimensionamento_allegro = 1;
        end
        if BPM_matrix(i,j) <= BPM_lento
            m_dimensionamento_lento = 1;
        end
    end    
    if (m_dimensionamento_lento + m_dimensionamento_allegro) == 2
        index_m_dimensionado = i;
        break;
    end
end
if m_dimensionamento_lento == 0
    index_m_dimensionado = size(m_array,2);
elseif m_dimensionamento_allegro == 0
    index_m_dimensionado = 1;
end

[index_l_lento] = findClosest(BPM_matrix, BPM_lento, index_m_dimensionado, size(l_array,2));
[index_l_allegro] = findClosest(BPM_matrix, BPM_allegro, index_m_dimensionado, size(l_array,2));

m_lento = m_array(index_m_dimensionado);
m_allegro = m_lento;
l_lento = l_array(index_l_lento);
l_allegro = l_array(index_l_allegro);


%%%%%%%%%%%%%%%
%plot do BPM
figure(1)
surfc(m_array, l_array, BPM_matrix);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%confirmar o dimensionamento acima

for i = 1:2
    if i == 1
        m = m_lento;
        l = l_lento;
    elseif i == 2
        m = m_allegro;
        l = l_allegro;
    end
    
    %recalcular matrizes e parametros necessarios para a simula�ao
    J = m*l^2+(M*L^2)/3;

    A =  [          0                    1      ;
         ((-k+g*(m*l+(M*L/2)))/J)    (-beta/J) ];

    B =  [        0          ;         1/J     ];

    C =  [        1                     0       ;
                  0                     1      ];

    D =  [        0          ;          0      ];

    T =  [        0                     0      ];

    x0 = [       pi/4                   0      ];  %condi�ao inicial

    stop_time = 60;

    sim9 = sim('system6');

    figure(i+1);
    plot(sim9.tout, sim9.y(:,1));
    
    wn = sqrt((k-g*(m*l+(M*L/2)))/J);       
    epsilon = (beta/(2*J*wn));
    
    %faz plot da envolvente teorica
    hold on
    plot(sim9.tout, (pi/4)*exp(-epsilon*wn*sim9.tout), 'r');
    plot(sim9.tout, -(pi/4)*exp(-epsilon*wn*sim9.tout), 'r');

    %estimar BPM 
    [peaks, locations] = findpeaks(sim9.y(:,1));
    
    T_estimada = (sim9.tout(locations(numel(locations))) - sim9.tout(locations(1)))/(numel(locations)-1);

    f_estimada = 1/T_estimada;
    
    BPM_estimado(i) = 2*f_estimada*60;
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%graficos

figure(1);
title('Varia��o de BPM para diferentes dimensionamentos de m e l');
xlabel('m');
ylabel('l');
zlabel('BPM');
shading interp;
colorbar;

figure(2);
xlabel('Tempo [s]');
ylabel('teta(t)');
title("Evolu��o do �ngulo de teta(t) (modelo linear) com: m = " + m_lento + " ; l = " + l_lento);
legend("BPM real =" + BPM_estimado(1) + "; BPM te�rico =" + BPM_lento, "Envolvente te�rica");
grid on;

figure(3);
xlabel('Tempo [s]');
ylabel('teta(t)');
title("Evolu��o do �ngulo de teta(t) (modelo linear) com: m = " + m_allegro + " ; l = " + l_allegro);
legend("BPM real =" + BPM_estimado(2) + "; BPM te�rico =" + BPM_allegro, "Envolvente te�rica");
grid on;

%% 
% Como podemos verificar pelas Figuras 1 e 2, a evolu��o de $\theta$(t) ao
% longo do tempo segue de forma bastante fiel a envolvente te�rica, como
% seria de esperar pois a estima��o do BPM est� bastante perto do desejado
% para ambas as cad�ncias, n�o se verificando nenhuma diferen�a de maior
% entre a estima��o do BPM lento para o allegro.


%% Fun��o findClosest

type('findClosest.m');


%% Exerc�cio 10

clear all
close all
clc

%inicializa��o de variaveis
g = 9.8;
L = 0.25;
M = 0.1;
k = 0.35;
beta = 0.001;
T = 0;

%condi�oes iniciais
teta0 = pi/4;
d_teta0 = 0;

%valores calculados em 9)
m_lento =  0.068844;
l_lento = 0.25;
m_allegro = m_lento;
l_allegro = 0.11633;

%definir BPMs que queremos (grupo 36 -> X = 3; Y = 6)
X = 3;
Y = 6;

BPM_lento = 50 + 2*(X+1);

BPM_allegro = 150 - 2*(Y+2);

%atribuir a m e l os valores dimensionados em 9)
for i = 1:2   
    if i == 1
        m = m_lento;
        l = l_lento;
    elseif i == 2
        m = m_allegro;
        l = l_allegro;
    end

    J = (3*m*l^2+M*L^2)/3;

    stop_time = 60;

    sim10 = sim('system10');

    figure(i)
    plot(sim10.tout, sim10.teta);
    hold on

    %estimar BPM
    [peaks, locations] = findpeaks(sim10.teta);

    T_estimada = (sim10.tout(locations(numel(locations))) - sim10.tout(locations(1)))/(numel(locations)-1);

    f_estimada = 1/T_estimada;

    BPM_estimado(i) = 2*f_estimada*60;
    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %queremos otimizar a posi�ao da massa m (l) para aproximar BPM 
    %do valor desejado
    l_array = linspace(l-0.001,l+0.001,20);
    
    for j = 1:size(l_array,2)
        l = l_array(j);
        J = (3*m*l^2+M*L^2)/3;
        stop_time = 60;
        sim10_testa_opti = sim('system10');
        [peaks_testa_opti, locations_testa_opti] = findpeaks(sim10_testa_opti.teta);
        T_testa_opti = (sim10_testa_opti.tout(locations_testa_opti(numel(locations_testa_opti))) - sim10_testa_opti.tout(locations_testa_opti(1)))/(numel(locations_testa_opti)-1);
        f_testa_opti = 1/T_testa_opti;
        BPM_testa_opti(1,j) = 2*f_testa_opti*60;
    end
    
    %encontra o valor de BPM mais proximo do desejado
    if i == 1
        l_index = findClosest(BPM_testa_opti, BPM_lento, 1, size(l_array,2));
    elseif i == 2
        l_index = findClosest(BPM_testa_opti, BPM_allegro, 1, size(l_array,2));
    end
    
    l = l_array(l_index);

    %volta a correr a simula�ao
    J = (3*m*l^2+M*L^2)/3;
    stop_time = 60;

    sim10_otimizado = sim('system10');

    figure(i)
    plot(sim10_otimizado.tout, sim10_otimizado.teta);

    %estimar BPM
    [peaks_otimizado, locations_otimizado] = findpeaks(sim10_otimizado.teta);

    T_estimada_otimizado = (sim10_otimizado.tout(locations_otimizado(numel(locations_otimizado))) - sim10_otimizado.tout(locations_otimizado(1)))/(numel(locations_otimizado)-1);

    f_estimada_otimizado = 1/T_estimada_otimizado;

    BPM_estimado_otimizado(i) = 2*f_estimada_otimizado*60;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%graficos
figure(1);
xlabel('Tempo [s]');
ylabel('teta(t)');
title("Evolu��o do �ngulo de teta(t) (modelo n�o linear) ");
leg1 = legend("BPM real =" + BPM_estimado(1) , "BPM otimizado =" + BPM_estimado_otimizado(1));
title(leg1, "BPM te�rico =" + BPM_lento);
grid on;

figure(2);
xlabel('Tempo [s]');
ylabel('teta(t)');
title("Evolu��o do �ngulo de teta(t) (modelo n�o linear)");
leg2 = legend("BPM real =" + BPM_estimado(2) , "BPM otimizado =" + BPM_estimado_otimizado(2));
title(leg2, "BPM te�rico =" + BPM_allegro);
grid on;

%%
% Como podemos verificar, ao usarmos o dimensionamento calculado em 9) para
% um modelo n�o linear do sistema verifica-se que a estima��o de BPM
% afastou-se um pouco do valor te�rico. De maneira a optimizarmos esta
% estima��o, assumindo que a massa m � a mesma, iremos voltar a simular o 
% modelo do sistema para um intervalo de valores na vizinhan�a do
% comprimento l calculado anteriormente e voltamos a proceder � mesma
% t�ctica do 9, isto �, calculamos o BPM real para todos os dimensionametos
% e calculamos qual o que tem menor dist�ncia do valor te�rico, sendo esse
% o dimensionamento escolhido para a optimiza��o.


%% Modelo SIMULINK 10

system10


%% Exerc�cio 11

clear
close all
clc

%inicializa��o de variaveis
g = 9.8;
k = 0.35;
L = 0.25;
M = 0.1;
beta = 0.001;

%condi�oes iniciais
teta0 = pi/4;
d_teta0 = 0;

%valores calculados em 9)
m_lento =  0.068844;
l_lento = 0.25;
m_allegro = m_lento;
l_allegro = 0.11633;

%definir BPMs que queremos (grupo 36 -> X = 3; Y = 6)
X = 3;
Y = 6;

BPM_lento = 50 + 2*(X+1);

BPM_allegro = 150 - 2*(Y+2);

T = 0.5;

for i = 1:2   
    if i == 1
        m = m_lento;
        l = l_lento;
    elseif i == 2
        m = m_allegro;
        l = l_allegro;
    end

    J = (3*m*l^2+M*L^2)/3;

    stop_time = 60;

    sim11 = sim('system11');
    
    figure(i)
    plot(sim11.tout, sim11.teta);
    hold on
    plot(sim11.tout, sim11.binario);
    
    %estimar BPM 
    [peaks, locations] = findpeaks(sim11.teta);
    
    T_estimada = (sim11.tout(locations(numel(locations))) - sim11.tout(locations(1)))/(numel(locations)-1);

    f_estimada = 1/T_estimada;
    
    BPM_estimado(i) = 2*f_estimada*60;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%graficos
figure(1);
xlabel('Tempo [s]');
ylabel('teta(t)');
leg = legend("Teta(t)", "Bin�rio");
title(leg, "BPM real =" + BPM_estimado(1) + "; BPM te�rico =" + BPM_lento);
title("Evolu��o do �ngulo teta(t) e do bin�rio externo aplicado");
grid on;

figure(2);
xlabel('Tempo [s]');
ylabel('teta(t)');
leg = legend("Teta(t)", "Bin�rio");
title(leg, "BPM real =" + BPM_estimado(2) + "; BPM te�rico =" + BPM_allegro);
title("Evolu��o do �ngulo teta(t) e do bin�rio externo aplicado");
grid on;

%%
% Como se pode verificar, a aplica��o de um bin�rio externo ao sistema
% afecta bastante a frequ�ncia de oscila��o pretendida, mas no entanto
% verifica-se que consegue evitar o decaimento do sistema oscilat�rio e com
% isso obt�m-se uma oscila��o constante ao longo do tempo. Isto deve-se ao
% facto de estarmos a contrariar as for�as de atrito e grav�ticas do
% sistema que t�m contribui��es para o amortecimento da oscila��o. Podemos
% concluir tamb�m que o nosso bin�rio externo, caso seja maior do que o
% suposto, acaba por acelerar a oscila��o do metr�nomo, causando um aumento 
% significativo na frequ�ncia de oscila��o, desfasando-se do pretendido.

%% Modelo SIMULINK 11

system11

%% Exerc�cio 12

clear
close all
clc

%inicializa��o de variaveis
g = 9.8;
k = 0.35;
L = 0.25;
M = 0.1;
beta = 0.001;

%condi�oes iniciais
teta0 = pi/4;
d_teta0 = 0;

%valores calculados em 9)
m_lento =  0.068844;
l_lento = 0.25;
m_allegro = m_lento;
l_allegro = 0.11633;

%definir BPMs que queremos (grupo 36 -> X = 3; Y = 6)
X = 3;
Y = 6;

BPM_lento = 50 + 2*(X+1);

BPM_allegro = 150 - 2*(Y+2);

T = 0.5;

for i = 1:2   
    if i == 1
        m = m_lento;
        l = l_lento;
    elseif i == 2
        m = m_allegro;
        l = l_allegro;
    end

    J = m*l^2+(M*L^2)/3;
    wn = sqrt((k-g*(m*l+(M*L/2)))/J);
    epsilon = (beta/(2*J*wn));
    G = 1/((wn^2)*J);
    
    Gs = G*tf(wn^2, [ 1 , (2*epsilon*wn) , wn^2 ] );
    
    bode(Gs);
    hold on
    
end

figure(1)
grid on
legend("BPM lento", "BPM allegro");

%%
% Como se pode verificar, para o BPM mais elevado o diagrama de bode est�
% deslocado para a direita em rela��o ao BPM mais pequeno, o que �
% expect�vel pois o BPM varia linearmente com a frequ�ncia. Tamb�m se
% verifica que o ganho para o BPM allegro � inferior.

%% Exerc�cio 13

clear
close all
clc

%inicializa��o de variaveis
g = 9.8;
k = 0.35;
L = 0.25;
M = 0.1;
beta = 0.001;
m = 0.07;
l = 0.12;

J = m*l^2+(M*L^2)/3;
wn_teorico = sqrt((k-g*(m*l+(M*L/2)))/J);


T = 0.5;

%definir as matrizes de estado
A =  [          0                    1      ;
     ((-k+g*(m*l+(M*L/2)))/J)    (-beta/J) ];

B =  [        0          ;         1/J     ];

C =  [        1                     0       ;
              0                     1      ];

D =  [        0          ;          0      ];

x0 = [        0          ;        pi/4     ];  % condi�oes iniciais

stop_time = 60;

%cria um array de frequencias
f_array = linspace(0, 10, 500);

%cria um array com amplitudes
a_array = zeros(1,size(f_array,2));

max_a = 0;
max_a_index = 1;

for i = 1:size(f_array,2)  
    f = f_array(i);

    sim13 = sim('system13');
    
    a_array(i) = mean(sim13.y(:,1));

    if(a_array(i) > max_a)
        max_a = a_array(i);
        max_a_index = i;
    end
        
end

wn = 2*pi*f_array(max_a_index);

m_estimado = (k-wn^2*(M*(L^2)/3)-(g*M*L/2))/((wn^2)*(l^2)+g*l);

fprintf("\n m real = ");
disp(m);
fprintf("\n m estimado = ");
disp(m_estimado);

figure(1)
plot(f_array, a_array);

grid on;
xlabel('Frequ�ncia [Hz]');
ylabel('Amplitude [V]');

%%
% Como se pode confirmar pela quest�o 11, apesar de o mecanismo de
% relojoaria que aplica um bin�rio externo resolver o problema do
% decaimento para zero do sistema, obtendo assim ondula��es constantes, o
% que � o que se pretende quando se usa um metr�nomo. Apesar disso, podemos
% verificar que a frequ�ncia de oscila��o est� bastante longe da
% pretendida, o que torna este mecanismo in�til para uma aplica��o pr�ctica
% pois se o metr�nomo n�o consegue oscilar � frequ�ncia pretendida com um
% erro m�nimo nao tem nenhum uso pr�ctico. Para contrariar isto desenvolve-se um sistema mec�nico para criar uma
% "balan�a". Infelizmento a nossa estima��o n�o est� muito precisa pois temos um erro superior a 100%, o que obviamente est� errado, mas esta solu��o poderia
% possivelmente contrariar o erro causado pela aplica��o do bin�rio
% externo.

%% Modelo SIMULINK 13

system13