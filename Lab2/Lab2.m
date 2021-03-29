% LAB 2 - MSIM
% Autor: Bernardo Rocha & Vasco Ara�jo
% N�mero de Aluno: 89867 & 90817
% Abril 2020; �ltima Revis�o: 29/04/2020                    

%% Exerc�cio 2

clear
close all
clc

t = linspace(-1, 1, 100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.1;

pb = generatePulse(t, beta);

plot(t, pb)
hold on
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.5;

pb = generatePulse(t, beta);

plot(t, pb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.9;

pb = generatePulse(t, beta);

plot(t, pb)


ylim([0 1.3])
xlim([-1 1])
xlabel('Tempo [s]')
ylabel('Pb(t)')
legend('Beta = 0.1', 'Beta = 0.5', 'Beta = 0.9')
title('Gera��o do impulso prot�tipo p_\beta(t)')

%% Fun��o generatePulse

type('generatePulse.m');


%% Exerc�cio 3

clear
close all
clc

%declara��o de vari�veis
T = 5;
alpha = 1;
U1 = 1;
U2 = 2;
n1 = 100;
n2 = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.1;
[ut, t] = generateU(T, alpha, beta, U1, U2, n1, n2);
plot(t, ut);
hold on
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.5;
[ut, t] = generateU(T, alpha, beta, U1, U2, n1, n2);
plot(t, ut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.9;
[ut, t] = generateU(T, alpha, beta, U1, U2, n1, n2);
plot(t, ut);


ylim([-U1-1 U2+1])
xlim([0 T])
xlabel('Tempo [s]')
ylabel('u(t)')
legend('Beta = 0.1', 'Beta = 0.5', 'Beta = 0.9')
title('Gera��o de u(t) a partir do impulso prot�tipo p_\beta(t)')

%% Fun��o generateU

type('generateU.m');

%% Exerc�cio 7

clear
close all
clc

%inicializar vari�veis
b = 0;
beta = 1;
n1 = 100;
n2 = 100;


for alpha = 0.5:0.5:1.5

    if alpha <= 1
        T = sqrt(2*(1+beta)*(1+alpha)/alpha);
    else
        T = sqrt(2*(1+beta)*(1+alpha));
    end
    T1 = T/(1+alpha);
    U1 = 2*(1+beta)/((T1^2)*(1+alpha));
    U2 = U1/alpha;
    [ut, t] = generateU(T,alpha,beta,U1,U2,n1,n2);
    figure(1);
    plot(t, ut);
    hold on;
    grid on;
    
    u_entrada.time = t';
    u_entrada.signals.values = ut';
    stop_time = T;
    sim7 = sim('system7');
    
    figure(2);
    plot(sim7.tout, sim7.y_saida);
    hold on;
    grid on;
    
    figure(3);
    plot(sim7.tout, sim7.d_y_saida);
    hold on;
    grid on;
    
    %sistema com b = 0.025
    b = 0.025;
    sim7_2 = sim('system7');
    
    figure(4);
    plot(sim7_2.tout, sim7_2.y_saida);
    hold on;
    grid on;
    
    figure(5);
    plot(sim7_2.tout, sim7_2.d_y_saida);
    hold on;
    grid on;
    
    %repor o valor de b
    b = 0;
    
    
end


figure(1)
xlabel('Tempo [s]')
ylabel('u(t)')
legend('alpha = 0.5', 'alpha = 1', 'alpha = 1.5', 'location', 'southeast');

title('Sinal de entrada do sistema - u(t)');

figure(2)
xlabel('Tempo [s]')
ylabel('y(t)')
legend('alpha = 0.5', 'alpha = 1', 'alpha = 1.5');
title('Evolu��o da posi��o do sistema - y(t) com b = 0')

figure(3)
xlabel('Tempo [s]')
ylabel('dy(t)')
legend('alpha = 0.5', 'alpha = 1', 'alpha = 1.5', 'location', 'southeast');
title('Evolu��o da velocidade do sistema - dy(t) com b = 0')

figure(4)
xlabel('Tempo [s]')
ylabel('y(t)')
legend('alpha = 0.5', 'alpha = 1', 'alpha = 1.5');
title('Evolu��o da posi��o do sistema - y(t) com b = 0.025');

figure(5)
xlabel('Tempo [s]')
ylabel('dy(t)')
legend('alpha = 0.5', 'alpha = 1', 'alpha = 1.5', 'location', 'southeast');
title('Evolu��o da velocidade sistema - dy(t) com b = 0.025');

%%
% Na al�nea 6 foi calculado que o valor de $\alpha$ que minimizava o T era
% $\alpha = 1$. Como se pode observar nas Figuras 2 e 4, o c�lculo te�rico
% confirma-se pois y(t) atinge o valor desejado de 0 mais rapidamente para
% estas condi��es.
% Com o sistema perturbado, verifica-se que este sistema n�o atinge o valor
% pretendido.
%% Sistema de SIMULINK 7

system7


%% Exerc�cio 8 

clear
close all
clc


test.y = linspace(-10,10,100);
test.y_dot = linspace(-10,10,100);

for n = 1:length(test.y)
    for m = 1:length(test.y_dot)
        u(n,m) = blockFunc8(test.y(n), test.y_dot(m));
    end
end

% curva de comuta�ao 2d
figure(1);
surfc(test.y, test.y_dot, u);
view([0 90])
xlabel('y');
ylabel('dy');
title(' u(y, dy) - Curva de comuta��o');
colorbar
shading flat;

%curva de comuta��o 3d
figure(2);
surfc(test.y, test.y_dot, u);
view([46 51])
xlabel('y');
ylabel('dy');
zlabel('u');
title('Gera��o de u(y, dy)');
shading interp;
colorbar;

%% Fun��o blockFun8

type('blockFunc8.m');

%%
% Nesta estrutura de controlo os nossos impulsos de entrada s�o definidos em
% fun��o da diferen�a entre a posi��o do bra�o e a refer�ncia dos sistema.
% Com isto conseguimos ter um maior controlo da posi��o do bra�o para que
% este se mova para a posi��o desejada.
% Temos que a nossa tens�o de entrada � influenciada tanto pela posi��o do 
% nosso bra�o como pela derivada da posi��o do mesmo.
% Neste sistema de controlo, podemos resumir a defini��o do sinal do impulso
% a duas perguntas: A posi��o est� acima ou abaixo da refer�ncia? 
% A derivada da posi��o � positiva ou negativa?
% Com a resposta a estas perguntas, conseguimos definir o sinal de u de 
% forma a que o bra�o se posicione na refer�ncia.
% Podemos usar um exemplo: Imaginando que o nosso y est� abaixo da nossa 
% refer�ncia 0, logo temos y < 0. Ao fazer a diferen�a deste n�mero 
% negativo com 0, passamos a ter um n�mero positivo � entrada de Subsys. 
% Logicamente, � sa�da de Subsys iremos ter um valor positivo novamente 
% (verificado atrav�s da substitui��o de um x positivo em f(x)).
% At� agora apenas utiliz�mos a resposta � primeira pergunta. 
% Questionamos agora: Ser� que a derivada da posi��o � positiva ou negativa? 
% Imaginemos que estamos perante o caso de uma derivada negativa. Isto quer 
% dizer que para al�m da posi��o do nosso bra�o se situar abaixo do desejado, 
% este tende a se afastar ainda mais, logo, o que n�s queremos � que este fa�a 
% precisamente o contr�rio e fique com uma derivada positiva, para que o valor 
% de y cres�a em dire��o � refer�ncia.
% Como para este exemplo temos um valor � sa�da de Subsys positivo e a nossa 
% derivada tem sinal negativo, a diferen�a entre estes dois ir� dar um valor 
% positivo � sa�da do somador. O mesmo acontecer� � sa�da do bloco Sign. 
% Com isto, alter�mos o declive da varia��o da posi��o y. Ap�s integrar este 
% valor positivo no 1� bloco integrador, constatamos que iremos ter uma 
% derivada positiva, que era aquilo que n�s quer�amos: ter o nosso y, que 
% est� abaixo da refer�ncia, a crescer e a aproximar-se de 0.

%% Exerc�cio 9

clear
close all
clc


%inicializa��o de vari�veis
K = 10000;
y0 = 1;
d_y0 = 0;

sim9 = sim('system9');

figure(1)
plot(sim9.tout, sim9.ut);
hold on
grid on

figure(2)
plot(sim9.tout, sim9.y);
hold on
grid on

figure(3)
plot(sim9.tout, sim9.d_y);
hold on
grid on

%fazer tra�ado de fase da pergunta 5
figure(4)
for i = -5:2:5
        y0 = i;
        d_y0 = i;
        sim9_i = sim('system9');
        plot(sim9_i.y, sim9_i.d_y);
        hold on;
        
end


figure(1)
xlabel('Tempo [s]')
ylabel('u(t)')
title('Sinal de entrada do sistema - u(t)');

figure(2)
xlabel('Tempo [s]')
ylabel('y(t)')
title('Evolu��o da posi��o do sistema - y(t)')

figure(3)
xlabel('Tempo [s]')
ylabel('dy(t)')
title('Evolu��o da velocidade do sistema - dy(t)')

figure(4)
xlabel('y(t)')
ylabel('dy(t)')
leg = legend('(-5, -5)', '(-3, -3)', '(-1, -1)', '(1, 1)', '(3, 3)','(5, 5)');
title(leg, 'Pontos iniciais (y,dy)');
title('Tra�ado de ')


%%
% Na Figura 1 observamos um ret�ngulo a partir do
% momento em que o estado do sistema converge para a origem. Este fen�meno,
% chamado $\textit{chattering}$, ocorre quando se tem oscila��es com frequ�ncia e amplitude
% finitas, causadas quando se tem din�micas do sistema a variar muito
% rapidamente, t�o r�pido que o pr�prio sistema n�o consegue acompanhar. 
% Este fen�meno � prejudicial pois pode causar: pouca precis�o
% no sistema, desgaste das partes mec�nicas em movimento e grandes
% perdas de calor no circuito.

%% Sistema de SIMULIK 9
 
system9


%Dentro do bloco Subsys encontra-se um bloco MATLAB Function com a seguinte
%fun��o:

% function y  = fcn(x)
% 
%     y = sign(x) * sqrt(2*abs(x));
% 
% end

%% Exerc�cio 10

clear
close all
clc

y1 = 4;

test.y = linspace(-10,10,100);
test.y_dot = linspace(-10,10,100);



for n = 1:length(test.y)
    for m = 1:length(test.y_dot)
        u(n,m) = blockFunc10(test.y(n), test.y_dot(m), y1);
    end
end


% curva de comuta�ao 2d
figure(1);
surfc(test.y, test.y_dot, u);
view([0 90])
grid on;
xlabel('y');
ylabel('dy');
title(' u(y, dy) - Curva de comuta��o com y_l = 4');
shading flat;
colorbar;

%curva de comuta�ao 3d
figure(2);
surfc(test.y, test.y_dot, u);
view([46 51])
xlabel('y');
ylabel('dy');
zlabel('u');
title('Gera��o de u(y, dy) com y_l = 4');
shading interp;
colorbar;

%% Fun��o blockFun10

type('blockFunc10.m');

%%
% Esta arquitectura de controlo garante que quando a posi��o do bra�o se 
% encontra perto da refer�ncia a sa�da do bloco "Subsys" pode assumir
% valores que n�o $\pm$ 1, o que faz com que deixe de haver efeito de
% $\textit{chattering}$,como podemos observar na Figura 1, em que entre 
% $\pm$ y_l a curva de comuta��o � linear.

%% Exerc�cio 11

clear
close all
clc

%inicializa��o de vari�veis
b = 0;

y1 = 1;
k1 = 1/y1;
k2 = sqrt(2*k1);
y0 = 1;
d_y0 = 0;


sim11 = sim('system11');

figure(1)
plot(sim11.tout, sim11.ut);
hold on
grid on

figure(2)
plot(sim11.tout, sim11.y);
hold on
grid on

figure(3)
plot(sim11.tout, sim11.d_y);
hold on
grid on

%fazer tra�ado de fase da pergunta 5
figure(4)
for i = -5:2:5
        y0 = i;
        d_y0 = i;
        sim11_i = sim('system11');
        plot(sim11_i.y, sim11_i.d_y);
        hold on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%resposta do sistema linear s� com a parte de cima do ramo da fun��o
sim11_1 = sim('system11_linear');

figure(1)
plot(sim11_1.tout, sim11_1.ut);
hold on
grid on

figure(2)
plot(sim11_1.tout, sim11_1.y);
hold on
grid on

figure(3)
plot(sim11_1.tout, sim11_1.d_y);
hold on
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%resposta do sistema em malha fechada sem chattering para compara��o
K = 10;      %reduzimos o ganho face ao utilizado na quest�o 9 para diminuir um pouco 
%as oscila��es de chattering para se poder comparar com o sistema sem chattering de forma mais fidedigna
y0 = 1;
d_y0 = 0;

sim9 = sim('system9');

figure(1)
plot(sim9.tout, sim9.ut);
hold on
grid on

figure(2)
plot(sim9.tout, sim9.y);
hold on
grid on

figure(3)
plot(sim9.tout, sim9.d_y);
hold on
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
xlabel('Tempo [s]')
ylabel('u(t)')
title('Sinal de entrada do sistema - u(t) com b = 0');
legend('Sistema malha fechada sem chattering', 'Sistema malha fechada sem chattering linear', 'Sistema malha fechada com chattering', 'location', 'southeast');

figure(2)
xlabel('Tempo [s]')
ylabel('y(t)')
title('Evolu��o da posi��o do sistema - y(t) com b = 0')
legend('Sistema malha fechada sem chattering', 'Sistema malha fechada sem chattering linear', 'Sistema malha fechada com chattering');

figure(3)
xlabel('Tempo [s]')
ylabel('dy(t)')
title('Evolu��o da velocidade do sistema - dy(t) com b = 0')
legend('Sistema malha fechada sem chattering', 'Sistema malha fechada sem chattering linear', 'Sistema malha fechada com chattering');

figure(4)
xlabel('y(t)')
ylabel('dy(t)')
leg = legend('(-5, -5)', '(-3, -3)', '(-1, -1)', '(1, 1)', '(3, 3)','(5, 5)');
title(leg, 'Pontos iniciais (y,dy)');
title('Tra�ado de fase do sistema malha fechada sem chattering')


%%
% Olhando para o sistema de controlo sem $\textit{chattering}$,
% verificamos que o sinal de controlo tem uma curva bastante mais suave e
% j� nao tem oscila��es, como seria de esperar. No entanto esta arquitectura
% tem um trade-off, que �  facto de o tempo que demora para a posi��o do
% bra�o tender para o equil�brio ser consideravelmente maior, pois no nosso
% exemplo o sistema com $\textit{chattering}$ estabilizava para t = 2.0211
% s, enquanto que o com $\textit{chattering}$ apenas estabiliza nos 7.5612
% s.
% Podemos tamb�m observar a resposta do sistema com um controlo linear,
% alterando a fun��o $\textit{f(x)}$ para apenas o ramo superior e retirando
% o bloco de satura��o e verificamos que, apesar de o sistema estabilizar �
% mesma, apresenta uma maior sobre-eleva��o no sinal de controlo e uma maior 
% oscila��o e tempo de estabiliza��o na posi��o e na velocidade.



%% Sistema de SIMULINK 11

system11


%Dentro do bloco Subsys encontra-se um bloco MATLAB Function com a seguinte
%fun��o:

% function y = fcn(x)
% y1 = 1;
% k1 = 1/y1;
% k2 = sqrt(2*k1);
%     
%     if abs(x) <= y1
%         y = (k1/k2)*x;
%     elseif abs(x) > y1    
%         y = sign(x) * (sqrt(2*abs(x))-(1/k2));
%     else
%         y = 0;
%     end
% 
% end

%% Exerc�cio 12

clear
close all
clc


%inicializar vari�veis
b = 0.025;

y1 = 1;
k1 = 1/y1;
k2 = sqrt(2*k1);
y0 = 1;
d_y0 = 0;

sim12 = sim('system11');

figure(1)
plot(sim12.tout, sim12.ut);
hold on
grid on

figure(2)
plot(sim12.tout, sim12.y);
hold on
grid on

figure(3)
plot(sim12.tout, sim12.d_y);
hold on
grid on

%fazer tra�ado de fase da pergunta 5
figure(4)
for i = -5:2:5
        y0 = i;
        d_y0 = i;
        sim12_i = sim('system11');
        plot(sim12_i.y, sim12_i.d_y);
        hold on;
end


figure(1)
xlabel('Tempo [s]')
ylabel('u(t)')
title('Sinal de entrada do sistema - u(t) com b = 0.025');

figure(2)
xlabel('Tempo [s]')
ylabel('y(t)')
title('Evolu��o da posi��o do sistema - y(t) com b = 0.025')

figure(3)
xlabel('Tempo [s]')
ylabel('dy(t)')
title('Evolu��o da velocidade do sistema - dy(t) com b = 0.025')

figure(4)
xlabel('y(t)')
ylabel('dy(t)')
legend('(-5, -5)', '(-3, -3)', '(-1, -1)', '(1, 1)', '(3, 3)','(5, 5)');
title('Tra�ado de fase')

%%
% Comparando com as simula��es efectuadas em 7) com b = 0.025, conclu�mos
% que o sistema em cadeia fechada � bastante mais robusto que o de cadeia
% aberta, pois em cadeia aberta a posi��o e a velocidade do sistema com
% atrito nunca chegavam ao ponto de equil�brio pretendido, enquanto que
% em cadeia fechada o sistema estabilizava para o equil�brio.

%% Exerc�cio 13

clear
close all
clc

%inicializar vari�veis
b = 0;

y1 = 0.1;
k1 = 1/y1;
k2 = sqrt(2*k1);
y0 = 0;
d_y0 = 0;

%1� referencia
sim13 = sim('system13');

figure(1)
plot(sim13.tout, sim13.y);
hold on
grid on
plot(sim13.tout, sim13.reference);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
xlabel('Tempo [s]')
ylabel('y(t)')
title('Evolu��o da posi��o do sistema')
legend('y(t)', 'Sinal de refer�ncia');
ylim([-1 5]);


%% Sistema de SIMULINK 13

system13

%Dentro do bloco Subsys encontra-se um bloco MATLAB Function com a mesma
%fun��o do sistema de Simulink system11

%%
% Como podemos observar, o sistema de controlo segue as refer�ncias n�o
% nulas com bastante fidelidade, apesar de ter um pouco de atraso a
% acompanhar a refer�ncia quando esta muda de amplitude bruscamente.
