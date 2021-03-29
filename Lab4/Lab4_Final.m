%%
% LAB 4 - MSIM
% Autor: Bernardo Rocha & Vasco Ara�jo
% N�mero de Aluno: 89867 & 90817
% Junho 2020; �ltima Revis�o: 03/06/2020   

%% Exerc�cio 2.a

close all
clear all
clc

MarkovChain = load('MarkovChain.mat');

P = MarkovChain.P;

[vectors, values] = eig(P');

old_diference = 10000;
val_u = 1;
%encontra o valor proprio mais perto de 1
for i = 1:size(values,1)
    diference = abs(values(i,i) - 1);
    if diference < old_diference
        old_diference = diference;
        value_index = i;
    end
end

%vector a normalizar
vec_n = vectors(:,value_index);

normalize = vec_n/sum(vec_n);

%se esta tudo bem soma das probabilidades tem que dar 1
prob_sum = sum(normalize)

%faz gr�fico de barras
figure(1)
bar(normalize);
title('Distribui��o de equil�brio da cadeia de Markov');
xlabel('N�mero do estado da cadeia de Markov');
ylabel('Probabilidade do token estar na posse do agente');
box on;
grid on;

%encontra estado mais provavel
high_prob = 0;
low_prob = 1;
%encontra o 1� estado mais provavel e o 1� estado menos provavel
for i = 1:size(normalize,1)
    if normalize(i) > high_prob
        high_prob = normalize(i);
    end
    if normalize(i) <= low_prob
        low_prob = normalize(i);
    end
end

%volta a correr para ver se ha mais estados com mesma probabilidade que o
%mais provavel e o menos provavel
h = 1;
l = 1;
for i = 1:size(normalize,1)
    if normalize(i) == high_prob
        highest_prob_index(h) = i;
        highest_prob(h) = normalize(i)
        highest_prob_anchor = i
        h = h + 1;
    end
    if normalize(i) == low_prob
        lowest_prob_index(l) = i;
        lowest_prob(l) = normalize(i)
        lowest_prob_anchor = i
    end
end



%%
% Olhando para o gr�fico de barras poderia pensar-se que os estados mais
% prov�veis eram os estados 7 e 19, no entanto correndo o c�digo para
% encontrar os maiores estados apenas no d� o estado 7. Isto � porque,
% analisando os valores de normalize verifica-se que o estado 7 tem uma
% probabilidade de 0.096493187402279 e o estado 19 uma probabilidade de
% 0.096493187402278, por isso o c�digo apenas retorna 7 como o estado
% com maior probabilidade. O mesmo acontece com os estados 8 e 17 para o
% c�lculo do estado com menor probabilidade, visto que o estado 17 tem uma 
% probabilidade um pouco menor. Assim sendo, se quiseremos ser rigorosos o
% estado com maior probabilidade � o estado 7 e o estado com menor 
% probabilidade � o estado 8. 



        

%% Exerc�cio 2.b

close all
clear all
clc

MarkovChain = load('MarkovChain.mat');

%par�metros da simula��o
N = size(MarkovChain.nodePos,1);        %numero de ancoras
n = 2;                                  %numero de dimensoes
sidelength = 100;  

M = 10000;                               %numero de observa��es

%posi�oes das �ncoras
a = [MarkovChain.nodePos(:,2) , MarkovChain.nodePos(:,3)]';
%posi��es das sources
x = MarkovChain.sourcePos';


D = squareform(pdist([x zeros(size(x)) a]'));
d = D(1,3:end);                     % Source-anchor distances
an = D(2,3:end);                    % Anchor norms
    
P0 = 100;                    
mu = 0;
sigma = 1e-2;                       

%inicializa matrizes
A = zeros(M,4);
b = zeros(M,1);

%ruido ni
ni_space = [-0.1:.0002:0.1];
ni = normpdf(ni_space, mu, sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c�digo da 2.a
P = MarkovChain.P;
[vectors, values] = eig(P');
old_diference = 10000;
val_u = 1;
%encontra o valor proprio mais perto de 1
for i = 1:size(values,1)
    diference = abs(values(i,i) - 1);
    if diference < old_diference
        old_diference = diference;
        value_index = i;
    end
end
%vector a normalizar
vec_n = vectors(:,value_index);
normalize = vec_n/sum(vec_n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%numero de vezes que cada �ncora tem o token
number_meas = normalize.*M; 
%aproxima��o para o n�mero inteiro mais pr�ximo
number_meas = round(number_meas);

%calcula as matrizes A e b
j = 1;
for i = 1:N
    meas = number_meas(i);
    if meas ~= 0
        P(j) = (P0./(d(i).^2))*exp(ni(j));
        A(j,1) = -2*P(j)*a(1,i);
        A(j,2) = -2*P(j)*a(2,i);
        A(j,3) = -1;
        A(j,4) = P(j);
        b(j) = -P(j)*(an(i)^2)';
        j = j + 1;
        meas = meas-1;
    end
end

z = A\b;
fprintf('Erro com %d mudan�as de agente', M);
x_estimated = z(1:n);
error1 = norm(x-x_estimated)

% RLS formulation (one-shot)
% RlsPar = struct('lam',1);
% [e,w,RlsPar] = qrrls(A,b,RlsPar);
%error1 = norm(z-w)
%x_estimated = [w(1) w(2)]';



figure(1)
plot(a'*[1; 1i],'o'); hold all
plot(x'*[1; 1i],'x');plot(x_estimated'*[1; 1i],'s'); hold off
axis(sidelength*[0 1 0 1]); axis('square')
title('Estima��o da posi��o da source');
legend('�ncoras', 'Posi��o real da source', 'Posi��o estimada da source com 10000 mudan�as');
xlabel('x');
ylabel('y');

figure(2)
hold all
plot(x'*[1; 1i],'x'); plot(x_estimated'*[1; 1i],'s'); 
axis([84.9999 85.0001 29.9999 30.0001]);
title('Estima��o da posi��o da source ampliada');
legend('Posi��o real da source', 'Posi��o estimada da source com 10000 mudan�as', 'Location','southwest');
xlabel('x');
ylabel('y');


%%
% Decidiu-se usar a primeira formula��o apresentada no ficheiro 'rssiloc.m' por
% ser a que apresentava um erro menor. Usando o algoritmo de m�nimos quedrados
% com base na solu��o matricial consegue-se um
% resultado bastante fiel para a fonte estimada, dado um n�mero elevado de
% transi��es. O erro nunca poderia ser zero devido ao ru�do gaussiano na
% medi��o da pot�ncia da fonte, por isso consideramos o nosso erro de
% 1.0658e-14 bastante bom. 




%% Exerc�cio 2.c

clear all
close all
clc



MarkovChain = load('MarkovChain.mat');
P = MarkovChain.P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.a (para usar mais a frente)
[vectors, values] = eig(P');
old_diference = 10000;
val_u = 1;
%encontra o valor proprio mais perto de 1
for i = 1:size(values,1)
    diference = abs(values(i,i) - 1);
    if diference < old_diference
        old_diference = diference;
        value_index = i;
    end
end
%vector a normalizar
vec_n = vectors(:,value_index);
normalize = vec_n/sum(vec_n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 20;
time = 500;

pi_array = zeros(N,time);

%inicializar a matriz pi_array_0
%cada coluna corresponde ao vetor inicial para cada simula��o
pi_array_0 = zeros(N,4);

%1� condi��o inicial -> come�a no estado 7
pi_array_0(7,1) = 1;
%2� condi��o inicial -> come�a no estado 8
pi_array_0(8,2) = 1;
%3� condi��o -> todos os estados t�m a mesma probabilidade
for i = 1:20
    pi_array_0(i,3) = 1/20;
end

for i = 1:3
    eq_time(i) = 0;
end


for j = 1:3
    pi_array(:,1) = pi_array_0(:,j);
    for i = 2:time
        pi_array(:,i) = pi_array(:,1)'*P^(i);
    end
    
    %para comparar com os resultados da decomposi��o em valores e vetores
    %proprios (2.a) vamos calcular o tempo que demora a atingir o
    %equil�brio
    for i = 1:time
        if pi_array(:,i) - normalize < 0.000001
            eq_time(j) = i;
            break;
        end
    end
    
    figure(j)
    t_array = 1:1:time;
    anc=repmat(1:1:20, length(t_array), 1);
    plot3(t_array, anc, pi_array);
end

%verificar que para cada t soma de pi(t)= 1
for j = 1:3
    for i = 1:time
        sum_pi(j,i) = sum(pi_array(:,i));
    end
end

disp('Tempo para atingir o equil�brio para cada condi��o:');
condicao1 = eq_time(1)
condicao2 = eq_time(2)
condicao3 = eq_time(3)


disp('Provar que para cada instante de t a sompa de pi(t) = 1:');
sum1 = sum(sum_pi(1,:))
sum2 = sum(sum_pi(2,:))
sum3 = sum(sum_pi(3,:))

for j = 1:3
    figure(j)
    if j == 1
        title('Condi��o inicial -> �ncora 7');
    elseif j == 2
        title('Condi��o inicial ->�ncora 8');
    elseif j == 3
        title('Condi��o inicial -> igual probabilidade para todas as �ncoras');
    end
    xlabel('Tempo');
    ylabel('�ncora');
    zlabel('Probabilidade de ter o token');
end

%%
% Decidiu-se fazer tr�s simula��es para estudar o impacto de $\pi$(0) na
% evolu��o do sistema. Escolhemos para a primeira condi��o o token estar
% come�ar na �ncora 7, pois como conclu�do na quest�o 2.a esta � a �ncora
% que mais vezes tem o token por ter uma localiza��o central e v�rias
% liga��es, logo intuitivamente seria esta a condi��o inicial que
% conduziria a um tempo de equil�brio menor, como foi comprovado. Para a
% segunda situa��o escolhemos o oposto da primeira, ou seja, colocou-se o
% token na �ncora 8 que � a que possui menos vezes o token, logo ser� esta
% a que demora um maior tempo para atingir o equil�brio. Na terceira
% condi��o, todas as �ncoras t�m uma igual probabilidade de possuir o
% token, logo o tempo para atingir o equil�brio � algo no meio da situa��o
% inicial 1 e 2.
% Para provar que para cada instante de $\pi$(t) a soma das suas
% probabilidades � 1 decidiu-se mostrar o somat�rio das probabilidades para
% todo o tempo , pois como o tempo � de 500, se o somat�rio for 500 � uma boa
% indica��o que est� tudo bem. Analisando a matriz sum_pi() pode-se
% confirmar isso, mas decidiu-se apresentar o resultado desta maneira para
% n�o saturar a Command Window.


%% Exerc�cio 2.d

clear all
close all
clc

MarkovChain = load('MarkovChain.mat');

worse_P = MarkovChain.P;

better_P = worse_P;

%refazer liga�oes de P para melhor
better_P(2,6) = 0.3; better_P(2,4) = 0.35; better_P(2,13) = 0.35;
better_P(6,2) = 0.25; better_P(6,1) = 0.25; better_P(6,15) = 0.25; better_P(6,11) = 0.25;
better_P(3,19) = 0.5; better_P(3,12) = 0.5;
better_P(18,8) = 0.3; better_P(18,14) = 0.35; better_P(18, 16) = 0.35;
better_P(8,9) = 0.35; better_P(8,12) = 0.35; better_P(8,18) = 0.3;
better_P(20,5) = 0.25; better_P(20,1) = 0.25; better_P(20,7) = 0.25; better_P(20,14) = 0.25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcular as probabilidades com a nova matriz de transi��es

[vectors, values] = eig(better_P');
old_diference = 10000;
val_u = 1;
%encontra o valor proprio mais perto de 1
for i = 1:size(values,1)
    diference = abs(values(i,i) - 1);
    if diference < old_diference
        old_diference = diference;
        value_index = i;
    end
end
%vector a normalizar
vec_n = vectors(:,value_index);
normalize = vec_n/sum(vec_n);
figure(1)
bar(normalize);
title('Distribui��o de equil�brio da cadeia de Markov melhorada');
xlabel('N�mero do estado da cadeia de Markov');
ylabel('Probabilidade do token estar na posse do agente');
box on;
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimar a localiza��o da source

sidelength = 100;
N = 20;
n = 2;
M = 100;
a = [MarkovChain.nodePos(:,2) , MarkovChain.nodePos(:,3)]';
x = MarkovChain.sourcePos';
D = squareform(pdist([x zeros(size(x)) a]'));
d = D(1,3:end);                     % Source-anchor distances
an = D(2,3:end);                    % Anchor norms
P0 = 100;                    
mu = 0;
sigma = 1e-2;                       
A = zeros(M,4);
b = zeros(M,1);
ni_space = [-0.1:.0002:0.1];
ni = normpdf(ni_space, mu, sigma);

number_meas = normalize.*M; 
number_meas = round(number_meas);
j = 1;
for i = 1:N
    meas = number_meas(i);
    if meas > 0
        better_P(j) = (P0./(d(i).^2))*exp(ni(j));
        A(j,1) = -2*better_P(j)*a(1,i);
        A(j,2) = -2*better_P(j)*a(2,i);
        A(j,3) = -1;
        A(j,4) = better_P(j);
        b(j) = -better_P(j)*(an(i)^2)';
        j = j + 1;
        meas = meas-1;
    end
end

z = A\b;
% RLS formulation (one-shot)
RlsPar = struct('lam',1);
[e,w,RlsPar] = qrrls(A,b,RlsPar);
fprintf('Erro com %d mudan�as de agente na cadeia melhorada', M);
x_estimated = z(1:n);
error1 = norm(x-x_estimated)
%error1 = norm(z-w)
%x_estimated = [w(1) w(2)]';

figure(2)
plot(a'*[1; 1i],'o'); hold all
plot(x'*[1; 1i],'x'); plot(x_estimated'*[1; 1i],'s'); hold off
axis(sidelength*[0 1 0 1]); axis('square')
title('Estima��o da posi��o da source com cadeia de Markov melhorada');
legend('�ncoras', 'Posi��o real da source', 'Posi��o estimada da source com 10000 mudan�as');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%refazer liga��es de P para pior
worse_P = worse_P;

worse_P(7,1) = 0.1; worse_P(7,20) = 0.7; worse_P(7,16) = 0.1; worse_P(7,19) = 0.1;
worse_P(19,3) = 0.1; worse_P(19,4) = 0.6; worse_P(19,13) = 0.2; worse_P(19,7) = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcular as probabilidades com a nova matriz de transi��es

[vectors, values] = eig(worse_P');
old_diference = 10000;
val_u = 1;
%encontra o valor proprio mais perto de 1
for i = 1:size(values,1)
    diference = abs(values(i,i) - 1);
    if diference < old_diference
        old_diference = diference;
        value_index = i;
    end
end
%vector a normalizar
vec_n = vectors(:,value_index);
normalize = vec_n/sum(vec_n);
figure(3)
bar(normalize);
title('Distribui��o de equil�brio da cadeia de Markov piorada');
xlabel('N�mero do estado da cadeia de Markov');
ylabel('Probabilidade do token estar na posse do agente');
box on;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimar a localiza��o da source

sidelength = 100;
N = 20;
n = 2;
M = 100;
a = [MarkovChain.nodePos(:,2) , MarkovChain.nodePos(:,3)]';
x = MarkovChain.sourcePos';
D = squareform(pdist([x zeros(size(x)) a]'));
d = D(1,3:end);                     % Source-anchor distances
an = D(2,3:end);                    % Anchor norms
P0 = 100;                    
mu = 0;
sigma = 1e-2;                       
A = zeros(M,4);
b = zeros(M,1);
ni_space = [-0.1:.0002:0.1];
ni = normpdf(ni_space, mu, sigma);

number_meas = normalize.*M; 
number_meas = round(number_meas);
j = 1;
for i = 1:N
    meas = number_meas(i);
    if meas > 0
        worse_P(j) = (P0./(d(i).^2))*exp(ni(j));
        A(j,1) = -2*worse_P(j)*a(1,i);
        A(j,2) = -2*worse_P(j)*a(2,i);
        A(j,3) = -1;
        A(j,4) = worse_P(j);
        b(j) = -worse_P(j)*(an(i)^2)';
        j = j + 1;
        meas = meas-1;
    end
end

z = A\b;
% RLS formulation (one-shot)
RlsPar = struct('lam',1);
[e,w,RlsPar] = qrrls(A,b,RlsPar);
fprintf('Erro com %d mudan�as de agente', M);
x_estimated2 = z(1:n);
error2 = norm(x-x_estimated2)
%error2 = norm(z-w)
%x_estimated2 = [w(1) w(2)]';

figure(4)
plot(a'*[1; 1i],'o'); hold all
plot(x'*[1; 1i],'x'); plot(x_estimated2'*[1; 1i],'s'); hold off
axis(sidelength*[0 1 0 1]); axis('square')
title('Estima��o da posi��o da source com cadeia de Markov piorada');
legend('�ncoras', 'Posi��o real da source', 'Posi��o estimada da source com 10000 mudan�as');
%%
% Olhando para a matriz de transmiss�o P � possivel concluir algumas
% coisas. Em primeiro, destaca-se a olho a exist�ncia de dois grandes
% subgrupos, isto �, o grupo do lado esquerdo da liga��o 7-19 e o grupo do 
% lado direito desta liga��o. Ambos estes grafos t�m uma probabilidade de
% 0.25 de transitarem para o outro, o que mostra como pode ser dif�cil para
% o token passar de um lado para o outro, pois estes grafos s�o os �nicos
% pontos de liga��o de um sub-grupo para o outro. Mesmo dentro destes
% sub-grupos o token n�o percorre os grafos de maneira equitativa. A ponto
% de refer�ncia referimos por exemplo que a probabilidade do token passar
% do grafo 3 para o 12 � de 0.1, ou que para passar do 1 para o 6 � de 0.2,
% sendo que estes dois grafos s�o a �nica liga��o para alguns grafos.
% Portanto, � necess�rio refazer a cadeia de Markov para o token poder
% passar por mais agentes em vez de passar bastante tempo na posse dos
% mesmos agentes.
% Observando as primeiras duas figuras conlu�mos que a distribui��o dos
% estados mais e menos prov�vel est� bastante mais equilibrada que no
% exerc�cio 2.a, e que o erro da cadeia melhorada tamb�m � um pouco melhor,
% Decidimos reduzir o n�mero de transi��es de 10000 (2.b) para 100, de forma
% a podermos observar diferen�a na estima��o da source, pois para 10000
% transi��es ambas as cadeias davam o mesmo erro.

%% Exerc�cio 3a e 3b
part3Exe(1000, 1000, 0, 1);

%%
% Conseguimos verificar que a distribui��o dos estados e a frequ�ncia dos
% mesmos corresponde ao obtido na sec��o 2. Obtiv�mos o mesmo gr�fico de
% barras para a distribui��o de equilibrio da cadeia de Markov.
% Conseguimos concluir tamb�m que quanto mais alto o n�mero de runs ou de
% instantes de tempo que se quer analisar, melhores ser�o os resultados, e
% o tempo de converg�ncia reduz-se. Isto � uma das caracter�sticas do
% algoritmo de estima��o de Monte Carlo.
%

%% Exerc�cio 3c
part3Exe(1000, 1000, 1, 1);  %lambda = 1

%% Functions - part3Exe
type('part3Exe.m')

%% Functions - execute
type('execute.m')

%% Functions - PowrGen
type('PowrGen.m')

        

