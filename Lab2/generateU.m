function [ut, t] = generateU(T,alpha,beta,U1,U2,n1,n2)
T_1 = T/(alpha + 1);
T_2 = (T_1*alpha);
miu_1 = T_1/(1 + beta);
miu_2 = T_2/(1 + beta);

n_total = n1 + (n2-1);
t = linspace(0, T, n_total);


pb1 = generatePulse((t-(T_1/2))/miu_1, beta);
u1 = -U1*pb1;

pb2 = generatePulse((t-T_1-(T_2/2))/miu_2, beta);
u2 = U2*pb2;

ut = u1 + u2;
end

