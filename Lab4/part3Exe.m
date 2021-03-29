function part3Exe(n_runs, steps, slow, lambda)
info = importdata('MarkovChain.mat');        %information about the problem
P = info.P;
info.nodePos(:,1) = [];
npos = info.nodePos';   %a
spos = info.sourcePos'; %x
n_anchor = length(info.P);                           %number of steps on the Markov
tok_node = ceil((rand()*19)+1);                      %generate random initial condition


[xe, esterror, error_zw, error, vec_w, statepath, poserror] = execute(spos, npos, P, tok_node, n_runs, steps, slow, lambda);
mode_state = mode(statepath(:));        %most frequent state

for i =1:length(npos)                   %creates ocurrences vector
    occur(i) = nnz(statepath==i);
end

occur = occur/sum(occur);

% Plots
figure(1)
plot(1:steps, error, '-b');
grid on
title('Qrrls Error');
xlabel('Time Steps');
ylabel('Error');

figure(2)
plot(1:steps, error_zw, '-r');
grid on
title('Z-W Error');
xlabel('Time Steps');
ylabel('Error');

figure(3)
plot(1:steps, poserror, '-g');
grid on
title('Source Position - W Error');
xlabel('Time Steps');
ylabel('Error');

figure(4)
bar(occur);
grid on
title('Ocurrence probability bar chart');
xlabel('States');
ylabel('Probability');
end

