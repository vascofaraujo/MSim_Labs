function [xe, poserr, error_zw, error, vec_w, statepath, poserror] = execute(spos, npos, P, tok_node, n_runs, steps, slow, lambda)

[Powr, an] = PowrGen(spos, npos);

h = waitbar(0,'Please wait...');
for n = 1:n_runs
    for c = 1:steps
        if slow == 1  %if slow mode activated
            spos = spos+0.1; %moves per interation
           [Powr, an] = PowrGen(spos, npos);  %new powr
        end
        tok_node = find(cumsum([P(tok_node,:)]) > rand,1,'first');          %next node
        statepath(n,c) = tok_node;
        % Localize source by least-squares
        A(c,:) = [-2*repmat(Powr(1,tok_node),[2 1]).*npos(:,tok_node); -1; Powr(1,tok_node)]';
        b(c,:) = (-Powr(1,tok_node)*norm(an(tok_node))^2)';
        
        z = A\b;
        xe(1,c) = z(1);
        xe(2,c) = z(2);
        poserr(c) = norm(spos-xe);             %estimation error
        
        % RLS formulation (incremental)
        RlsPar = struct('lam',lambda);
        
        [e,w,RlsPar] = qrrls(A(c,:),b(c),RlsPar);
        error(c) = abs(e);
        vec_w(1,c) = w(1);
        vec_w(2,c) = w(2);
        error_zw(c) = norm(z-w);
        poserror(c) = norm(spos-w(1:2));
        
    end
    waitbar(n/n_runs);
end
close(h);

end

