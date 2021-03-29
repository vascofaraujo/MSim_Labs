function [u] = blockFunc10(y, y_dot, y1)

k1 = 1/y1;
k2 = sqrt(2*k1);

in_subsys = -y;

if abs(in_subsys) <= y1
    out_subsys = (k1/k2)*in_subsys;
else 
    out_subsys = sign(in_subsys)*(sqrt(2*abs(in_subsys))-(1/k2));
end

eq =  out_subsys - y_dot;

% SIGN BLOCK

u = sign(eq);


end