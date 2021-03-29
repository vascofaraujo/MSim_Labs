function [u] = blockFunc8(y, y_dot)

in_subsys = -y;

out_subsys = sign(in_subsys)*sqrt(2*abs(in_subsys));

eq =  out_subsys - y_dot;

% SIGN BLOCK

u = sign(eq);


end

