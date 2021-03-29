

function [Pb] = generatePulse(t, beta)
    
    for i = 1: length(t)
        if t(i) <= - beta/2 - 0.5                               %1ª secçao
            Pb(i) = 0;
        elseif t(i) > - beta/2 - 0.5 && t(i) <= -0.5              %2ª secção
            Pb(i) = (2*(t(i)^2)+2*t(i)+0.5)/(beta^2) + (2*t(i)+1)/beta + 0.5;
        elseif t(i) > -0.5 && t(i) < beta/2 - 0.5               %3ª secção
            Pb(i) = -(2*(t(i)^2)+2*t(i)+0.5)/(beta^2) + (2*t(i)+1)/beta + 0.5;
        elseif t(i) >= beta/2 - 0.5 && t(i) <= - beta/2 + 0.5       %4ª secção
            Pb(i) = 1;
        elseif t(i) > - beta/2 + 0.5 && t(i) < 0.5                %5ª secção
            Pb(i) = (-2*(t(i)^2)+2*t(i)-0.5)/(beta^2) + (-2*t(i)+1)/beta + 0.5;
        elseif t(i) >= 0.5 && t(i) < beta/2 + 0.5                %6ª secção
            Pb(i) = (2*(t(i)^2)-2*t(i)+0.5)/(beta^2) + (-2*t(i)+1)/beta + 0.5;
        elseif t(i) >= beta/2 + 0.5                           %7ª secção
            Pb(i) = 0;
        end
    end
end
       


