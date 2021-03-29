

function [Pb] = generatePulse(t, beta)
    
    for i = 1: length(t)
        if t(i) <= - beta/2 - 0.5                               %1� sec�ao
            Pb(i) = 0;
        elseif t(i) > - beta/2 - 0.5 && t(i) <= -0.5              %2� sec��o
            Pb(i) = (2*(t(i)^2)+2*t(i)+0.5)/(beta^2) + (2*t(i)+1)/beta + 0.5;
        elseif t(i) > -0.5 && t(i) < beta/2 - 0.5               %3� sec��o
            Pb(i) = -(2*(t(i)^2)+2*t(i)+0.5)/(beta^2) + (2*t(i)+1)/beta + 0.5;
        elseif t(i) >= beta/2 - 0.5 && t(i) <= - beta/2 + 0.5       %4� sec��o
            Pb(i) = 1;
        elseif t(i) > - beta/2 + 0.5 && t(i) < 0.5                %5� sec��o
            Pb(i) = (-2*(t(i)^2)+2*t(i)-0.5)/(beta^2) + (-2*t(i)+1)/beta + 0.5;
        elseif t(i) >= 0.5 && t(i) < beta/2 + 0.5                %6� sec��o
            Pb(i) = (2*(t(i)^2)-2*t(i)+0.5)/(beta^2) + (-2*t(i)+1)/beta + 0.5;
        elseif t(i) >= beta/2 + 0.5                           %7� sec��o
            Pb(i) = 0;
        end
    end
end
       


