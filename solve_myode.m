function [t,y] = solve_myode(A,T,y0,beta)
N=size(A,1);
tspan = [0 T];
% y0 = .5*ones(N,1);

% Call the ODE solver ode15s.
[t,y] = ode15s(@Amult,tspan,y0);

    % Define the ODE function.
    function dydt = Amult(t,y)
    dydt = A*y+beta;
    end
end