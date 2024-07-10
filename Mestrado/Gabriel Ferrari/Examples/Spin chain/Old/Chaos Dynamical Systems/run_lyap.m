[T,Res]=lyapunov(3,@duffing_map,@ode45,0,0.5,200,[0 1 1],0);
plot(T,Res);
title('Dynamics of Lyapunov exponents');
xlabel('Time'); ylabel('Lyapunov exponents');
