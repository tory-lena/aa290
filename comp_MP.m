function [x, u, tau]=comp_MP(Prob, t, x_con, m, n, N, incr )

%% solve

Prob = ProbCheck(Prob,'snopt');

Prob.SOL.optPar(30)=80000; %-------- was active!!!
% Prob.SOL.optPar(10)=10^(-3);

% Prob.SOL.optPar(27)=10^(-20);
% constraint, iabort function

Result = snoptTL(Prob)
sol=Result.x_k;
tf=sol(end);

[x, u, o, V, tau]=check_res(sol, tf, x_con, m, N, n, t, incr); %output x0 form ode missing - plot!!

%% plot

%landing zone
lz=[-.5 -.5 0; .5 -.5 0; .5 .5 0; -.5 .5 0; -.5 -.5 0; -.5 -.5 .1; .5 -.5 .1; .5 .5 .1; -.5 .5 .1];

%fig
figure()
plot3(x(:, 10), -x(:, 11), -x(:, 12), 'b',...
    sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), 'ko',...
    lz(:, 1), lz(:, 2), lz(:, 3), 'r--')
grid on

figure()
quiver3(sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), o(:, 1), o(:, 2), o(:, 3), 'b')
hold on
quiver3(sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), V(:, 1), V(:, 2), V(:, 3), 'r')

% figure()
% plot3(x(:, 10), -x(:, 11), -x(:, 12), 'b',...
%     sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), 'ko',...
%     lz(:, 1), lz(:, 2), lz(:, 3), 'r--')
% grid on
% hold on
% plot3(x0(:, 10), -x0(:, 11), -x0(:, 12), 'b')


end