function [x, u, tau, x_sol, u_sol]=comp_MP(Prob, t, x_con, m, n, N, incr )

%% solve

Prob = ProbCheck(Prob,'snopt');

 Prob.SOL.optPar(30)=80000; %-------- was active!!!
 T=10^(-6);
 Prob.SOL.optPar(10)=sqrt(T);
 Prob.SOL.optPar(41)=T; clear T
% constraint, iabort function

Result = snoptTL(Prob)
sol=Result.x_k;
tf=sol(end);

[x, u, o, V, tau, x0, x_sol, u_sol]=check_res(sol, tf, x_con, m, N, n, t, incr); %output x0 form ode missing - plot!!
T=(flip(t)+1)/2*tf;

%% plot

%landing zone
lz=[-.5 -.5 0; .5 -.5 0; .5 .5 0; -.5 .5 0; -.5 -.5 0; -.5 -.5 .1; .5 -.5 .1; .5 .5 .1; -.5 .5 .1];
a=400;

%fig
figure()
plot3(x(:, 10), -x(:, 11), -x(:, 12), 'b',...
    sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), 'ko',...
    lz(:, 1), lz(:, 2), lz(:, 3), 'r--')
axis([-a a -a a -a a])
grid on

figure()
quiver3(sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), o(:, 1), o(:, 2), o(:, 3), 'b')
hold on
quiver3(sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), V(:, 1), V(:, 2), V(:, 3), 'r')
axis([-a a -a a -a a])

figure()
plot3(x(:, 10), -x(:, 11), -x(:, 12), 'b',...
    sol(10:m:m*(N+1)), -sol(11:m:m*(N+1)), -sol(12:m:m*(N+1)), 'ko',...
    lz(:, 1), lz(:, 2), lz(:, 3), 'r--')
grid on
hold on
plot3(x0(:, 10), -x0(:, 11), -x0(:, 12), 'b')
axis([-a a -a a -a a])


%Title=['Thrust', 'p in %', 'eps', 'M_y'];
figure()
for j=1:n
    subplot(4,1,j)
    plot(tau, u(:,j),'b-',  T, u_sol(:, j),'b*')
    title(['u_' num2str(j) ])
end

end