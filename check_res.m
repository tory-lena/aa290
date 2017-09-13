function [x, u, o, V, tau]=check_res(sol, tf, x_con, m, N, n, t, incr)

if tf==-1
tf=sol(end);
end

X=sol(1:m*(N+1)); U=sol(m*(N+1)+1:(m+n)*(N+1));
for j=1:m
    x(:, j)=X(j:m:end)';
end
for k=1:n
    u(:, k)=U(k:n:end)';
end
clear i j k
x_sol=x;
u_sol=u;

%% satisfy constraints?

flag=1;
for i=1:N+1
    if all(x(i, :)'-x_con(:, 1)>=0) && all(x(i, :)'-x_con(:, 2)<=0)
        continue
    else
        flag=0; %disp('not satistfying constraints')
        break
    end
end
clear i

if flag==1
    disp('valid path')
else
    disp('invalid path')
end


%% plot

%get orientation
for i=1:N+1
    
phi=x(i,7);
theta=x(i,8);
psi=x(i,9);

u=x(i,4);
v=x(i,5);
w=x(i,6);

T2=[cos(theta), 0, -sin(theta);...
    0,1,0;...
    sin(theta), 0,cos(theta)];
T3=[cos(psi), sin(psi), 0;...
    -sin(psi), cos(psi), 0; ...
    0,0,1];
T1=[1,0,0;...
    0,cos(phi), sin(phi);...
    0,-sin(phi), cos(phi)];
T=T1*T2*T3;

a=pi;
M=[1 0 0;...
    0 cos(a) sin(a);...
    0 -sin(a) cos(a)];

o(i,:)=M*inv(T)*[1;0;0];
V(i,:)=M*inv(T)*[u;v;w];
end

%% retrace path with computed control input

% x_ref(1,:)=x(1,:);
% for i=1:N+1
%     dx(i, :)=ode_290(x_ref(i,:), u(i,:))';
%     x_ref(i+1,:)=x_ref(i,:)+dx(i, :);
% end


%% give entire smoothened path

dt=tf/incr;
tau=0:dt:tf;
s=(2*tau-tf)/tf;

L = lagrange(length(s)-1,N,s,t);

clear x u
for i=1:m
    x(:, i)=flip(X(i:m:end)'*L)';
end
for i=1:n
    u(:,i)=flip(U(i:n:end)'*L)';
end


%% ode
% x0(1,:)=x(1,:);
% clear k
% for k=1:incr
% [t_r, x_r]= ode45(@(t_r, x_r)ode_290(x_r, u(k, :)), [tau(k) tau(k+1)], x0(k, :));
% x0(k+1,:)=x_r(end, :);
% end

end