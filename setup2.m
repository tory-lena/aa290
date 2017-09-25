function [name, tf, Q, R, m, n, N, x_con, u_con, xyz_0, xyz_f, v_i]=setup2()

m=12;
n=4;%6;
N=40 %30

Q=zeros(m);
R=zeros(n); %R(3,3)=10;

tf=-1;
name='Aircraft';
c=343; % m/s

run_l=.305*265; %feet to m -> .305*(-)
run_w=.305*30;  %on each side of centerline
run_h=.305*50;  

a=2;

%F_max=15*[0, 1; -1, 1; -1, 1];
pm=[-ones(3, 1), ones(3,1)];
M_max=10000*[-1,1];
max_rate=[.5*pi*pm, .1*pi*pm];
max_angle=pi*[-.5,.5, -.01, .01; -30/180, 30/180, -3/180, 3/180; -1, 1, -.01, .01]; %theta (-5,25)
max_pos=[-10^(a+1), 10^(a+1), -run_l, 0; -10^(a+1), 10^(a+1), -run_w, run_w; -10^(a+1), 0, -run_h, 0]; %z is upside down!! [-Inf, inf]
max_vel=c*[-.9 .9 -.1 .1; -.3 .3 -.1 .1; -.5 .5 -.1 .1];%[-1 c;-c c; -c c];%[--1 .9*c; -.3*c .3*c; -.5*c .5*c]

xyz_0=10^(a)*[-2 0 1]; xyz_0=rot(xyz_0);% -5; -10, -2
xyz_f=[0 0 0]; xyz_f=rot(xyz_f);
v_i=.05*c*[1 0 0]; v_i=rot(v_i); %.2 Ma .05 Ma

u_con=[0, 100000, 0, 100; -pi/2, pi/2, -pi/4, pi/4; 0, 1, 0, 1; M_max, M_max];%[F_max; M_max] ->*10;
x_con=[max_rate; max_vel; max_angle; max_pos];% p q r u v w phi theta psi x y z

end

%% def
function y=rot(x)
a=pi;
M=[1 0 0;...
    0 cos(a) sin(a);...
    0 -sin(a) cos(a)];
y=(M*x')';
end