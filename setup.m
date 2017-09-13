function [name, tf, Q, R, m, n, N, x_con, u_con, xyz_0, xyz_f, v_i]=setup()

m=12;
n=4;%6;
N=200 %30

Q=zeros(m);
R=zeros(n); %R(3,3)=10;

tf=-1;
name='Aircraft';
c=343; % m/s

%F_max=15*[0, 1; -1, 1; -1, 1];
pm=[-ones(3, 1), ones(3,1)];
M_max=1000*[-1,1];
max_rate=.5*pi*pm;
max_angle=pi*[-.5,.5; -45/180, 25/180; -1, 1]; %theta (-5,25)
max_pos=20000*[-1, 1; -1, 1; -1, 0]; %z is upside down!! [-Inf, inf]
max_vel=c*[-.9 .9; -.3 .3; -.5 .5];%[-1 c;-c c; -c c];%[--1 .9*c; -.3*c .3*c; -.5*c .5*c]

xyz_0=10^(3)*[-15 15 5]; xyz_0=rot(xyz_0);% [-10 , 1, 1];%1000*[-10 0 1];
xyz_f=[0 0 0]; xyz_f=rot(xyz_f);
v_i=.8*c*[1 0 0]; v_i=rot(v_i); %.8 Ma

u_con=[0, 10000; -pi/2, pi/2;0, 1; M_max];%[F_max; M_max];
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