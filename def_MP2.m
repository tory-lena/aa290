function [Prob, t] = def_MP2(Name, tf, Q, R, m, n, N, x_con, u_con, xyz_0, xyz_f, v_i)

if tf<=0
    E=abs(tf);
end
 
%state const
x_U=x_con(:,2); 
x_L=x_con(:,1);

x_Urun=x_con(:, 4); 
x_Lrun=x_con(:, 3);

%control const
u_U=u_con(:,2);
u_L=u_con(:,1);

u_Urun=u_con(:, 4); 
u_Lrun=u_con(:, 3);

%time
max_t=200;

%joint const
run_n=N/10;
xu_U=[kron(ones(N+1-run_n, 1),x_U); kron(linspace(1,.1,run_n)', x_Urun); kron(ones(N+1-run_n, 1), u_U); kron(ones(run_n, 1), u_Urun); max_t];
xu_L=[kron(ones(N+1-run_n, 1),x_L);kron(linspace(1,.1,run_n)', x_Lrun); kron(ones(N+1-run_n, 1), u_L); kron(ones(run_n, 1), u_Lrun); 0]; %check here impose Upper and lower for x direction (landing strip at run_n)
% xu_U((N-run_n+1)*m+12)=xu_L((N-run_n+1)*m+12)+1;
% xu_U((N-run_n+1)*m+10)=xu_L((N-run_n+1)*m+10)+10;

% get t, w
[t, w]=chebyshev(N);
t=flip(t);

% get D
D=ChebyshevDiffMatrix(N, t);
D=kron(D, eye(m));

%cost function
Q_w=kron(diag(w),Q);
R_w=kron(diag(w),R);
F=blkdiag(Q_w, R_w);

f=@(xu) E*xu(end)+(xu(end)/2)*xu(1:end-1)'*F*xu(1:end-1) ;
G=@(xu) xu(end)*F*xu(1:end-1);
H=@(xu) xu(end)*F;

%dyn constraints
c= @(xu) dyn_con(xu, D, m, n, N, tf);
c_U=zeros(m*(N+1),1);
c_L=c_U;

%initial guess
yaw_f=.0*pi;%+pi
yaw_i=-.0*pi;
r=round(N/4);

%xu0=[zeros((m+n)*(N+1), 1); 10];

xu0=[ones((m+n)*(N+1), 1)*.5;3];

xu0(1:m:m*(N+1))=zeros((N+1),1);
xu0(2:m:m*(N+1))=-xu0(2:m:m*(N+1));
xu0(3:m:m*(N+1))=(yaw_f/50)*linspace(0, 1, (N+1))';

xu0(4:m:m*(N+1))=v_i(1)*linspace(1, 0, (N+1))';
%xu0(5:m:m*(N+1))=;
xu0(6:m:m*(N+1))=20*[linspace(0,1, 2*r)'; linspace(1,0, (N+1)-2*r)'];

xu0(7:m:m*(N+1))=zeros((N+1),1);
%xu0(8:m:m*(N+1))=[zeros(1, r), .2*pi*ones(1, (N+1)-2*r),zeros(1, r)]';
xu0(9:m:m*(N+1))=yaw_i*[linspace(1, 0, (N+1)-r), ones(1, r)]';

xu0(10:m:m*(N+1))=[ones((N+1), 1)*xyz_0(1)+(xyz_f(1)-xyz_0(1))*[linspace(0, 1,(N+1)-r),ones(1,r)]'];% linspace(.8, 1, r)]'];
xu0(11:m:m*(N+1))=[ones((N+1), 1)*xyz_0(2)+(xyz_f(2)-xyz_0(2))*[zeros(1,r), linspace(0, 1,(N+1)-2*r), ones(1,r)]'];%, linspace(1, 1, r)]'];
xu0(12:m:m*(N+1))=[ones(4, 1)*xyz_0(3); ones((N+1)-4, 1)*xyz_0(3)-(xyz_0(3)-xyz_f(3))*[linspace(0, 1, (N+1)-8), ones(1,4)]'];

xu0(m*(N+1)+1:n:end-1)=.5*(u_U(1)-u_L(1))*ones((N+1),1);
xu0(m*(N+1)+2:n:end-1)=.5*(u_U(2)-u_L(2))*ones((N+1),1);
xu0(m*(N+1)+3:n:end-1)=.5*(u_U(3)-u_L(3))*ones((N+1),1);
xu0(m*(N+1)+4:n:end-1)=.5*(u_U(4)-u_L(4))*ones((N+1),1);


%equality constraints
Aeq=zeros(28, length(xu0)); %p q r u v w phi theta psi x y z 
Aeq(1:12, 1:12)=eye(12);
Aeq(13:24, m*N+1: m*N+12)=eye(12);
Aeq(25:28, m*(N+1)+n*N+1:m*(N+1)+n*N+4)=eye(4);
u_f=[0 0 0 0];u_f_dev=[10 pi*.25 .01 5];

mean=[0 0 0, v_i(1) v_i(2) v_i(3), 0 0 yaw_i, xyz_0, 0 0 0, 0 0 0, 0 0 yaw_f, xyz_f, u_f]'; 
dev=[0.1 0.1 0.1, .1 .1 .1, .05 .05 .05, .005 .005 .005, 0.1 0.1 0.1, 1 1 1, .01 .01 .01, 1 5 .1, u_f_dev]';
b_L=mean-dev;
b_U=mean+dev;

%PB defintion

%options.stop=iAbort();

Prob = conAssign(@(xu) xu(end), [], [], [], xu_L, xu_U, Name, xu0, ...
                     [], 0, Aeq, b_L, b_U, c, [], [], [], c_L, c_U, ...
                     [],[],[],[]); %n(f, G, H OR  %@(xu) xu(end), [], []

end