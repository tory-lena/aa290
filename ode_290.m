function dx=ode_290(x, u, t, ut)

if ~exist('t', 'var')&& ~exist('ut', 'var')
    U=u;
else
    for k=1:length(u(1, :))
        U(:,k)=interp1(ut, u(:,k), t);
    end
end
clear u %cause u is vel in body x direction

dx=zeros(12,1);

p=x(1);
q=x(2);
r=x(3);
u=x(4);
v=x(5);
w=x(6);
phi=x(7);
theta=x(8);
psi=x(9);
X=x(10);
Y=x(11);
Z=x(12);

%constants
g=9.8;
area=88.95;

m=26421*.435;%6350;
w=54.70*.305;%18.92;
%h=2.2;
l=29.56*.305;%11.63;

I=157.1*diag([48120, 31640, 98475]); %10^5*diag([.6067, 1.7562, 2.2860]);%I=find_intertia(.5*w,h,l);I=m*I; %[Ix, 0, Ixz; 0, Iy, 0; Ixz, 0, Iz];
mg=m*g;

a=4;% distance ebtween centerline and turbine
Ix=I(1,1); Iy=I(2,2); Iz=I(3,3); Ixz=0;%(Ix+Iz)*.001;

rho_inf= 1.225;%At sea level and at 15 �C air has a density of approximately 1.225 kg/m3
p_inf= 1;%bar
T=15+273.5; %K

% R and T
Rot=(1/cos(theta))*[cos(theta), sin(theta)*sin(phi), sin(theta)*cos(phi);...
    0,cos(theta)*cos(phi), -cos(theta)*sin(phi);...
    0,sin(phi), cos(phi)];

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

%F=T*[0;0;m*g]; %add lift and drag, F_aero(pressure)
%Q=zeros(3,1); %resulting moms

eps=U(:, 2);% eps angle btw thrust and body x-axis?


c_l=find_cl(theta);
L=[0,0,-.5*rho_inf*(u^2+v^2+w^2)*c_l*area]'; L=T*L;
c_d=find_cd(theta, c_l);
D=[-.5*rho_inf*(u^2+v^2+w^2)*c_d*(.1*area), 0, 0]'; D=T*D;

% c_l=.5+(theta*(180/pi))*.1;
% L=[0,0,-.5*rho_inf*(u^2+v^2+w^2)*c_l*area]';
% c_d=.9;
% D=[-.5*rho_inf*(u^2+v^2+w^2)*c_d*(area*.2), 0, 0]';

F=L+D+[-mg*sin(theta) + U(:,1).*cos(eps); ...
    mg*cos(theta)*sin(phi); ...
    mg*cos(theta)*cos(phi) - U(:,1).*sin(eps)]; %U(1:3)    %L + D + [-mg*sin(theta) + U(1)*cos(eps); mg*cos(theta)*sin(phi); mg*cos(theta)*cos(phi) - U(1)*sin(eps)]
u3=U(:, 3);
M=[(ones(length(u3), 1)-2*u3).*U(:, 1).*(sin(eps))*a, U(:, 4), (ones(length(u3), 1)-2*u3).*U(:, 1).*cos(eps)*a]';

dx(1:3)= inv(I)*[M(1)-(Iz-Iy)*q*r - Ixz*q*p; ...
    M(2)-(Ix-Iz)*p*r-Ixz*(r^2-p^2);...
    M(3)-(Iy-Ix)*p*q+Ixz*q*r];

dx(4)=F(1)./m+(r*v-q*w);
dx(5)=F(2)./m+(p*w-u*r);
dx(6)=F(3)./m+(u*q-p*v);

dx(7:9)=Rot*[p;q;r];

dx(10:12)=inv(T)*[u;v;w];

end

function cl=find_cl(theta)
max=17;
if theta/pi*180 <=max
    cl=(theta/pi*180+5)*(.75/10);%.5+(theta*(180/pi))*.1;
else
    cl=1.7-(theta/pi*180-max)*(.7/13);%1.7-((theta-max)*(180/pi))*.1;
end
end

function cd=find_cd(theta, cl)
cd=(.2/100)*(theta/pi*180-5)^(2/3)+.01;
% if theta<=11;
%     coeff=.5*theta/pi*180;
% else
%     coeff=5-.005*(theta/pi*180)^2;
% end %Phanton 2 - cd0= 0.021!
% cd=cl/coeff;
end




