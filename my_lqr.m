function [K, t]=my_lqr(xd, ud, tau)
Q=zeros(12,12);% Q(10:12, 10:12)=eye(3);
R=.01*diag([1 0.01 0.01 .5]);
Qi=eye(12);
I=157.1*diag([48120, 31640, 98475]); 

[H, G]=HG(I);
[T P_flipped]=ode45(@(T,P_flipped)mRiccati(T, P_flipped, Q, R, xd, ud, tau, H, G), [tau(end) tau(1)], Qi);
P=flipud(P_flipped); 
t=flipud(T);

%get K
for i=1:length(t)
    Pi=reshape(P(i, :), [12 12]);
    s=t(i);
    for j=1:4
        ut(j)=interp1(tau, ud(:, j), s);
    end
    Bi=G(ut(1), ut(2), ut(3)); Ki=inv(R)*Bi'*Pi; K(i, :)=Ki(:)';
end
end

function dedt = mRiccati(t, e, Q, R, xd, ud, tau, H, G)

m=12; n=4;

clear i
for i=1:m
    if i<=n
        ut(i)=interp1(tau, ud(:, i), t);
    end
    xt(i)=interp1(tau, xd(:, i), t);
end
%xt=e+xd;
%ut=ud;%+v;??
%H(xt(1), xt(2), xt(3), xt(4), xt(5), xt(6), xt(7), xt(8), xt(9), xt(10), xt(11), xt(12), ut(1), ut(2), ut(3), ut(4), cl , cd, I(1,1), I(2,2), I(3,3), I(1,3));
%G(xt(1), xt(2), xt(3), xt(4), xt(5), xt(6), xt(7), xt(8), xt(9), xt(10), xt(11), xt(12), ut(1), ut(2), ut(3), ut(4), cl , cd, I(1,1), I(2,2), I(3,3), I(1,3));

cl=find_cl(xt(8)); cd=find_cd(xt(8));
A=real(H(xt(9), cd, cl, xt(1), xt(7), xt(2), xt(3), xt(8), xt(4), xt(5), xt(6)));
B=G(ut(1), ut(2), ut(3));

e=reshape(e, size(A));%Convert from "n^2"-by-1 to "n"-by-"n"
dedt = -(A'*e + e*A - e*B*inv(R)*B'*e + Q);
dedt = dedt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end


function cl=find_cl(theta)
max=17;
if theta/pi*180 <=max
    cl=(theta/pi*180+5)*(.75/10);
else
    cl=1.7-(theta/pi*180-max)*(.7/13);
end
end

function cd=find_cd(theta)
cd=(.2/100)*(theta/pi*180-5)^(2/3)+.01;
end