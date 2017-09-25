function [k, x, count]=simctrl(x_star, u_star, tau, err, P, t, x_curr)

clear k
count=[0];

for k=1:(length(tau)-1)
    clear K
    if ~exist('P', 'var') && ~exist('t', 'var') && ~exist('x_curr', 'var')
        x0(1,:)=x_star(1,:);
        u=u_star;
    else
        x0(1,:)=x_curr;
        for i=1:48 %4*12
            K(1, i)=interp1(t, P(:, i), tau(k));
        end
        K=reshape(K, [4 12]);
        u(k, :)=-K*(x0(k,:)'-x_star(k, :)')+u_star(k, :)';
    end
    
    if k>=2 && e(k-1)>=err
        start=x_star(k,:);
        count=[count, k];
    else
        start=x0(end,:);
    end
    
    
    [t_r, x_r]= ode45(@(t_r, x_r)ode_290(x_r, u(k, :) ), [tau(k) tau(k+1)], start);%x0(end, :));%x_star(k, :));
    x0(k+1,:)=x_r(end, :);
    e(k)=norm(x0(k+1, 10:12)-x_star(k+1, 10:12));
%     if e(k)>=err*.1
%         s=find_closest(x0(k,:), x_star); k=s; clear s
%     end
%         break
% %     elseif any(abs(x0(k, :)-x_star(k, :))>=9)
% %         break
%     end
    if x0(k, 12)>=2 && k>=30
        break
    end
        
    
end
%t=tau(k);
x=x0(1:k, :);
x=real(x0);
count=count(2:end);
disp(['Has been resetted' num2str(length(count)) 'times.'])

end

function j=find_closest(x, x_star)
x=x(1, 10:12); x_r=x_star(:, 10:12);

for i=1:length(x_r)
    err(i)=norm(x-x_r(i, :));
end
[E, j]=min(err);

end