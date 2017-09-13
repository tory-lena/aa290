function c=dyn_con(xu, D, m, n, N, tf)


%% filter state and control var
X=xu(1:m*(N+1)); U=xu(m*(N+1)+1:(m+n)*(N+1));
if tf<=0%length(xu)>(m+n)*(N+1)
    tf=xu(end);
end

for i=0:N
    for j=1:m
    x(i+1, j)=X(i*m+j);
    end
    for k=1:n
    u(i+1, k)=U(i*n+k);
    end
end

%% dynamic sys constraints

for i=0:N
    dx(i*m+1:m*(i+1))=ode_290(x(i+1, :)', u(i+1, :));
end
clear i

c=(2/tf)*D*X-dx';

% c_av=mean(abs(c));
% if c_av<=.0001
%     error_paths(xu, N, m)
% end

end