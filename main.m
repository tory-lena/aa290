clear all

%% set up
[name, tf, Q, R, m, n, N, x_con, u_con, xyz_0, xyz_f, v_i]=setup();

%% solve
[Prob, s] =def_MP(name, tf, Q, R, m, n, N, x_con, u_con, xyz_0, xyz_f, v_i);
[x_star, u_star, tau]=comp_MP(Prob, s, x_con, m, n, N, 1000);

%% run ----------------------------------  controller
%initialize
[k, x]= simctrl(x_star, u_star, tau);
x_curr=x; t_curr=tau(k);
clear x

tf=tau(end);
step_l=k;
T=tau(step_l+1);

N=6; Tp=T;
while any(abs(x_curr(end, 10:12))>=.5)
    [Prob_i, s_i] =def_NMPC2('one_step', m, n, N, x_con, u_con, x_curr(end, :), ...
         x_star(k+1:end, :), tau(k+1:end), Tp); 
    [x_e, u_e, t_e, flag]=comp_NMPC(Prob_i, s_i, Tp, x_con, m, n, N, 100);

    %use control
    [t, x]=ode45(@(t, x)ode_290(x, u_e, t, t_e), [t_curr t_curr+T], x_curr(end, :));
    x_curr=vertcat(x_curr, x); t_curr=t_curr+T;

    k=interp3(x_star(:, 10), x_star(:, 11), x_star(:,12), linspace(1, length(tau), length(tau)), x_curr(end, 10), x_curr(:, 11), x_curr(end, 12));
    k=round(k);
    clear Prob_i s_i x_e u_e t_e
end

% %recompute path - compensate deviation
% for i=1:12
%     N=6; Tp=2*T;
%     [Prob_i, s_i] =def_NMPC('one_step', tf, m, n, N, x_con, u_con, x_curr(end, :), ...
%         x_star((step_l)*i+1:end, :), tau((step_l)*i+1:end), t_curr+Tp, t_curr);
%     [x_e, u_e, t_e, flag]=comp_NMPC(Prob_i, s_i, Tp, x_con, m, n, N, step_l);
%     if flag==0
%         clear flag
%         [Prob_i, s_i] =def_NMPC('one_step', tf, m, n, N-2, x_con, u_con, x_curr(end, :), ...
%             x_star((step_l)*i+1:end, :), tau((step_l)*i+1:end), t_curr+Tp, t_curr);
%         [x_e, u_e, t_e, flag]=comp_NMPC(Prob_i, s_i, Tp, x_con, m, n, N-2, step_l);
%         if flag==0
%             break
%         end
%     end
%     
%     %use control
%     [t, x]=ode45(@(t, x)ode_290(x, u_e, t, t_e), [t_curr t_curr+T], x_curr(end, :));
%     x_curr=vertcat(x_curr, x); t_curr=t_curr+T;
% end


% hor=[40, 20];
% x_curr=online(m, n, x_con, u_con, x, u, tau, hor);


%% like what u see?
%file_name=''; % do not use .txt, name only!
%save_f(file_name, x, u);
