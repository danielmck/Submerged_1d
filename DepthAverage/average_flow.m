opengl software
vec = load("p_driven_slowing.txt");

N = 200; % number of discretisation points in z for each of tau, u
t_points = size(vec,1);
h = 4e-3; % layer height (m)
d=1.43e-4; % grain diameter (m)
phi_c=0.6; % Volume fraction
eta_f = 0.0010016; % Pa s
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 9.5; % deg
alpha = 0.001; % 1/Pa
dz = h/(N-0.5); % z spacing
z_pe = linspace(dz/2,h,N);
z_u = linspace(0,h-dz/2,N);

v_scale = sqrt(g.*h);
p_scale = rho_f.*g.*h;
t_scale = sqrt(h./g);
z_scale = h;
density_ratio = rho_p./rho_f;

p_e = vec(:,1:N)./p_scale;
% p_p = (p_b-p_e);
% phi = alpha*p_p;
u_f = vec(:,N+1:2*N)./v_scale;
u_p = vec(:,2*N+1:end)./v_scale;

u_p_ave = (sum(u_p,2)-(u_p(:,1)+u_p(:,end))/2)/N;
u_f_ave = (sum(u_f,2)-(u_f(:,1)+u_f(:,end))/2)/N;
p_e_ave = (sum(p_e,2)-(p_e(:,1)+p_e(:,end))/2)/N;
is_flow = (abs(u_p)>1e-5);
u_p_depth = sum(is_flow,2)+1;
u_p_depth(u_p_depth>N) = N;
u_p_flow_ave= zeros(N,1);
u_f_flow_ave= zeros(N,1);
p_e_flow_ave= zeros(N,1);
for j=1:t_points
    u_p_flow_ave(j) = (sum(u_p(j,1:u_p_depth(j)))-(u_p(j,1)+u_p(j,u_p_depth(j))/2))/N;
    u_f_flow_ave(j) = (sum(u_f(j,1:u_p_depth(j)))-(u_f(j,1)+u_f(j,u_p_depth(j))/2))/N;
    p_e_flow_ave(j) = (sum(p_e(j,1:u_p_depth(j)))-(p_e(j,1)+p_e(j,u_p_depth(j))/2))/N;
end
u_p_depth = u_p_depth .* h/N;

ani_num = 1;


    
nr_fr = t_points;
t_vals = (1:nr_fr)*0.1/t_scale(ani_num);
figure(1)
animation_name = 'depth_average_comp.gif';
for i = 2 : nr_fr
    subplot(2,1,1) 
    
    plot(p_e(i,:), z_pe/h);
    hold on
    plot(p_e_ave(i)*ones(200,1),z_pe/h,':k');
    plot([0 0 p_e_flow_ave(i) p_e_flow_ave(i)],[0 1-u_p_depth(i)/h 1-u_p_depth(i)/h 1],'--r')
    hold off
    legend("t="+num2str(t_vals(i)));
    xlim([0 1]);
    ylim([0 1]);
    xlabel("p_e");
    ylabel("z");

    subplot(2,1,2)
    plot(u_p(i,:), z_pe/h);
    hold on
    plot(u_p_ave(i)*ones(200,1),z_pe/h,':k');
    plot([0 0 u_p_flow_ave(i) u_p_flow_ave(i)],[0 1-u_p_depth(i)/h 1-u_p_depth(i)/h 1],'--r')
    hold off
    xlim([0,0.5])
    ylim([0 1])
    xlabel('u_p');
    ylabel("z");
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
      imwrite(imind,cm,animation_name,'gif', 'LoopCount',Inf);
    else
      imwrite(imind,cm,animation_name,'gif','WriteMode','append','DelayTime', 0.1);
    end
end