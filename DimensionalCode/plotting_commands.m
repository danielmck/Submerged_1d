opengl software
vec = load("dil9_5deg_vsmall.txt");
N=200;
h = 4e-3;
d=1.43e-5;
dz = h/(N-0.5); % z spacing
phi_c=0.6; % Volume fraction
eta_f = 0.0010016; % Pa s
z = linspace(dz/2,h,N);
g=9.81; % m/s^2
rho_f = 1000; % kg/m^3
rho_p = 2500; % kg/m^3
theta = 9; % deg
alpha = 0.001; % 1/Pa

rho = rho_p*phi_c + rho_f*(1-phi_c);
rho_d = rho_p*phi_c - rho_f*(1-phi_c);
buoyancy = -(rho_p-rho_f)*g*phi_c*cosd(theta);

dz = h/(N-0.5); % z spacing
z_pe = linspace(dz/2,h,N);
z_u = linspace(0,h-dz/2,N);
p_b = (rho_p-rho_f)*g*phi_c*cosd(theta)*(h-z_pe);
sus_frac = 0.5;

p_e = vec(:,1:N)';
phi = vec(:,N+1:2*N)';
u_f = vec(:,2*N+1:3*N)';
u_p = vec(:,3*N+1:end)';
% pe_vals = s_frac*p_b;
% pp_vals = (1-s_frac)*p_b;
% phi_vals = alpha*pp_vals;
% wp_vals = pe_vals;
% for j = 1:6
%     pp_vals(:,j) = p_b'-pe_vals(:,j);
%     phi_vals(:,j) = alpha*pp_vals(:,j);
%     test = diff(pe_vals(:,j))./dz;
%     wp_vals(:,j) = [test' 0]'./beta_fn(phi_vals(:,j));
% end
% beta = beta_fn(phi_vals);
% plot(beta(:,6),phi_c+phi_vals(:,6));
% subplot(3,1,1);

% plot(p_e(:,5),z);
% hold on
% plot(p_e(:,6),z);
% plot(p_e(:,7),z);
% plot(p_e(:,8),z);
% plot(p_e(:,9),z);
% plot(p_e(:,10),z);
% plot(0.46*p_b,z,"--k");
% legend("28s","30s","32s","34s","36s","40s","Critical Fluid Pressure");
% xlabel("Excess Pressure");
% ylabel("z");
% hold off

% plot(u_p(:,5),z);
% hold on
% plot(u_p(:,6),z);
% plot(u_f(:,7),z);
% plot(u_p(:,8),z);
% plot(u_p(:,9),z);
% plot(u_p(:,10),z);
% legend("28s","30s","32s","34s","36s","40s");
% xlabel("u_p");
% ylabel("z");
% hold off

% plot(u_f(:,5),z);
% hold on
% plot(u_f(:,6),z);
% plot(u_f(:,7),z);
% plot(u_f(:,8),z);
% plot(u_f(:,9),z);
% plot(u_f(:,10),z);
% legend("28s","30s","32s","34s","36s","40s");
% xlabel("u_f ");
% ylabel("z");
% hold off

% Define number of frames
nr_fr = 1001;
% Initialize matrix using 'moviein'
% frames = moviein(nr_fr); 

% Generate frames with any plotting function. 
% We use a cosine with variable frequency.



figure(1)
filename = 'dilative_evolution_vsmall.gif';
v_scale = sqrt(g.*h);
p_scale = rho_f.*g.*h;
t_scale = sqrt(h./g);
z_scale = h;
t_vals = (1:nr_fr)*10/t_scale;

% break_t = (0.085+177*0.000025)/t_scale;
% plot(u_p(:,178)/v_scale, z_u'/h);
% ylim([0 1])
% xlabel("u_p");
% ylabel("z");
% legend("t="+num2str(break_t));
% saveas(gcf,"breaking_point_up.png")

for i = 1 : nr_fr    
    subplot(2,2,1)
    plot(p_e(:,i)./p_scale, z_u/h);
    legend("t="+num2str(t_vals(i)));
    xlim([-0.5 0]);
    ylim([0 1]);
    xlabel("p_e");
    ylabel("z");
    
    subplot(2,2,2)
    plot(phi(:,i), z_u/h);
    xlim([-0.01,0])
    ylim([0 1])
    xlabel('\phi-\phi_c');
    ylabel("z");
    
    subplot(2,2,3)
    plot(u_f(:,i)/v_scale, z_u/h);
    xlim([0,0.03])
    ylim([0 1])
    xlabel("u_f");
    ylabel("z");
    
    subplot(2,2,4)
    plot(u_p(:,i)/v_scale, z_u/h);
    xlim([0,0.03])
    ylim([0 1])
    xlabel("u_f");
    ylabel("z");
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
      imwrite(imind,cm,filename,'gif', 'LoopCount',Inf);
    else
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0.01);
    end
end




% u_diff = u_f_vals - u_p_vals;
% save('u_f_out.txt', 'u_f_vals','-ascii')
% subplot(3,1,1);
% plot(u_diff(1,:),z);
% 
% hold on
% plot(u_diff(3,:),z);
% plot(u_diff(5,:),z);
% plot(u_diff(7,:),z);
% plot(u_diff(9,:),z);
% plot(u_diff(11,:),z);
% xlabel("drag term");
% ylabel("z");
% hold off
% legend("0s","2s","4s","6s","8s","10s");
% drag_mult = (1-phi_c)^2.*beta_fn(alpha*(1-sus_frac)*p_b);
% grav = rho_f*(1-phi_c)*g*sind(theta);
% tau_f = [0 eta_f*diff([diff(u_f_vals(:,1)') 0])./(dz^2)];
% drag_term = drag_mult.*(u_f_vals(:,1)'-u_p_vals(:,1)');
% % plot(tau_f-drag_term+grav,z_u);
% 
% hold on
% tau_f = [0 eta_f*diff([diff(u_f_vals(:,2)') 0])./(dz^2)];
% drag_term = drag_mult.*(u_f_vals(:,2)'-u_p_vals(:,2)');
% plot(tau_f+grav,z_u);
% tau_f = [0 eta_f*diff([diff(u_f_vals(:,3)') 0])./(dz^2)];
% drag_term = drag_mult.*(u_f_vals(:,3)'-u_p_vals(:,3)');
% plot(tau_f+grav,z_u);
% tau_f = [0 eta_f*diff([diff(u_f_vals(:,4)') 0])./(dz^2)];
% drag_term = drag_mult.*(u_f_vals(:,4)'-u_p_vals(:,4)');
% plot(tau_f+grav,z_u);
% tau_f = [0 eta_f*diff([diff(u_f_vals(:,5)') 0])./(dz^2)];
% drag_term = drag_mult.*(u_f_vals(:,5)'-u_p_vals(:,5)');
% plot(tau_f+grav,z_u);
% tau_f = [0 eta_f*diff([diff(u_f_vals(:,10)') 0])./(dz^2)];
% drag_term = drag_mult.*(u_f_vals(:,10)'-u_p_vals(:,10)');
% plot(tau_f+grav,z_u);
% xlabel("u_p");
% ylabel("z");
% legend("2s","4s","6s","8s","10s");

% dupdz = diff(u_p_vals,1,2)./dz;
% dufdz = diff(u_f_vals,1,2)./dz;
% 
% subplot(3,1,2);
% plot(eta_f*[0 diff(dufdz(1,:))]./dz,z);
% 
% hold on
% plot(eta_f*[0 diff(dufdz(3,:))]./dz,z);
% plot(eta_f*[0 diff(dufdz(5,:))]./dz,z);
% plot(eta_f*[0 diff(dufdz(7,:))]./dz,z);
% plot(eta_f*[0 diff(dufdz(9,:))]./dz,z);
% plot(eta_f*[0 diff(dufdz(11,:))]./dz,z);
% final_tau_f = eta_f*[0 diff(dufdz(2,:))]./dz;
% xlabel("\tau_f");
% ylabel("z");
% hold off
% legend("0s","2s","4s","6s","8s","10s");
% saveas(gcf,"tau_f.png");
% % 
% I = 2.*d.*dupdz(1,:).*sqrt(rho_p)./sqrt(pp_vals+0.000001);
% % subplot(3,1,3);
% plot([0 diff(pp_vals.*mu_I_fn(I))]./dz,z);
% 
% hold on
% grav_p = rho_p*phi_c*g*sind(theta);
% for j=[2,3,4,5,10]
%     I = 2.*d.*[diff(u_f_vals(:,j)') 0]./dz.*sqrt(rho_p)./sqrt(pp_vals+0.000001);
%     drag_term = drag_mult.*(u_f_vals(:,j)'-u_p_vals(:,j)');
%     tau_p = [0 diff(pp_vals.*mu_I_fn(I))]./dz;
%     plot(tau_p(5:end)+drag_term(5:end)+grav_p,z(5:end));
% end    
% xlabel("\tau_p");
% ylabel("z");
% hold off
% legend("0s","2s","4s","6s","8s","10s");
% saveas(gcf,"tau_p.png");

% subplot(2,1,1);
% plot(u_p_vals(:,1),z);
% 
% hold on
% plot(u_p_vals(:,2),z);
% plot(u_p_vals(:,3),z);
% plot(u_p_vals(:,4),z);
% plot(u_p_vals(:,5),z);
% plot(u_p_vals(:,6),z);
% xlabel("u_p");
% ylabel("z");
% hold off

% for j = 1:6
%     pp_vals(:,j) = p_b'-pe_vals(:,j);
%     phi_vals(:,j) = alpha*pp_vals(:,j);
%     test = diff(pe_vals(:,j))./dz;
%     wp_vals(:,j) = [test' 0]'./beta_fn(phi_vals(:,j));
% end
% plot(pe_vals(:,1),z);
% hold on
% plot(pe_vals(:,2),z);
% plot(pe_vals(:,3),z);
% plot(pe_vals(:,4),z);
% plot(pe_vals(:,5),z);
% plot(pe_vals(:,6),z);
% legend("0s","0.5s","1s","1.5s","2s","2.5s");
% xlabel("Excess Pressure");
% ylabel("z");
% hold off
% 
% subplot(2,2,2)
% plot(pp_vals(:,1),z);
% hold on
% plot(pp_vals(:,2),z);
% plot(pp_vals(:,3),z);
% plot(pp_vals(:,4),z);
% plot(pp_vals(:,5),z);
% plot(pp_vals(:,6),z);
% % legend("0s","10s","20s","30s","40s","50s");
% xlabel("Particle Pressure");
% ylabel("z");
% hold off
% subplot(2,2,3)
% plot(phi_c+phi_vals(:,1),z);
% hold on
% plot(phi_c+phi_vals(:,2),z);
% plot(phi_c+phi_vals(:,3),z);
% plot(phi_c+phi_vals(:,4),z);
% plot(phi_c+phi_vals(:,5),z);
% plot(phi_c+phi_vals(:,6),z);
% % legend("0s","10s","20s","30s","40s","50s");
% xlabel("Volume Fraction");
% ylabel("z");
% hold off
% subplot(2,1,2)
% plot(wp_vals(:,1),z);
% hold on
% plot(wp_vals(:,2),z);
% plot(wp_vals(:,3),z);
% plot(wp_vals(:,4),z);
% plot(wp_vals(:,5),z);
% plot(wp_vals(:,6),z);
% legend("0s","10s","20s","30s","40s","50s");
% xlabel("Paricle z velocity");
% ylabel("z");
% hold off

function beta_val=beta_fn(phihat)
    phi_c=0.6; % Volume fraction
    eta_f = 0.0010016; % Pa s
    d=1.43e-4; % grain diameter (m)
    beta_val = 150*(phi_c + phihat).^2.*eta_f./((1-phi_c-phihat).^3.*d^2);
end

function mu_val = mu_I_fn(I)
    mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
    mu2_I=0.557; % 
    I_0 = 0.069;
    mu_val = tanh(1e8*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
end