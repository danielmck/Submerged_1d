% fname = ["Ive_da_13_deep_quad.txt"]; %,"Ive_da_13_deep_quad_alpha4.txt","Ive_da_13_deep_quad_alpha5.txt"];
% fname = ["Ive_da_5_deep_alpha3.txt","Ive_da_5_deep_alpha4.txt","Ive_da_5_deep_alpha5.txt"];
% fname = ["Ive_da_4_deep_9_2_start.txt","Ive_da_4_deep_9_start.txt"];
fname = ["Ive_da_9_deep_13_start.txt"];
plot_titles = ["$\alpha = 10^{-3}$","$\alpha = 10^{-4}$","$\alpha = 10^{-5}$"];

sim_num = size(fname,2);
sim_list = cell(1,sim_num);
n_times = zeros(sim_num,1);
custom_times = zeros(sim_num,1);

N=[200,200,200]';
h = [4e-2,4e-2,4e-2]';
d= 1.43e-5*ones(sim_num,1);

phi_c= 0.585*ones(sim_num,1);
eta_f = 0.0010016*ones(sim_num,1);
g=9.81*ones(sim_num,1); % m/s^2
theta = 25*ones(sim_num,1); % deg
alpha = [1e-4,1e-4]'; % 1/Pa
rho_f = 1000*ones(sim_num,1);
density_ratio = 2.5*ones(sim_num,1);
t_step = 10*ones(sim_num,1);


N = 200; % number of discretisation points in z for each of tau, u
h0 = 4e-2; % layer height (m)
d=1.43e-5; % grain diameter (m)

mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 10^8;

fric_ang = 0.65;

kappa = ((1-phi_c).^3.*d.^2)./(150.*phi_c.^2);

buoyancy = -rho_f.*(density_ratio-1).*g.*phi_c.*cosd(theta);

v_scale = sqrt(g.*h0);
p_scale = rho_f.*g.*h0;
t_scale = sqrt(h0./g);
z_scale = h0;
rho = phi_c.*density_ratio+1-phi_c;

d_dl = d./z_scale;
eta_f_dl = eta_f./(p_scale.*t_scale);
alpha_dl = alpha.*p_scale;
kappa_dl = kappa./(z_scale).^2;

t_vals = zeros(sim_num,5001);
h = zeros(sim_num,5001);
phi = zeros(sim_num,5001);
u = zeros(sim_num,5001);
pe = zeros(sim_num,5001);

cd DA_Results
for k=1:sim_num
    data_file = load(fname(k));
    custom_time = (size(data_file,2) == 5);
    custom_times(1,k) = custom_time;
    n_times(1,k) = size(data_file,1);
    if (custom_time)
        t_vals(k,1:n_times(1,k)) = data_file(:,1);
    else
        t_vals(k,1:n_times(1,k)) = linspace(0,(size(data_file,1)-1)*t_step,size(data_file,1));
    end
    h(k,1:n_times(1,k)) = data_file(:,custom_time+1);
    phi(k,1:n_times(1,k)) = data_file(:,custom_time+2);
    u(k,1:n_times(1,k)) = data_file(:,custom_time+3);
    pe(k,1:n_times(1,k)) = data_file(:,custom_time+4)'-cosd(theta(k)).*h(k,1:n_times(1,k));
end
cd ..

rho = density_ratio.*phi+(1-phi);
pp = (rho-1).*h.*cosd(theta)-pe;

Iv = 3.*u.*eta_f_dl./(h.*(pp));
Iv_base = 3.*u.*eta_f_dl./(h.*(pp));
Iv_phi = (phi_c-phi).^2./(phi.^2);

tan_psi = phi-phi_c./(1+sqrt(Iv));
tau_zx = pp.*mu_Iv_fn(Iv_base)+(1-phi).*eta_f_dl.*2.*u./h;
D = -2.*kappa_dl./(eta_f_dl.*h).*pe;

dhdt = (rho-1)./rho.*D;

dphidt = -phi.*D./h;
dudt = sind(theta).*h-tau_zx./rho;
dudt_approx = sind(theta) - cosd(theta).*tand(13)+tand(13).*pe./(rho-1);
dudt_approx2 = sind(theta) - pp.*tand(13)./(rho-1);
diffusion = -3.*kappa_dl./(alpha_dl.*eta_f_dl.*h.^2).*pe;
dilatancy = -3.*u./(h.*alpha_dl).*(tan_psi);
dpbdt = -3.*kappa_dl./(alpha_dl.*eta_f_dl.*h.^2).*pe+cosd(theta).*dhdt/4-3.*u./(h.*alpha_dl).*(tan_psi);

dIvdt_u = dudt.*3.*eta_f_dl./pp;
dIvdt_pp = 1./pp.*Iv.*dpbdt;
dIvdt = dudt.*3.*eta_f_dl./pp+1./pp.*Iv.*dpbdt;
dIvdt_brack = dudt+u./pp.*dpbdt;
dIvdt_phi = 2.*phi_c.*(phi-phi_c)./(phi.^3).*dphidt;

dIvdt_fast = -3.*eta_f_dl.*tand(13)./(rho-1)+sind(theta)./u(:,1).*Iv+(phi_c-phi(:,1))./(alpha_dl.*eta_f_dl).*Iv.^2-phi(:,1)./(alpha_dl.*eta_f_dl).*Iv.^(5/2);

pp_approx = sind(theta).*(rho-1)./tand(13)+((rho-1).*cosd(13)-sind(theta).*(rho-1)./tand(13)).*exp(-3.*eta_f_dl.*phi.^2.*tand(13)./((rho-1).*(phi-phi_c).^2).*t_vals);

tan_psi_approx = alpha_dl.*eta_f_dl./u./((phi-phi_c)./phi).^2.*(sind(theta)-pp.*tand(13)./(rho-1));
dpbdt_approx = -3*eta_f_dl./((phi-phi_c)./phi).^2.*(sind(theta)-pp.*tand(13)./(rho-1));
dpbdt_approx2 = -3*eta_f_dl./((phi-phi_c)./phi).^2.*(dudt);
dpbdt_approx3 = -pp./Iv.*(dIvdt_u-dIvdt_phi);
tan_psi_ss = alpha_dl.* eta_f_dl./(u(:,1).*Iv).*dudt;
Iv_0 = Iv_phi.*((rho-1).*cosd(13)./(rho.*cosd(theta)-cosd(13)));
tan_psi_0 = (phi-phi_c).*(1-sqrt((rho-1).*cosd(13)./(rho.*cosd(theta)-cosd(13))));
tan_psi_approx = tan_psi_ss + (tan_psi_0 - tan_psi_ss).*exp(-phi.*Iv_0.^(3/2)./(2.*eta_f_dl.*alpha_dl).*t_vals);
pp_approx_early = (rho-1).*cosd(13).*(phi(1)-phi_c).^2./(tan_psi_approx + (phi_c-phi(1))).^2;

Iv_comp = 3.*u.*eta_f_dl./(rho.*cosd(theta)-cosd(13));

fast_scale = eta_f_dl.*alpha_dl./(phi.*Iv.^(3/2));
medium_scale = u./sind(9);
%% 

% SetPaperSize(10,10);
% figure
hold on
for i=1:sim_num
    t_start = 0;
    t_stop = 500000;
    t_begin = max(sum(t_vals(i,1:n_times(1,i))<t_start)-1,1);
    t_end = min(sum(t_vals(i,1:n_times(1,i))<t_stop)+1,n_times(1,i));
%     plot(t_vals(i,t_begin:t_end),(phi(i,t_begin:t_end)-phi_c(1)),'DisplayName',"Exact Value")
    plot(t_vals(i,t_begin:t_end),(u(i,t_begin:t_end)),'DisplayName',"Approximation")
    
%     plot(t_vals(i,t_begin:t_end),(u_p(i,t_begin:t_end)),'DisplayName',"Exponential Approximtaion")
%     plot(t_vals(i,t_begin:t_end),(dIvdt_pp(i,t_begin:t_end)),'DisplayName',"Approximtaion")
%     plot(t_vals(i,t_begin:t_end),(tand(theta(i))-tand(13))./(rho(i,t_begin:t_end)-1),'DisplayName',"Approximtaion")
%     plot(linspace(0,0.06,1000),mu_Iv_fn(linspace(0,0.06,1000)))
end
% plot(t_vals(i,t_begin:t_end),(rho(2,1)-1)*cosd(13).*ones(1,t_end-t_begin+1),"--k",'DisplayName',"Previous Slope Value")
legend('Location', "best",'UserData', 8);
xlabel("$t$");
ylabel('$p_p$');%,'Position',[-35 0.0045]);
% xlim([t_start t_stop])
% ylim([0 0.01])
box on
% title('Exact Value and Approximation of $p_p$ during the Short Timescale')

ax = gca;
%      set(gca,'ytick',[])
% ax.YAxis.Exponent = 0;
% yticks(-1e-4:2e-5:0)
% ytickformat('%5.0e');
% ax.XAxis.Exponent = 0;

% PrintFig('Ive_da_fast_pp')
% plot(linspace(0,1000,101),dudt_approx(1:101))
% plot(linspace(0,5000,5001),tand(theta)*rho./())

function mu_val = mu_Iv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    reg_param = 10^8;
    phi_c = 0.585;
    mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
end