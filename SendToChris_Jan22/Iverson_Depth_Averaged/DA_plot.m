fname = ["Ive_da_5_deep_alpha4.txt"];
plot_titles = ["$\alpha = 10^{-3}$","$\alpha = 10^{-4}$","$\alpha = 10^{-5}$"];

sim_num = size(fname,2);
sim_list = cell(1,sim_num);
n_times = zeros(sim_num,1);
custom_times = zeros(sim_num,1);

sim_type = strings(1,sim_num);
% Need to specify the type of simulation:
% dilatancy - "dil"
% diffusion only - "pdriv"
% constant pressure profile - "pcon"
% constant velocity profile - "ucon"

% Need to create lists storing the parameters in positions that corrrelate to
% the data in sim_list
d=zeros(1,sim_num);
dz = zeros(1,sim_num);
phi_c=zeros(1,sim_num);
eta_f_dl = zeros(1,sim_num);
g=9.81*ones(1,sim_num); % m/s^2
theta = zeros(1,sim_num); % deg
alpha_dl = zeros(1,sim_num); % 1/Pa
density_ratio = zeros(1,sim_num);

record = readtable('Results/result_record.csv');

% Extracts data from result_record.csv
for k=1:sim_num   
    in_table = strcmp(record.Name, names(k));
    cd Results
    sim_list{k,1} = load(names(k));
    custom_times(1,k) = (size(data_file,2) == 5);
    n_times(1,k) = size(sim_list{k,1},1);
    cd ../
    sim_type(k) = record.sim_type(in_table);
    h0(k) = record.h0(in_table);
    d(k) = record.d(in_table);
    phi_c(k) = record.phi_c(in_table);
    density_ratio(k) = record.rho_r(in_table);
    eta_f_dl(k) = record.eta_f(in_table);
    theta(k) = record.theta(in_table);
    alpha_dl(k) = record.alpha(in_table);
end

d=1.43e-5; % grain diameter (m)

mu1_Iv = 0.32;
mu2_Iv = 0.7;
Iv_0 = 0.005;

reg_param = 10^8;

fric_ang = 0.65;

kappa = ((1-phi_c).^3.*d.^2)./(150.*phi_c.^2); % Matrix Permeability

buoyancy = -rho_f.*(density_ratio-1).*g.*phi_c.*cosd(theta);

v_scale = sqrt(g.*h0); % Typical Scales
t_scale = sqrt(h0./g);
z_scale = h0;
rho = phi_c.*density_ratio+1-phi_c;

d_dl = d./z_scale; % Dimensionless Quantities
kappa_dl = kappa./(z_scale).^2;

large_sim_size = max(n_times);

t_vals = zeros(sim_num,large_sim_size);
h = zeros(sim_num,large_sim_size);
phi = zeros(sim_num,large_sim_size);
u = zeros(sim_num,large_sim_size);
pe = zeros(sim_num,large_sim_size);

for k=1:sim_num
    data_file = sim_list{k,1};
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

% Reconstructs all of the quantities used in the simulation
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
u_approx = Iv_phi.*pp_approx./(3.*eta_f_dl);

tan_psi_approx = alpha_dl.*eta_f_dl./u./((phi-phi_c)./phi).^2.*(sind(theta)-pp.*tand(13)./(rho-1));
dpbdt_approx = -3*eta_f_dl./((phi-phi_c)./phi).^2.*(sind(theta)-pp.*tand(13)./(rho-1));
dpbdt_approx2 = -3*eta_f_dl./((phi-phi_c)./phi).^2.*(dudt);
dpbdt_approx3 = -pp./Iv.*(dIvdt_u-dIvdt_phi);
tan_psi_ss = alpha_dl*pp/3./u.^2.*dudt;
Iv_0 = Iv_phi.*((rho-1).*cosd(13)./(rho.*cosd(theta)-cosd(13)));
tan_psi_0 = (phi-phi_c).*(1-sqrt((rho-1).*cosd(13)./(rho.*cosd(theta)-cosd(13))));
tan_psi_approx = tan_psi_ss + (tan_psi_0 - tan_psi_ss).*exp(-phi.*Iv_0.^(3/2)./(2.*eta_f_dl.*alpha_dl).*t_vals);
pp_approx_early = (rho-1).*cosd(13).*(phi(1)-phi_c).^2./(tan_psi_approx + (phi_c-phi(1))).^2;

Iv_comp = 3.*u.*eta_f_dl./(rho.*cosd(theta)-cosd(13));

tan_psi_ss2 = alpha_dl.* eta_f_dl./(u(:,1).*Iv).*dudt-kappa_dl./(eta_f_dl.*u).*pe;
pp_eq = (3.*eta_f_dl.*u)./Iv_phi;
pp_exact = (3.*eta_f_dl.*u).*(phi-tan_psi).^2./((phi-phi_c)+tan_psi).^2;
pp_prime = ((alpha_dl.*pp_eq)./(3.*u.^2).*(sind(theta)-pp_eq.*mu_Iv_fn(Iv_phi)./(rho-1))-kappa_dl./(u.*eta_f_dl).*((rho-1).*cosd(theta)-pp_eq))./(Iv_phi./pp_eq-(alpha_dl)./(3.*u.^2).*(sind(theta)-2.*pp_eq.*mu_Iv_fn(Iv_phi)./(rho-1))-kappa_dl./(u.*eta_f_dl));

fast_scale = eta_f_dl.*alpha_dl./(phi.*Iv.^(3/2));
medium_scale = u./sind(9);
%% 

% SetPaperSize(10,10);
% figure

% Loads data from a reduced simulation
hold on
red_vals = load("../Iverson_DA/DA_Results/Ive_da_red_4_deep_9_2_start.txt");
red_times = red_vals(:,1);
red_h = red_vals(:,2);
red_phi = red_vals(:,3);
red_u = red_vals(:,4);
% red_pe = red_vals(:,5)-cosd(theta);
red_pp = (red_phi.^2./(red_phi-phi_c).^2).*3.*red_u.*eta_f_dl./red_h;

% Loads data from a simulation with a variable theta
theta_vals = load("../Iverson_DA/DA_Results/Ive_da_theta_4_deep_9_2_start.txt");
theta_start = 9.2;
theta_stop = 4;
theta_delta_t = 10;
theta_times = theta_vals(:,1);
theta_change = max(theta_start-theta_times.*(theta_start-theta_stop)/theta_delta_t,theta_stop);
theta_h = theta_vals(:,2);
theta_phi = theta_vals(:,3);
theta_u = theta_vals(:,4);
theta_pe = theta_vals(:,5)-cosd(theta_change).*theta_h;
theta_pp = (density_ratio*theta_phi+(1-theta_phi)-1).*theta_h.*cosd(theta_change)-theta_pe;
theta_Iv = 3.*theta_u.*eta_f_dl./(theta_h.*(theta_pp));
theta_tan_psi = theta_phi-phi_c./(1+sqrt(theta_Iv));

for i=1:sim_num
    t_start = 500;
    t_stop = 1500;
    t_begin = max(sum(t_vals(i,1:n_times(1,i))<t_start)-1,1);
    t_end = min(sum(t_vals(i,1:n_times(1,i))<t_stop)+1,n_times(1,i));
    plot(t_vals(i,t_begin:t_end),dpbdt(i,t_begin:t_end),"blue",'DisplayName',plot_titles(i))
%     plot(theta_times,theta_tan_psi)

end
% plot(t_vals(i,t_begin:t_end),(rho(2,1)-1)*cosd(13).*ones(1,t_end-t_begin+1),"--k",'DisplayName',"Previous Slope Value")
legend('Location', "best",'UserData', 8);
xlabel("$t$");
ylabel('$\tan \psi$');
xlim([t_start t_stop])
% ylim([0 0.01])
box on
title('Exact Value and Approximation of $\tan \psi$ during the Fast Timescale')

ax = gca;
ax.YAxis.Exponent = 0;

% PrintFig('Ive_da_early_tan_psi')

function mu_val = mu_Iv_fn(Iv)
    mu1_Iv = 0.32;
    mu2_Iv = 0.7;
    Iv_0 = 0.005;
    reg_param = 10^8;
    phi_c = 0.585;
    mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
end