% fname = ["Ive_da_13_acc.txt"]; %,"Ive_da_13_deep_quad_alpha4.txt","Ive_da_13_deep_quad_alpha5.txt"];
% fname = ["Ive_da_5_deep_alpha3.txt","Ive_da_5_deep_alpha4.txt","Ive_da_5_deep_alpha5.txt"];
% fname = ["Ive_da_4_deep_9_2_start.txt"]; %,"Ive_da_4_deep_9_start.txt"];
fname = ["sin_u_period.txt"];
% plot_titles = ["$\alpha = 10^{-3}$","$\alpha = 10^{-4}$","$\alpha = 10^{-5}$"];

sim_num = size(fname,2);
sim_list = cell(1,sim_num);
n_times = zeros(sim_num,1);
custom_times = zeros(sim_num,1);

h0 = [1e-1]';
d= 1e-4*ones(sim_num,1);

phi_c= 0.585*ones(sim_num,1);
eta_f = 0.0010016*ones(sim_num,1);
g=9.81*ones(sim_num,1); % m/s^2
theta = 0*ones(sim_num,1); % deg
alpha = [1e-5]'; % 1/Pa
rho_f = 1000*ones(sim_num,1);
density_ratio = 2.5*ones(sim_num,1);
t_step = 10*ones(sim_num,1);

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
pb = zeros(sim_num,5001);
pe = zeros(sim_num,5001);

cd Results
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
    phi(k,1:n_times(1,k)) = data_file(:,custom_time+2)./h(k,1:n_times(1,k))';
    u(k,1:n_times(1,k)) = data_file(:,custom_time+3);
    pb(k,1:n_times(1,k)) = data_file(:,custom_time+4);
    pe(k,1:n_times(1,k)) = pb(k,1:n_times(1,k))-cosd(theta(k)).*h(k,1:n_times(1,k));
end
cd ..

rho = density_ratio.*phi+(1-phi);
pp = (rho-1).*h.*cosd(theta)-pe;

Iv = 3.*abs(u).*eta_f_dl./(h.*(pp));
Iv_base = 3.*abs(u).*eta_f_dl./(h.*(pp));
Iv_phi = (phi_c-phi).^2./(phi.^2);

phi_Iv = phi_c./(1+sqrt(Iv));
tan_psi = phi-phi_Iv;
tau_zx = pp.*mu_Iv_fn(Iv_base)+(1-phi).*eta_f_dl.*2.*u./h;
D = -2.*kappa_dl./(eta_f_dl.*h).*pe;

dhdt = (rho-1)./rho.*D;

dphidt = -phi.*D./h;
dhphidt = -phi.*D./rho;
dudt = sind(theta).*h-tau_zx./rho;
dudt_approx = sind(theta) - cosd(theta).*tand(13)+tand(13).*pe./(rho-1);
dudt_approx2 = sind(theta) - pp.*tand(13)./(rho-1);
diffusion = -3.*kappa_dl./(alpha_dl.*eta_f_dl.*h.^2).*pe;
dilatancy = -3.*abs(u)./(h.*alpha_dl).*(tan_psi);
dpbdt = diffusion+cosd(theta).*dhdt/4+dilatancy;

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
tan_psi_ss = alpha_dl.*pp/3./u.^2.*dudt;
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
scale_vec = [t_scale,z_scale,1,v_scale,p_scale];
scale_data = data_file.*scale_vec;
dl_name = fname(k);
dl_name = split(dl_name,".txt");
dl_name = dl_name(1);
scale_name = dl_name+"_dim.txt";
save("Results/"+scale_name,"scale_data","-ascii")

%% 

% SetPaperSize(10,10);
% figure
% hold on
% red_vals = load("../Iverson_DA/DA_Results/Ive_da_red_4_deep_9_2_start.txt");
% red_times = red_vals(:,1);
% red_h = red_vals(:,2);
% red_phi = red_vals(:,3);
% red_u = red_vals(:,4);
% % red_pe = red_vals(:,5)-cosd(theta);
% red_pp = (red_phi.^2./(red_phi-phi_c).^2).*3.*red_u.*eta_f_dl./red_h;
% 
% theta_vals = load("../Iverson_DA/Results/Ive_da_theta_4_deep_9_2_start.txt");
% 
% theta_start = 9.2;
% theta_stop = 4;
% theta_delta_t = 10;
% theta_times = theta_vals(:,1);
% theta_change = max(theta_start-theta_times.*(theta_start-theta_stop)/theta_delta_t,theta_stop);
% theta_h = theta_vals(:,2);
% theta_phi = theta_vals(:,3);
% theta_u = theta_vals(:,4);
% theta_pe = theta_vals(:,5)-cosd(theta_change).*theta_h;
% theta_pp = (density_ratio*theta_phi+(1-theta_phi)-1).*theta_h.*cosd(theta_change)-theta_pe;
% theta_Iv = 3.*theta_u.*eta_f_dl./(theta_h.*(theta_pp));
% theta_tan_psi = theta_phi-phi_c./(1+sqrt(theta_Iv));

smooth_span = 402.5;

Iv_smooth = zeros(sim_num,5001);
phi_smooth = zeros(sim_num,5001);
phi_Iv_smooth = zeros(sim_num,5001);
pe_smooth = zeros(sim_num,5001);
for i=1:sim_num
    for k=1:n_times(1,i)
        t_now = t_vals(i,k);
        n_below = min(sum((t_vals(i,:)-t_now > -smooth_span/2) & (t_vals(i,:)<t_now))+1,k-1);
        n_above = min(sum(t_vals(i,:)-t_now < smooth_span/2 & (t_vals(i,:)>t_now))+1,n_times(1,i)-k);
        pe_ave = 0;
        Iv_ave = 0;
        phi_Iv_ave = 0;
        phi_ave = 0;
        for l=-n_below+1:n_above
            pe_ave = pe_ave+(phi(i,k+l)+phi(i,k+l-1))/2*(t_vals(i,k+l)-t_vals(i,k+l-1));
            Iv_ave = Iv_ave+(Iv(i,k+l)+Iv(i,k+l-1))/2*(t_vals(i,k+l)-t_vals(i,k+l-1));
            phi_ave = phi_ave+(phi(i,k+l)+phi(i,k+l-1))/2*(t_vals(i,k+l)-t_vals(i,k+l-1));
            phi_Iv_ave = phi_Iv_ave+(phi_Iv(i,k+l)+phi_Iv(i,k+l-1))/2*(t_vals(i,k+l)-t_vals(i,k+l-1));
        end
        pe_ave = pe_ave/(t_vals(i,k+n_above)-t_vals(i,k-n_below));
        Iv_ave = Iv_ave/(t_vals(i,k+n_above)-t_vals(i,k-n_below));
        phi_ave = phi_ave/(t_vals(i,k+n_above)-t_vals(i,k-n_below));
        phi_Iv_ave = phi_Iv_ave/(t_vals(i,k+n_above)-t_vals(i,k-n_below));
        pe_smooth(i,k) = pe_ave;
        Iv_smooth(i,k) = Iv_ave;
        phi_smooth(i,k) = phi_ave;
        phi_Iv_smooth(i,k) = phi_Iv_ave;
    end
%     SetPaperSize(8,8);  
    if (i == 1)
        colour = "blue";
    else
        colour = "red";
    end
    t_start = 0;
    t_stop = 650;
    t_begin = max(sum(t_vals(i,1:n_times(1,i))<t_start)-1,1);
    t_end = min(sum(t_vals(i,1:n_times(1,i))<t_stop)+1,n_times(1,i));
    pe_ave = 0;
    phi_pe_ave = 0;
    pe_mag = 0;
    phi_pe_mag = 0;
    for k=2:t_end
        pe_ave = pe_ave + (pe(k)+pe(k-1))/2*(t_vals(k)-t_vals(k-1));
        phi_pe_ave = phi_pe_ave + (pe(k)*phi(k)+pe(k-1)*phi(k-1))/2*(t_vals(k)-t_vals(k-1));
        pe_mag = pe_mag + (abs(pe(k))+abs(pe(k-1)))/2*(t_vals(k)-t_vals(k-1));
        phi_pe_mag = phi_pe_mag + (abs(pe(k).*phi(k))+abs(pe(k-1).*phi(k-1)))/2*(t_vals(k)-t_vals(k-1));
    end
    pe_ave = pe_ave/(t_vals(t_end)-t_vals(t_begin));
    phi_pe_ave = phi_pe_ave/(t_vals(t_end)-t_vals(t_begin));
    pe_mag = pe_mag/(t_vals(t_end)-t_vals(t_begin));
    phi_pe_mag = phi_pe_mag/(t_vals(t_end)-t_vals(t_begin));
    hold on
    plot(t_vals(i,t_begin:t_end), phi(t_begin:t_end),'DisplayName',"Smoothed $\phi$")
%     plot(t_vals(i,t_begin:t_end),phi_Iv_smooth(i,t_begin:t_end),'DisplayName',"Smoothed $\phi(I_v)$")
%     plot(t_vals(i,t_begin:t_end),phi_c./(1+sqrt(Iv_smooth(i,t_begin:t_end))),'DisplayName',"Exact Value")
    
%     plot(t_vals(i,t_begin:t_end),dIvdt_phi(i,t_begin:t_end),"red",'DisplayName',"Exact Value")
%     plot(t_vals(i,t_begin:t_end),(3.*u(i,t_begin:t_end).*eta_f_dl./Iv_phi(i,t_begin:t_end)),'DisplayName',"Approximation")
%     plot(t_vals(i,t_begin:t_end),(pp_approx_early(i,t_begin:t_end)),"red",'DisplayName',"Approximation")
%     plot(theta_times,theta_tan_psi)


%     plot(t_vals(i,t_begin:t_end),(u_p(i,t_begin:t_end)),'DisplayName',"Exponential Approximtaion")
%     plot(t_vals(i,t_begin:t_end),(dIvdt_pp(i,t_begin:t_end)),'DisplayName',"Approximtaion")
%     plot(t_vals(i,t_begin:t_end),(tand(theta(i))-tand(13))./(rho(i,t_begin:t_end)-1),'DisplayName',"Approximtaion")
%     plot(linspace(0,0.06,1000),mu_Iv_fn(linspace(0,0.06,1000)))
end
% plot(t_vals(i,t_begin:t_end),(rho(2,1)-1)*cosd(13).*ones(1,t_end-t_begin+1),"--k",'DisplayName',"Previous Slope Value")
legend('Location', "best",'UserData', 8);
xlabel("$t$");
ylabel('$p_e$ (Pa)');%,'Position',[-35 0.0045]);
xlim([t_start t_stop])
% ylim([0 0.01])
box on
% title('Exact Value and Approximation of $\tan \psi$ during the Fast Timescale')
% exp_graph(gcf,"u_sinus_phi_Iv_smooth.pdf")
ax = gca;
%      set(gca,'ytick',[])
% ax.YAxis.Exponent = 0;
% yticks(-1e-4:2e-5:0)
% ytickformat('%5.0e');
ax.YAxis.Exponent = 0;
