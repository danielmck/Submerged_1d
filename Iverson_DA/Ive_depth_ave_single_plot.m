fname = ["Ive_da_5deg_Fr_10_mud.txt"]; % 
% plot_titles = ["$\alpha = 10^{-3}$","$\alpha = 10^{-4}$","$\alpha = 10^{-5}$"];

sim_num = size(fname,2);
sim_list = cell(1,sim_num);
n_times = zeros(sim_num,1);
custom_times = zeros(sim_num,1);

h0 = [3e-1]';
d= 1e-4*ones(sim_num,1);

phi_c= 0.585*ones(sim_num,1);
eta_f = 0.0010016*ones(sim_num,1);
g=9.81*ones(sim_num,1); % m/s^2
theta = 5*ones(sim_num,1); % deg
theta0 = 10*ones(sim_num,1); % deg
alpha = [1e-5]'; % 1/Pa
rho_f = 1000*ones(sim_num,1);
density_ratio = 2.5*ones(sim_num,1);
% t_step = 10*ones(sim_num,1);

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
beta_dl = eta_f_dl/kappa_dl;


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
dudt_approx = sind(theta) - cosd(theta).*tand(theta0)+tand(theta0).*pe./(rho-1);
dudt_approx2 = sind(theta) - pp.*tand(theta0)./(rho-1);
diffusion = -3.*kappa_dl./(alpha_dl.*eta_f_dl.*h.^2).*pe;
dilatancy = -3.*abs(u)./(h.*alpha_dl).*(tan_psi);
dpbdt = diffusion+cosd(theta).*dhdt/4+dilatancy;

Iv_init = Iv(1);
tan_psi_init = tan_psi(1);
phi_init = phi(1);
Iv_eq_init = Iv_phi(1);
u_init = u(1);
rho_init = rho(1);
pp_fast_ss = 3*eta_f_dl*u_init/Iv_eq_init;
pb_fast_ss = rho_init*cosd(theta)-pp_fast_ss;
tan_psi_ss_new = alpha_dl.*eta_f_dl/Iv_eq_init/u_init.*dudt+2/3/u_init/beta_dl*(pb_fast_ss-cosd(theta));
tan_psi_approx_new = tan_psi_ss_new + (tan_psi_init-tan_psi_ss_new).*exp(-phi_init*Iv_eq_init^(3/2)/2/eta_f_dl/alpha_dl*t_vals);
pp_fast_new = 3*u_init*eta_f_dl*phi_init^2./(tan_psi_approx_new+(phi_c-phi_init)).^2;
fast_ts = 2*eta_f_dl*alpha_dl/(phi_init*Iv_eq_init^(3/2));

pp_med_init = 3*u_init*eta_f_dl/Iv_eq_init;
pp_med_ss = sind(theta).*rho./mu_Iv_fn(Iv_phi);
pe_med_ss = (rho_init-1)*cosd(theta)-pp_med_ss;
u_med_ss = pp_med_ss.*Iv_phi./3/eta_f_dl;
pp_approx_new = pp_med_ss + (pp_med_init-pp_med_ss).*exp(-3*eta_f_dl*mu_Iv_fn(Iv_eq_init)/Iv_eq_init/rho_init*t_vals);
pf_approx_new = rho_init*cosd(theta)-pp_approx_new;
u_approx_new = pp_approx_new*Iv_eq_init/3/eta_f_dl;
medium_ts = Iv_eq_init*rho_init/(3*eta_f_dl*mu_Iv_fn(Iv_eq_init));

dphidt_approx = 2.*kappa_dl./(eta_f_dl.*h).*phi./h.*pe_med_ss;
phi_approx = zeros(1,n_times(k));
phi_app=phi_init;
phi_approx(1) = phi_app;
for j = 2:n_times(k)
    if (t_vals(j)>8)
        phi_app = phi_app+(dphidt_approx(j)+dphidt_approx(j-1))*(t_vals(j)-t_vals(j-1))/2;
    end    
        phi_approx(j) = phi_app;
   if (u(j)<1.4e-2)
       phi_approx(j) = phi(j);
   end
end
t_depo = t_vals(sum(phi_approx<phi_c));


slow_ts = (phi_c-phi_init)./phi_init./(2.*kappa_dl).*(eta_f_dl).*pe_med_ss;

dIvdt_u = dudt.*3.*eta_f_dl./pp;
dIvdt_pp = 1./pp.*Iv.*dpbdt;
dIvdt = dudt.*3.*eta_f_dl./pp+1./pp.*Iv.*dpbdt;
dIvdt_brack = dudt+u./pp.*dpbdt;
dIvdt_phi = 2.*phi_c.*(phi-phi_c)./(phi.^3).*dphidt;

dIvdt_fast = -3.*eta_f_dl.*tand(13)./(rho-1)+sind(theta)./u(:,1).*Iv+(phi_c-phi(:,1))./(alpha_dl.*eta_f_dl).*Iv.^2-phi(:,1)./(alpha_dl.*eta_f_dl).*Iv.^(5/2);

pp_approx = sind(theta).*(rho-1)./tand(theta0)+((rho-1).*cosd(theta0)-sind(theta).*(rho-1)./tand(theta0)).*exp(-3.*eta_f_dl.*phi.^2.*tand(theta0)./((rho-1).*(phi-phi_c).^2).*t_vals);
u_approx = Iv_phi.*pp_approx./(3.*eta_f_dl);

tan_psi_approx = alpha_dl.*eta_f_dl./u./((phi-phi_c)./phi).^2.*(sind(theta)-pp.*tand(theta0)./(rho-1));
dpbdt_approx = -3*eta_f_dl./((phi-phi_c)./phi).^2.*(sind(theta)-pp.*tand(13)./(rho-1));
dpbdt_approx2 = -3*eta_f_dl./((phi-phi_c)./phi).^2.*(dudt);
dpbdt_approx3 = -pp./Iv.*(dIvdt_u-dIvdt_phi);
tan_psi_ss = alpha_dl.*pp/3./u.^2.*dudt;
Iv_0 = Iv_phi.*(cosd(theta0)./(cosd(theta)));
tan_psi_0 = (phi-phi_c).*(1-sqrt(cosd(theta0)./(cosd(theta))));
tan_psi_approx = tan_psi_ss + (tan_psi_0 - tan_psi_ss).*exp(-phi.*Iv_0.^(3/2)./(2.*eta_f_dl.*alpha_dl).*t_vals);
pp_approx_early = (rho-1).*cosd(13).*(phi(1)-phi_c).^2./(tan_psi_approx + (phi_c-phi(1))).^2;

Iv_comp = 3.*u.*eta_f_dl./(rho.*cosd(theta)-cosd(13));

tan_psi_ss2 = alpha_dl.* eta_f_dl./(u(:,1).*Iv).*dudt-kappa_dl./(eta_f_dl.*u).*pe;
pp_eq = (3.*eta_f_dl.*u)./Iv_phi;
pp_exact = (3.*eta_f_dl.*u).*(phi-tan_psi).^2./((phi-phi_c)+tan_psi).^2;
pp_prime = ((alpha_dl.*pp_eq)./(3.*u.^2).*(sind(theta)-pp_eq.*mu_Iv_fn(Iv_phi)./(rho-1))-kappa_dl./(u.*eta_f_dl).*((rho-1).*cosd(theta)-pp_eq))./(Iv_phi./pp_eq-(alpha_dl)./(3.*u.^2).*(sind(theta)-2.*pp_eq.*mu_Iv_fn(Iv_phi)./(rho-1))-kappa_dl./(u.*eta_f_dl));

fast_scale = eta_f_dl.*alpha_dl./(phi.*Iv.^(3/2));
medium_scale = u./sind(theta0);

Iv_phi = (phi_c-phi).^2./phi.^2;
pp_phi = rho.*h.*sind(theta)./mu_Iv_fn(Iv_phi);
pf_phi = rho*cosd(theta).*h-pp_phi;

D_phi = -2.*kappa_dl./(eta_f_dl.*h).*(pf_phi-cosd(theta)*h);
dphidt_phi = -phi.*D_phi./h;
%%
scale_vec = [t_scale,z_scale,1,v_scale,p_scale];
scale_data = data_file.*scale_vec;
dl_name = fname(k);
dl_name = split(dl_name,".txt");
dl_name = dl_name(1);
scale_name = dl_name+"_dim.txt";
save("Results/"+scale_name,"scale_data","-ascii")

%% 

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
    SetPaperSize(7.8,7.8);  
    if (i == 1)
        colour = "blue";
    else
        colour = "red";
    end
    
    
    t_start = 0;
    t_stop = 8;
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
    
    [pos,ind] = max(phi_approx);
    
    C = plasma(5);
    plot(t_vals(i,t_begin:t_end), Iv(t_begin:t_end),'DisplayName',"Model", 'color', '#fc8d62','LineWidth',1.3) %[0.127,0.201,0.127])
%     plot(t_vals(i,t_begin:t_end), phi_approx(t_begin:t_end),'DisplayName',"Approximation", 'color', '#8da0cb','LineWidth',1.3)
%     plot(t_vals(i,t_begin:t_end), tan_psi_approx(t_begin:t_end),'DisplayName',"Smoothed $\phi$")

end
% plot(t_vals(i,t_begin:t_end),(rho(2,1)-1)*cosd(13).*ones(1,t_end-t_begin+1),"--k",'DisplayName',"Previous Slope Value")
legend('Location', "best",'UserData', 8);
xlabel("$t$");
ylabel('$\phi$');%,'Position',[-35 0.0045]);
xlim([t_start t_stop])
% ylim([0 0.01])
box on
% title('Exact Value and Approximation of $\tan \psi$ during the Fast Timescale')

ax = gca;
%      set(gca,'ytick',[])
% ax.YAxis.Exponent = 0;
% yticks(-1e-4:2e-5:0)
% ytickformat('%5.0e');
% ax.YAxis.Exponent = 0;
% exp_graph(gcf,"Ive_da_phi_mud_slow.pdf")