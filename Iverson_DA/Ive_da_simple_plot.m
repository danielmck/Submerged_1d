function Ive_da_simple_plot
    fname = ["Ive_da_0deg_30init_change_0.txt"];

    sim_num = size(fname,2);
    sim_list = cell(1,sim_num);
    n_times = zeros(sim_num,1);
    custom_times = zeros(sim_num,1);
    
    h0 = 0.0388; % layer height (m)
    d= 1e-4*ones(sim_num,1);
    phi_c= 0.585*ones(sim_num,1);
    eta_f = 0.0010016*ones(sim_num,1);
    g=9.81*ones(sim_num,1); % m/s^2
    theta = 0*ones(sim_num,1); % deg
    alpha = [1e-5]'; % 1/Pa
    rho_f = 1000*ones(sim_num,1);
    density_ratio = 2.5*ones(sim_num,1);

    theta0=13;

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

    cd Results
    for k=1:sim_num
        data_file = load(fname(k));
        n_times(1,k) = size(data_file,1);
        t_vals(k,1:n_times(1,k)) = data_file(:,1);
        h(k,1:n_times(1,k)) = data_file(:,2);
        phi(k,1:n_times(1,k)) = data_file(:,3);
        u(k,1:n_times(1,k)) = data_file(:,4);
        pb(k,1:n_times(1,k)) = data_file(:,5);
    end
    cd ..

    rho = rho_f*(density_ratio*phi_c+(1-phi_c));
    pe = (pb-rho_f*g*cosd(theta)*h);
    pp = rho*h*g*cosd(theta)-pb;
    Iv = 3*u.*eta_f_dl./(h.*(pp));
    tan_psi = phi-phi_c./(1+sqrt(Iv));
    tau_zx = pp.*mu_Iv_fn(Iv)+(1-phi)*eta_f_dl*2.*u./h;
    D = -2*kappa_dl./(eta_f_dl.*h).*pe;

    dhdt = (rho-1)./rho*D;
%         dhdt=0;
    dphidt = -phi.*D./h;
    dudt = sind(theta).*h-tau_zx./rho;
    
    dudt_approx = sind(theta) - cosd(theta).*tand(theta0)+tand(theta0).*pe./(rho-1);
    dudt_approx2 = sind(theta) - pp.*tand(theta0)./(rho-1);
    diffusion = -3.*kappa_dl./(alpha_dl.*eta_f_dl.*h.^2).*pe;
    dilatancy = -3.*u./(h.*alpha_dl).*(tan_psi);
    dpbdt = -3.*kappa_dl./(alpha_dl.*eta_f_dl.*h.^2).*pe+cosd(theta).*dhdt/4-3.*u./(h.*alpha_dl).*(tan_psi);

    % Early timescale
    dIvdt_u = dudt.*3.*eta_f_dl./pp;
    dIvdt_pp = 1./pp.*Iv.*dpbdt;
    dIvdt = dudt.*3.*eta_f_dl./pp+1./pp.*Iv.*dpbdt;
    dIvdt_brack = dudt+u./pp.*dpbdt;
    dIvdt_phi = 2.*phi_c.*(phi-phi_c)./(phi.^3).*dphidt;
    
    init_Iv = newt_solve_crit_Iv(theta0,density_ratio,1);
    u_0 = init_Iv.*(rho-1)*cosd(theta0)/3/eta_f_dl;
    Iv_0 = u_0(1)/(rho*cosd(theta)-cosd(theta0))*3*eta_f_dl;
    tan_psi_0 = (phi-phi_c)+phi*sqrt(Iv_0);
    tan_psi_ss = alpha_dl.*(rho-1).*cosd(theta)/3./u.^2.*(sind(theta)-tand(theta0)/cosd(theta));
    tan_psi_eq = -(sind(theta)-tand(theta0)/cosd(theta)).*3.*eta_f_dl./((rho-1).*cosd(theta0))/(-init_Iv^2./eta_f_dl/((alpha_dl)));
    tan_psi_approx = tan_psi_eq + (tan_psi_0 - tan_psi_eq).*exp(-phi.*Iv_0.^(3/2)./(2.*eta_f_dl.*alpha_dl).*t_vals);
    pp_approx_early = (rho-1).*cosd(theta0).*(phi(1)-phi_c).^2./(tan_psi_approx + (phi_c-phi(1))).^2;
    
    pp_med_ss = sind(theta).*(rho-1)./tand(theta0);
    pp_approx = pp_med_ss+((rho-1).*cosd(theta0)-pp_med_ss).*exp(-3.*eta_f_dl.*phi.^2.*tand(13)./((rho-1).*(phi-phi_c).^2).*t_vals);
    u_approx = init_Iv.*pp_approx./(3.*eta_f_dl);
    
    fast_scale = 2*eta_f_dl.*alpha_dl./(phi.*Iv.^(3/2));
    medium_scale = u./sind(9);
    
%     SetPaperSize(7,7)
    C = viridis(3);
    
    for i=1:sim_num
        t_start = 0;
        t_stop = 500;
        t_begin = max(sum(t_vals(i,1:n_times(1,i))<t_start)-1,1);
        t_end = min(sum(t_vals(i,1:n_times(1,i))<t_stop)+1,n_times(1,i));
        plot(t_vals(i,t_begin:t_end),u(i,t_begin:t_end),'color',C(1,:),'DisplayName',"Depth Averaged Model",'LineWidth',1)
        hold on
    %     plot(t_vals(i,t_begin:t_end),dIvdt_phi(i,t_begin:t_end),"red",'DisplayName',"Exact Value")
    %     plot(t_vals(i,t_begin:t_end),(3.*u(i,t_begin:t_end).*eta_f_dl./Iv_phi(i,t_begin:t_end)),'DisplayName',"Approximation")
%         plot(t_vals(i,t_begin:t_end),(u_approx(i,t_begin:t_end)),'LineStyle','--','color',C(2,:),'DisplayName',"Approximation",'LineWidth',1)
    end
    ylabel('$u$')
    xlabel('$t$')
    xlim([0,500])
    legend("Location","best")
    figname = "Ive_da_5deg_u_medium.pdf";
%     exp_graph(gcf,figname)
%     movefile(figname,"../Iverson_Closure/Figures/SecondYearReport")
end