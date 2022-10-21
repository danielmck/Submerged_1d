function simple_plotting
    names = ["Ive_5deg_13init_deposit.txt"]; %"Ive_comp_4_deep_1_flux.txt"];%,"Ive_comp_4_deep_9_2_start.txt""Ive_comp_4_deep_9_1_start.txt","
    % names = ["Ive_comp_5_deep_custom_time_IC.txt"]; %,"Ive_comp_5_deep_custom_time_IC.txt"];%
    % Loads the data to compare and puts them in a list of simulations

    sim_num = size(names,2);
    sim_list = cell(1,sim_num);
    n_times = zeros(sim_num,1);
    % for j=1:sim_num
    %     n_times(1,j) = size(sim_list{1,j},1);
    % end
    sim_type = strings(1,sim_num);
    % Need to specify the type of simulation:
    % dilatancy - "dil"
    % diffusion only - "pdriv"
    % constant pressure profile - "pcon"
    % constant velocity profile - "ucon"
    % Need to create lists storing the parameters in positions that corrrelate to
    % the data in sim_list

    N=zeros(1,sim_num);
    h = zeros(1,sim_num);
    d=zeros(1,sim_num);
    dz = zeros(1,sim_num);
    phi_c=zeros(1,sim_num);
    eta_f_dl = zeros(1,sim_num);
    g=9.81*ones(1,sim_num); % m/s^2
    theta = zeros(1,sim_num); % deg
    alpha_dl = zeros(1,sim_num); % 1/Pa
    s_frac = [0.6,0.6,0.6];
    density_ratio = zeros(1,sim_num);
    t_step = zeros(1,sim_num);

    record = readtable('Results/result_record.csv');

    for k=1:sim_num
        in_table = strcmp(record.Name, names(k));
        cd Results
        sim_list{k,1} = load(names(k));
        cd ../
        sim_type(k) = record.sim_type(in_table);
        N(k) = record.N(in_table);
        h(k) = record.h(in_table);
        d(k) = record.d(in_table);
        dz(k) = h(k)/(N(k)-0.5);
        phi_c(k) = record.phi_c(in_table);
        density_ratio(k) = record.rho_r(in_table);
        eta_f_dl(k) = record.eta_f(in_table);
        theta(k) = record.theta(in_table);
        alpha_dl(k) = record.alpha(in_table);
        t_step(k) = record.t_step(in_table);
        n_times(k) = size(sim_list{k,1},1);
    end

    v_scale = sqrt(g.*h);
    t_scale = sqrt(h./g);
    z_scale = h;
    rho = density_ratio.*phi_c+1-phi_c;

    d_dl = d./z_scale;
    dz_dl = dz./z_scale;

    % Creates cells to store the matrices of values

    buoyancy = zeros(sim_num);
    t_vals = cell(sim_num,1);
    z_pe = cell(sim_num,1);
    z_u = cell(sim_num,1);
    p_b = cell(sim_num,1);
    p_e = cell(sim_num,1);
    p_p = cell(sim_num,1);
    phi = cell(sim_num,1);
    u_f = cell(sim_num,1);
    u_p = cell(sim_num,1);
    Iv = cell(sim_num,1);
    Iv_phi = cell(sim_num,1);
    phi_Iv = cell(sim_num,1);
    mu_I = cell(sim_num,1);
    mu_Iv = cell(sim_num,1);
    beta_pe = cell(sim_num,1);
    beta_u = cell(sim_num,1);
    tau_p = cell(sim_num,1);
    tau_f = cell(sim_num,1);
    drag_term_p = cell(sim_num,1);
    drag_term_f = cell(sim_num,1);
    dpdz = cell(sim_num,1);
    dufdz = cell(sim_num,1);
    d2ufdz2 = cell(sim_num,1);
    dupdz = cell(sim_num,1);
    dilatancy = cell(sim_num,1);
    diffusion_term = cell(sim_num,1);
    tan_psi = cell(sim_num,1);
    dpdt = cell(sim_num,1);
    dphidt = cell(sim_num,1);
    dupdt = cell(sim_num,1);
    dufdt = cell(sim_num,1);
    flow_depth = cell(sim_num,1);

    for i = 1:sim_num
        buoyancy(i) = -(density_ratio(i)-1)*phi_c(i)*cosd(theta(i));
        z_pe{i,1} = linspace(1/(2*N(i)),1,N(i))';
        z_u{i,1} = linspace(0,1-1/(2*N(i)),N(i))';
        p_b{i,1} = (density_ratio(i)-1)*phi_c(i)*cosd(theta(i))*(1-z_pe{i,1});

        vec = sim_list{i,1};
        if (mod(size(vec,2),N(i)))
            t_vals{i,1} = vec(:,1);
            vec = vec(:,2:end);
        else
            t_vals{i,1} = linspace(0,(n_times(i)-1)*t_step(i),n_times(i));
        end
        %     tlength=[sum(t_vals{i,1}<1000)];
        %     t_vals{i,1} = t_vals{i,1}(1:tlength+1);
        %     vec = vec(1:tlength+1,:);
        %     n_times(i) = tlength+1;

        p_e{i,1} = vec(:,1:N(i))';
        p_p{i,1} = p_b{i,1}-p_e{i,1};
        phi{i,1} = vec(:,N(i)+1:2*N(i))';
        u_f{i,1} = vec(:,2*N(i)+1:3*N(i))';
        u_p{i,1} = vec(:,3*N(i)+1:end)';

%         flow_depth{i,1} = sum(u_p{i,1} > 1e-3,1)/N(i);
%         beta_pe{i,1} = 150*(phi{i,1}).^2.*eta_f_dl(i)./((1-phi{i,1}).^3.*d_dl(i)^2);
%         beta_u{i,1} = interp1(z_pe{i,1},beta_pe{i,1},z_u{i,1},'linear','extrap');
%         % If the velocity evolves we need to define these quantities
%         dpdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_e{i,1}.*ones(N(i),n_times(i)),1,1))./dz_dl(i);
%         dufdz{i,1} = vertcat(diff(u_f{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
%         d2ufdz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dufdz{i,1},1,1))./dz_dl(i);
% 
%         dupdz{i,1} = vertcat(diff(u_p{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
%         Iv{i,1} = eta_f_dl(i).*abs(dupdz{i,1})./(p_p{i,1}+1e-8);
%         tau_f{i,1} = eta_f_dl(i).*d2ufdz2{i,1}./(1-phi_c(i));
% 
%         mu_Iv{i,1} = mu_Iv_fn(Iv{i,1});
%         tau_p{i,1} = 1/(phi_c(i).*density_ratio(i)).*sign(dupdz{i,1}).*vertcat(zeros(1,n_times(i)),diff(mu_Iv{i,1}.*(p_p{i,1}.*ones(1,n_times(i))),1,1)./dz_dl(i));
%         drag_term_p{i,1} = (1-phi_c(i))^2.*beta_u{i,1}./(density_ratio(i)*phi_c(i)).*(u_f{i,1}-u_p{i,1});
%         drag_term_f{i,1} = (1-phi_c(i)).*beta_u{i,1}.*(u_f{i,1}-u_p{i,1});
%         dupdt{i,1} = vertcat(zeros(1,n_times(i)),tau_p{i,1}(2:end,:)+sind(theta(i))+drag_term_p{i,1}(2:end,:));
%         dufdt{i,1} = vertcat(zeros(1,n_times(i)),tau_f{i,1}(2:end,:)+sind(theta(i))-drag_term_f{i,1}(2:end,:));
% 
%         dpdt{i,1} = 1./(alpha_dl(i)).*vertcat(diff(1./beta_pe{i,1}.*dpdz{i,1}),zeros(1,n_times(i)))./dz_dl(i);
%         phi_Iv{i,1} = phi_c(i)./(1+sqrt(abs(Iv{i,1})));
%         Iv_phi{i,1} = (phi{i,1}-phi_c(i)).^2./phi{i,1}.^2;
%         tan_psi{i,1} = phi{i,1}-phi_Iv{i,1};
%         dilatancy{i,1} = -1./(alpha_dl(i)).*dupdz{i,1}.*(phi{i,1}-phi_Iv{i,1});
%         diffusion_term{i,1} = dpdt{i,1};
%         dphidt{i,1} = -dpdt{i,1}.*alpha_dl(i).*phi_c(i);
%         dpdt{i,1} = dpdt{i,1} + dilatancy{i,1};
    end
%% 
%     theta0 = 13;
%     % Fast Timescale
%     Iv_theta0 = newt_solve_crit_Iv(theta0,density_ratio,1);
%     pp_theta0 = (rho-1)*cosd(theta0)*(1-z_pe{i,1});
%     pp0 = (rho*cosd(theta)-cosd(theta0))*(1-z_pe{i,1});
%     Iv0 = Iv_theta0./pp0.*pp_theta0;
%     tanpsi0 = phi{i,1}(:,1)-phi_c./(1+sqrt(Iv0));
%     tau_s = phi{i,1}(:,1).*Iv_theta0^(3/2)./(2*eta_f_dl*alpha_dl);
%     tan_psi_fast = tanpsi0.*exp(-tau_s.*t_vals{i,1}');
%     pp_fast = (rho-1)*cosd(theta0).*(phi{i,1}(:,1)-phi_c(1)).^2./(tan_psi_fast-(phi{i,1}(:,1)-phi_c(1))).^2.*(1-z_pe{i,1});
%     pe_fast = p_b{i,1}-pp_fast;
%         
%     % Medium Timescale
%     dudz_ss_grad = rho(i).*sind(theta(i))./mu_Iv{i,1}.*(1-z_pe{i,1});
%     dudz_ss = rho(i).*sind(theta(i).*(phi{i,1}(100,1)-phi_c(i)).^2)./(phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1)).*(1-z_pe{i,1});
%     fric_mult = (phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1))./(phi{i,1}(100,1)-phi_c(i)).^2;
%     dudz_init = (rho(i)-1).*cosd(theta0).*(phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i));
%     cos_const_term = (phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i)).*((rho(i)-1).*cosd(theta0)-rho(i).*sind(theta(i))./mu_Iv{i,1}(100,1));
%     dudz_medium = dudz_ss;
%     for k =1:10
%         dudz_medium = dudz_medium + 8./(2*k-1).^2./pi.^2.*cos_const_term.*exp(-fric_mult(:,1)./rho(i).*pi^2/4.*(2*k-1).^2.*t_vals{i,1}').*cos(pi/2.*(2*k-1).*z_pe{i,1});
%     end
%     pp_medium = dudz_medium.*eta_f_dl.*phi{i,1}(1,100)^2./(phi{i,1}(1,100)-phi_c(i))^2;
%     pe_medium = p_b{i,1} - pp_medium;
%     
%     % Critical Pressure
%     crit_grad = -(density_ratio(i)*phi_c(i)+(1-phi_c(i)))*sind(theta(i))/0.32;
%     crit_pe = p_b{i,1}+crit_grad*(1-z_pe{i,1});
    
    nline=1;
    C = viridis(4);
    plot_vec = p_e{i,1};
    t_init = 11000;
    t_step = 0.1;
    t_maxes = linspace(t_init,t_init+t_step*(nline-1),nline);
%     SetPaperSize(8,8)
    hold on
    for j=linspace(1,nline,nline)
%         t_vals{j,1}(t1(j))
        t_max = t_maxes(j);
%         t_max = [50000,100000];
%         t_val = (k-1)*t_step(i);
        t1=[sum(t_vals{i,1}<t_max)+1];
%         dil_smooth = smooth(dilatancy{i,1}(1:end,t1));
%         plot(dil_smooth,z_pe{i,1}(1:end),'color',C(1,:),'DisplayName',"Dilatant Term")
%         plot(diffusion_term{i,1}(1:end,t1),z_pe{i,1}(1:end),'color',C(2,:),'DisplayName',"Diffusive Term")
        plot(plot_vec(1:end,t1),z_pe{i,1}(1:end),'color',C(j,:),'DisplayName',"t="+num2str(t_max))  
%         plot(smooth(plot_vec(1:end,t1)),z_pe{i,1}(1:end),'color',C(3,:),'DisplayName',"Full Profile")
%         plot(pe_fast(1:end,t1),z_pe{i,1}(1:end),'LineStyle','--','color',C(j,:),'DisplayName',"Approximation")
%         plot(dudz_medium(1:end,t1),z_pe{i,1}(1:end),'LineStyle','--','color',C(j,:),'DisplayName',"Approximation")
    end
%     plot(crit_pe,z_pe{i,1}(1:end),'LineStyle','--','color','r','DisplayName',"$p_{crit}$")
    legend('Location','northeast');
%     xlabel('$\frac{\partial p_e}{\partial t}$')
    xlabel('$\phi$')
    ylabel('$z$')
    title("t="+num2str(t_max))
%     xlim([0,0.045]);
    figname = "Ive_5deg_slow_cancelation.pdf";
%     exp_graph(gcf,figname)
%     movefile(figname,"Figures/SecondYearReport")
    function beta_val=beta_fn(phi)
        beta_val = 150*phi.^2.*eta_f_dl./((1-phi).^3.*d_dl^2);
    end

    function mu_val = mu_I_fn(I)
        mu1_I=0.342; % \frac{u+u_d}{\rho_f \phi_c}
        mu2_I=0.557; % 
        I_0 = 0.069;
        mu_val = tanh(reg_param*I).*(mu1_I+(mu2_I-mu1_I)./(1+I_0./abs(I)));
    end

    function mu_val = mu_Iv_fn(Iv)
        mu1_Iv = 0.32;
        mu2_Iv = 0.7;
        Iv_0 = 0.005;

        reg_param = 1*10^7;
        mu_val = tanh(reg_param*Iv).*(mu1_Iv+(mu2_Iv-mu1_Iv)./(1+Iv_0./abs(Iv))+Iv+5/2*phi_c*sqrt(Iv));
    end
end