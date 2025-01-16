function simple_plotting
    names = ["Crit_pe_test.txt"]; %Ive_5deg_13initIve_5deg_13init.txt"Ive_comp_4_deep_1_flux.txt"];%,"Ive_comp_4_deep_9_2_start.txt""Ive_comp_4_deep_9_1_start.txt","
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

        flow_depth{i,1} = sum(u_p{i,1} > 1e-3,1)/N(i);
        beta_pe{i,1} = 150*(phi{i,1}).^2.*eta_f_dl(i)./((1-phi{i,1}).^3.*d_dl(i)^2);
        beta_u{i,1} = interp1(z_pe{i,1},beta_pe{i,1},z_u{i,1},'linear','extrap');
        % If the velocity evolves we need to define these quantities
        dpdz{i,1} = vertcat(zeros(1,n_times(i)),diff(p_e{i,1}.*ones(N(i),n_times(i)),1,1))./dz_dl(i);
        dufdz{i,1} = vertcat(diff(u_f{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
        d2ufdz2{i,1} = vertcat(zeros(1,n_times(i)),diff(dufdz{i,1},1,1))./dz_dl(i);

        dupdz{i,1} = vertcat(diff(u_p{i,1},1,1),zeros(1,n_times(i)))./dz_dl(i);
        Iv{i,1} = eta_f_dl(i).*abs(dupdz{i,1})./(p_p{i,1}+1e-8);
        tau_f{i,1} = eta_f_dl(i).*d2ufdz2{i,1}./(1-phi_c(i));

        mu_Iv{i,1} = mu_Iv_fn(Iv{i,1});
        tau_p{i,1} = 1/(phi_c(i).*density_ratio(i)).*sign(dupdz{i,1}).*vertcat(zeros(1,n_times(i)),diff(mu_Iv{i,1}.*(p_p{i,1}.*ones(1,n_times(i))),1,1)./dz_dl(i));
        drag_term_p{i,1} = (1-phi_c(i))^2.*beta_u{i,1}./(density_ratio(i)*phi_c(i)).*(u_f{i,1}-u_p{i,1});
        drag_term_f{i,1} = (1-phi_c(i)).*beta_u{i,1}.*(u_f{i,1}-u_p{i,1});
        dupdt{i,1} = vertcat(zeros(1,n_times(i)),tau_p{i,1}(2:end,:)+sind(theta(i))+drag_term_p{i,1}(2:end,:));
        dufdt{i,1} = vertcat(zeros(1,n_times(i)),tau_f{i,1}(2:end,:)+sind(theta(i))-drag_term_f{i,1}(2:end,:));

        dpdt{i,1} = 1./(alpha_dl(i)).*vertcat(diff(1./beta_pe{i,1}.*dpdz{i,1}),zeros(1,n_times(i)))./dz_dl(i);
        phi_Iv{i,1} = phi_c(i)./(1+sqrt(abs(Iv{i,1})));
        Iv_phi{i,1} = (phi{i,1}-phi_c(i)).^2./phi{i,1}.^2;
        tan_psi{i,1} = phi{i,1}-phi_Iv{i,1};
        dilatancy{i,1} = -1./(alpha_dl(i)).*dupdz{i,1}.*(phi{i,1}-phi_Iv{i,1});
        diffusion_term{i,1} = dpdt{i,1};
        dphidt{i,1} = -dpdt{i,1}.*alpha_dl(i).*phi_c(i);
        dpdt{i,1} = dpdt{i,1} + dilatancy{i,1};
    end
%% 
    theta0 = 10;
    % Fast Timescale
    Iv_theta0 = newt_solve_crit_Iv(theta0,density_ratio,1);
    pp_theta0 = (rho-1)*cosd(theta0)*(1-z_pe{i,1});
    pp0 = (rho*cosd(theta)-cosd(theta0))*(1-z_pe{i,1});
    Iv0 = Iv_theta0./pp0.*pp_theta0;
    tanpsi0 = phi{i,1}(:,1)-phi_c./(1+sqrt(Iv0));
    tau_s = phi{i,1}(:,1).*Iv_theta0^(3/2)./(2*eta_f_dl*alpha_dl);
    tan_psi_fast = tanpsi0.*exp(-tau_s.*t_vals{i,1}');
    pp_fast = (rho-1)*cosd(theta0).*(phi{i,1}(:,1)-phi_c(1)).^2./(tan_psi_fast-(phi{i,1}(:,1)-phi_c(1))).^2.*(1-z_pe{i,1});
    pe_fast = p_b{i,1}-pp_fast;
        
    % Medium Timescale
    dudz_ss_grad = rho(i).*sind(theta(i))./mu_Iv{i,1}.*(1-z_pe{i,1});
    dudz_ss = rho(i).*sind(theta(i).*(phi{i,1}(100,1)-phi_c(i)).^2)./(phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1)).*(1-z_pe{i,1});
    fric_mult = (phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv{i,1}(100,1))./(phi{i,1}(100,1)-phi_c(i)).^2;
    dudz_init = (rho(i)-1).*cosd(theta0).*(phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i));
    cos_const_term = (phi{i,1}(100,1)-phi_c(i)).^2./(phi{i,1}(100,1).^2.*eta_f_dl(i)).*((rho(i)-1).*cosd(theta0)-rho(i).*sind(theta(i))./mu_Iv{i,1}(100,1));
    dudz_medium = dudz_ss;
    for k =1:10
        dudz_medium = dudz_medium + 8./(2*k-1).^2./pi.^2.*cos_const_term.*exp(-fric_mult(:,1)./rho(i).*pi^2/4.*(2*k-1).^2.*t_vals{i,1}').*cos(pi/2.*(2*k-1).*z_pe{i,1});
    end
    pp_medium = dudz_medium.*eta_f_dl.*phi{i,1}(1,100)^2./(phi{i,1}(1,100)-phi_c(i))^2;
    pe_medium = p_b{i,1} - pp_medium;
    
    
    Fr = 5;
    % Fast Timescale
    Iv_theta0_new = (phi_c-phi{i,1}(1,1))^2/phi{i,1}(1,1)^2;
%     pp_theta0_new = (rho-1)*cosd(theta)*(1-z_pe{i,1});
%     pp0 = rho*cosd(theta)*(1-z_pe{i,1});
    Iv0_new = dupdz{i,1}(1,1)*eta_f_dl/p_p{i,1}(1,1);
    tanpsi0_new = phi{i,1}(:,1)-phi_c./(1+sqrt(Iv0_new));
    tau_s_new = phi{i,1}(:,1).*Iv_theta0_new^(3/2)./(2*eta_f_dl*alpha_dl);
    tan_psi_fast_new = tanpsi0_new.*exp(-tau_s_new.*t_vals{i,1}');
    pp_fast_new = eta_f_dl*dupdz{i,1}(:,1).*phi{i,1}(:,1).^2./(tan_psi_fast_new-(phi{i,1}(:,1)-phi_c(1))).^2;
    pe_fast_new = p_b{i,1}-pp_fast_new;
        
    % Medium Timescale
    dudz_ss_grad_new = rho(i).*sind(theta(i))./mu_Iv{i,1}.*(1-z_pe{i,1});
    dudz_ss_new = rho(i).*sind(theta(i).*(phi{i,1}(100,1)-phi_c(i)).^2)./(phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv_fn(Iv_theta0_new)).*(1-z_pe{i,1});
    fric_mult_new = (phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv_fn(Iv_theta0_new))./(phi{i,1}(100,1)-phi_c(i)).^2;
    cos_const_term_new = dupdz{i,1}(1,1)-rho(i).*sind(theta(i).*(phi{i,1}(100,1)-phi_c(i)).^2)./(phi{i,1}(100,1).^2.*eta_f_dl(i).*mu_Iv_fn(Iv_theta0_new));
    dudz_medium_new = dudz_ss_new;
    for k =1:10
        dudz_medium_new = dudz_medium_new + 8./(2*k-1).^2./pi.^2.*cos_const_term_new.*exp(-fric_mult_new(:,1)./rho(i).*pi^2/4.*(2*k-1).^2.*t_vals{i,1}').*cos(pi/2.*(2*k-1).*z_pe{i,1});
    end
    pp_medium_new = dudz_medium_new.*eta_f_dl./Iv_theta0_new;
    pe_medium_new = p_b{i,1} - pp_medium_new;
    % Critical Pressure
    crit_grad_new = -(density_ratio(i)*phi_c(i)+(1-phi_c(i)))*sind(theta(i))/0.32;
    crit_pe_new = (rho-1)*cosd(theta)*(1-z_pe{i,1})+crit_grad_new*(1-z_pe{i,1});
    
    nline=5;
    n_subs = 1;
    % C = viridis(nline+1);
    plot_vec = u_p{i,1};
%     plot_cell = {p_e{i,1},phi{i,1}};
    t_init = 0;
    t_step = 0.5;
    t_maxes = linspace(t_init,t_init+t_step*(nline-1),nline);
%     t_maxes(1)=2;
    SetPaperSize(15,10);
    plot(plot_vec(:,end),z_u{i,1})

    f = zeros(1,nline);
    a = zeros(1,2);
    hold on
%     a(1) = plot([NaN],[NaN],'color','k','DisplayName',"Full");
%     a(2) = plot([NaN],[NaN],'color','k','LineStyle','--','DisplayName',"Approx.");
%     for k = 1:n_subs
%         gca=subplot(1,3,3);
%         hold on
% %         plot_vec = plot_cell{1,k};
%         for j=linspace(1,nline,nline)
%     %         t_vals{j,1}(t1(j))
%             t_max = t_maxes(j);
%     %         t_max = [50000,100000];
%     %         t_val = (k-1)*t_step(i);
%             t1=[sum(t_vals{i,1}<t_max)+1];
%             f(j) = plot(plot_vec(1:end-1,t1),z_pe{i,1}(1:end-1),'color',C(j,:),'DisplayName',"t="+num2str(t_max));
% %             plot(pe_medium_new(:,t1),z_pe{i,1},'LineStyle','--','color',C(j,:),'HandleVisibility','off');
% %             plot(tan_psi_fast_new(1:end,t1),z_pe{i,1}(1:end),'LineStyle','--','color',C(j,:),'DisplayName',"Approximation")
%         end
% %         if k == 1
%         ylabel('$z$','Interpreter','latex')
%         xlabel('$\phi$','Interpreter','latex')
%         gca.Position(2)=0.22;
%         gca.Position(4)=0.7;
% %         else\frac{\partial u}{\partial z}
% %             xlabel('$\phi$','Interpreter','latex')
% %         end
%     end
    % plot(crit_pe_new(1)*ones(N,1),z_pe{i,1}(1:end),'LineStyle','--','color','r','DisplayName',"$p_{crit}$")
%     legend([a,f],'Location','southoutside'); %["Full Profile","Approximation","t=1","t=2"]
%     legend([f(1),f(2)],{"Full Profile","Approximation"},'Location','east');
%     leg=legend('Position',[0.65 0.65 0.1 0.2],'Interpreter','latex');
%     xlabel('$\frac{\partial u}{\partial z}$')
%     xlabel('$p_e$')
%     ylabel('$z$')
%     annotation('arrow',[0.25 0.45],...
%     [0.3 0.3]);
%     annotation('textbox',...
%     [0.477793248945147 0.245454545454545 0.387185654008439 0.168181818181818],...
%     'String',{'Increasing $t$'},...
%     'FitBoxToText','off');
%     title("t="+num2str(t_max))
    ylim([0.4,0.5]);
    figname = "Quad_pe_u_.pdf";
%     exportgraphics(gcf,figname,"Resolution",300)
    exp_graph(gcf,figname)
    clf()
%     movefile(figname,"Figures/BAMC_Talk")


end