    full_name = "Ive_5deg_Fr_10_mud.txt";
    da_name = "Ive_da_5deg_Fr_10_mud_early.txt";


    [phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
    
    
    cd ../Iverson_Closure/Results
    record = readtable('result_record.csv');
    in_table = strcmp(record.Name, full_name);
    full_vec = load(full_name);
    cd ../../Iverson_DA
    sim_type = record.sim_type(in_table);
    N = record.N(in_table);
    h = record.h(in_table);
    d = record.d(in_table);
    density_ratio = record.rho_r(in_table);
    eta_f_dl = record.eta_f(in_table);
    theta = record.theta(in_table);
    alpha_dl = record.alpha(in_table);
    t_step = record.t_step(in_table);
    full_n_times = size(full_vec,1);
    dz = 1/(N-0.5); % z spacing
    z_pe = linspace(dz/2,1,N);
    p_tot = (density_ratio-1)*phi_c*cosd(theta)*(1-z_pe);

    cd Results
    da_vec = load(da_name);
    cd ../
    
    if (mod(size(full_vec,2),N))
        t_vals = full_vec(:,1);
        full_vec = full_vec(:,2:end);
    else
        t_vals = linspace(0,(full_n_times-1)*t_step,full_n_times);
    end
    %     tlength=[sum(t_vals<1000)];
    %     t_vals = t_vals(1:tlength+1);
    %     vec = vec(1:tlength+1,:);
    %     n_times = tlength+1;
    n_times = length(t_vals);
    p_e = full_vec(:,1:N)';
    p_p = p_tot'.*ones(1,size(t_vals,1))-p_e;
    phi = full_vec(:,N+1:2*N)';
    u_f = full_vec(:,2*N+1:3*N)';
    u_p = full_vec(:,3*N+1:end)';
    dupdz = vertcat(diff(u_p,1),zeros(1,size(t_vals,1)))*(N-0.5);
    iv= dupdz*eta_f_dl./max(p_p,1e-8);
    tan_psi = phi-phi_c./(1+sqrt(iv));
    
    u_p_full_ave = depth_average(u_p,N,n_times);
    u_f_full_ave = depth_average(u_f,N,n_times);
    u_full_ave = phi_c*u_p_full_ave+(1-phi_c)*u_f_full_ave;
    phi_full_ave = depth_average(phi,N,n_times);
    tan_psi_full_ave = depth_average(tan_psi,N,n_times);
    pe_full_base = p_e(1,:);
    pp_full_base = p_p(1,:);
    
    pe_full_quart = p_e(N/4,:);
    pp_full_quart = p_p(N/4,:);
    
    da_file = load(da_name);
    t_vals_da = da_file(:,1);
    h_da = da_file(:,2);
    phi_da = da_file(:,3)./h_da;
    u_da = da_file(:,4);
    pb_da = da_file(:,5);
    pe_da = pb_da-h_da*cosd(theta);
    pp_da = (density_ratio-1)*phi_c*cosd(theta)*h_da-pe_da;
    pe_da_quart = pe_da*3/4;
    pp_da_quart = pp_da*3/4;
    
    iv_da = 3*u_da*eta_f_dl./pp_da./h_da;
    tan_psi_da = phi_da-phi_c./(1+sqrt(iv_da));
%%
    SetPaperSize(7.8,7.8)
%     SetPaperSize(10,10)
    hold on
    t_min = 0;
    t_max = 2;
    
    full_min = max(sum(t_vals<t_min),1);
    full_max = min(sum(t_vals<t_max)+1,size(t_vals,1));
    
    da_min = max(sum(t_vals_da<t_min),1);
    da_max = min(sum(t_vals_da<t_max)+1,size(t_vals_da,1));
    
    plot(t_vals(full_min:full_max), tan_psi_full_ave(full_min:full_max),'DisplayName',"Full model", 'color', '#fc8d62','LineWidth',1.3) %[0.127,0.201,0.127])
    plot(t_vals_da(da_min:da_max), tan_psi_da(da_min:da_max),'DisplayName',"DA model", 'color', '#8da0cb','LineWidth',1.3)
    legend('Location', "best",'UserData', 8);
    xlabel("$t$");
    ylabel('$\tan\bar{\psi}$');%,'Position',[-35 0.0045]);
    xlim([t_min t_max])
    exp_graph(gcf,"Comp_tanpsi_early_mud.pdf")
    