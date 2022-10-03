function rate_comp_plot
    z_model = "Ive_5deg_13init_deposit.txt";
    N = 200;
    vec = load("Results/"+z_model);
    theta = 5;
    t_z = vec(:,1);
    n_t_z = size(t_z,1);
    vec = vec(:,2:end);
    p_e_z = vec(:,1)';
    phi_z = depth_average(vec(:,N+1:2*N)',N,n_t_z);
    u_f_z = depth_average(vec(:,2*N+1:3*N)',N,n_t_z);
    u_p_z = depth_average(vec(:,3*N+1:end)',N,n_t_z);

    da_model = "Ive_da_5deg_13init_long.txt";
    data_file = load("../Iverson_DA/Results/"+da_model);
    t_da = data_file(:,1);
    n_t_da = size(t_da,1);
    h_da = data_file(:,2);
    phi_da = data_file(:,3);
    u_da = data_file(:,4);
    pb_da = data_file(:,5);
    pe_da = pb_da - cosd(theta)*h_da;
 %%   
    t_min = 500;
    t_max = 5000;
    
    start_ind_z = max(sum(t_z<t_min),1);
    end_ind_z = min(sum(t_z>t_min)+1,n_t_z);
    
    start_ind_da = max(sum(t_da<t_min),1);
    end_ind_da = min(sum(t_da>t_min)+1,n_t_da);
    
    SetPaperSize(8,8)
    C = viridis(4);
    
    hold on
    plot(t_z(start_ind_z:end_ind_z),u_p_z(start_ind_z:end_ind_z),'color',C(1,:),'DisplayName',"Full Model")
    plot(t_da(start_ind_da:end_ind_da),u_da(start_ind_da:end_ind_da),'color',C(3,:),'DisplayName',"Depth Averaged Model")
    xlim([500,5000])
    xlabel("$p_e$")
    legend("Location","best")
    ylabel("$\bar{u}$")
    figname = "Ive_5deg_u_comp_slow.pdf";
    exp_graph(gcf,figname)
    movefile(figname,"Figures/SecondYearReport")
end