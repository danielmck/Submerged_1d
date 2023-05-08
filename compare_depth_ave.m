function compare_depth_ave
% Compares the depth averaged models with the full models and different
% closure relations
    theta = 5.0;
    n_times = 5001;
    N=200;
    z_pe = linspace(1/(2*N),1,N)';
    z_u = linspace(0,1-1/(2*N),N)';
    q=floor(theta);
    r=round(mod(10*theta,10));
    
    if (r==0)
        deg_string = num2str(q)+"_deep";
    else
        deg_string = num2str(q)+"_"+num2str(r)+"_deep";
    end
    Rau_name = "Rauter_Closure/EqnOfState_Results/Rauter_"+deg_string+".txt";
    Ive_name = "Iverson_Closure/Results/Ive_comp_"+deg_string+".txt";
    da_name = "Iverson_DA/DA_Results/Ive_da_"+deg_string+"_v3.txt";
    da_name3 = "Iverson_DA/DA_Results/Ive_da_"+deg_string+"_v4.txt";
    
    h = 4e-2;
    eta_f = 0.0010016; % Pa s
    g=9.81; % m/s^2
    rho_f = 1000;
    shear_lim = 5;
    phi_rlp = 0.53;
    phi_rcp = 0.63;
    phi_c=0.585;
    density_ratio=2.5;
    a = 20;
    p_scale = rho_f.*g.*h;
    t_scale = sqrt(h./g);
    eta_f_dl = eta_f/(p_scale*t_scale);
    a_dl = a/p_scale;
    shear_lim_dl = shear_lim*t_scale;
    
    z_pe = linspace(1/(2*N),1,N)';
    z_u = linspace(0,1-1/(2*N),N)';
    p_b = (density_ratio-1)*phi_c*cosd(theta)*(1-z_pe);
    
    Rau_file = load(Rau_name);
    Rau_up = Rau_file(:,2*N+1:end);
    Rau_phi = Rau_file(:,1);
    
    dz = 1/(N-0.5);
    Rau_dupdz = diff(Rau_up,1,2)/dz;
    Rau_dupdz = horzcat(Rau_dupdz,zeros(n_times,1));
    Rau_dupdz_base = Rau_dupdz(:,1);
    Rau_is_flow = ((Rau_up > 1e-4) | (Rau_dupdz > 2e-4));
    Rau_h_vals = sum(Rau_is_flow,2)/N;
    Rau_erosion_start = sum(Rau_h_vals>0.999);
    Rau_u_vals = depth_average(Rau_up',N,n_times);
    
    Rau_pc = a_dl*(phi_c+Rau_phi-phi_rlp)./(phi_rcp-phi_c-Rau_phi);
    phi_m = phi_c+(phi_rcp-phi_c).*(abs(Rau_dupdz_base)<shear_lim_dl).*(shear_lim_dl-abs(Rau_dupdz_base)).^2/shear_lim_dl.^2;
    Rau_pi = eta_f_dl.*Rau_dupdz_base./((phi_m./(phi_c+Rau_phi)-1).^2);
    Rau_pe = (density_ratio-1)*phi_c*cosd(theta) - Rau_pc - Rau_pi;
    
    Ive_file = load(Ive_name);
    Ive_up = Ive_file(:,3*N+1:end);
    Ive_pe_all = Ive_file(:,1:200);
    Ive_pp_all = p_b' - Ive_pe_all;
    Ive_pe_base = Ive_file(:,1);
    Ive_pp_base = (density_ratio-1)*phi_c*cosd(theta)-Ive_pe_base;
    Ive_dupdz = diff(Ive_up,1,2)/dz;
    Ive_dupdz = horzcat(Ive_dupdz,zeros(n_times,1));
    Ive_dupdz_base = Ive_dupdz(:,1);
    Ive_Iv = Ive_dupdz.*eta_f_dl./(Ive_pp_all);
    Ive_Iv_base = Ive_dupdz_base.*eta_f_dl./((density_ratio-1)*phi_c*cosd(theta)-Ive_pe_base);
    Ive_is_flow = ((Ive_up > 1e-4) | (Ive_dupdz > 2e-4));
    Ive_h_vals = sum(Ive_is_flow,2)/N;
    Ive_erosion_start = sum(Ive_h_vals>0.999);
    Ive_u_vals = depth_average(Ive_up',N,n_times);
    
    DA_file = load(da_name);
    DA_h_vals = DA_file(:,1);
    DA_phi_vals = DA_file(:,2);
    DA_u_vals = DA_file(:,3);
    DA_pe_vals = DA_file(:,4)-cosd(theta)*DA_h_vals;
    DA_dupdz_base = 3*DA_u_vals;
    DA_Iv_base = 3*DA_u_vals*eta_f_dl./(2*h*((density_ratio-1)*h*cosd(theta)-DA_pe_vals));
    
    
    DA_file3 = load(da_name3);
    DA_h_vals3 = DA_file3(:,1);
    DA_phi_vals3 = DA_file3(:,2);
    DA_u_vals3 = DA_file3(:,3);
    DA_pe_vals3 = DA_file3(:,4)-cosd(theta)*DA_h_vals3;
    DA_pp_vals3 = ((density_ratio-1)*DA_phi_vals3.*DA_h_vals3*cosd(theta)-DA_pe_vals3);
    DA_dupdz_base3 = 5/2*DA_u_vals3;
    DA_Iv_base3 = 5*DA_u_vals3*eta_f_dl./(2*DA_h_vals3.*DA_pp_vals3);

    %% 
    t_max = 100;
    hold on
%     SetPaperSize(10,10)
%     plot(linspace(1,t_max,t_max+1),Rau_u_vals(1:t_max+1),'DisplayName',"Rauter Model")
    plot(linspace(0,10*t_max,t_max+1),Ive_u_vals(1:t_max+1),'DisplayName',"Iverson Model")
    plot(linspace(0,10*t_max,t_max+1),DA_u_vals(1:t_max+1),'DisplayName',"DA Model")
    plot(linspace(0,10*t_max,t_max+1),DA_u_vals3(1:t_max+1),'DisplayName',"DA Model with Adapted ICs")
%     plot(Ive_pe_all(floor(t_max),:),z_pe);
%     plot(DA_pe_vals3(floor(t_max)),0,'x');
%     top_grad = 2*Ive_pe_all(floor(t_max),100);
%     plot(top_grad*(1-z_pe),z_pe,'--')

%     plot(Rau_up(850,:),z_u,'DisplayName',"Rauter Model");
%     plot(Ive_up(850,:),z_u,'DisplayName',"Iverson Model");
%     ylim([0,4])
    xlabel('$t$');
    ylabel('Basal $p_e$');
    title('Full Iverson Model vs Depth Averaged Models with two sets of ICs');
    legend('Location', "best",'UserData', 18);
%     PrintFig('Ive_da_adapt_IC_pre_depo_pe')
    
end