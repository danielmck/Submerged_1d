function chop_end
    filename = 'master_pres_h.txt';
    master_file = load("Results/"+filename);
    xi_in = master_file(1,:);
    y_in = master_file(2:end,:);
    record = readtable('Results/full_record.csv');

    in_table = strcmp(record.Name, filename);
    wave_type = record.wave_type(in_table);
    theta = record.theta(in_table); 
    Fr = record.Fr(in_table);
    tau0 = record.tau0(in_table);
    full_model=true;
    d = record.d(in_table);
    alpha = record.alpha(in_table);
    u_w = record.u_w(in_table);
    lambda = record.lambda(in_table);
    xi_crit = record.crit_xi(in_table);
    
    [phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0,false,true);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
%     p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;

%     h_stop = tau0/(rho*g*cosd(theta))/(tand(theta)-(rho-rho_f)/rho*mu1_Iv);


    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;
    
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    g_dl = g*t_scale/v_scale; 
    
    Q1 = y_in(1,:);
    h = y_in(2,:);
    u = u_w - Q1./h;
    m = y_in(3,:);
    phi = y_in(4,:)./Q1;
    rho = rho_p_dl*phi+rho_f_dl*(1-phi);
    chi = (rho_f_dl+3*rho)/(4*rho);
    pb = y_in(5,:) + rho.*g_dl*cosd(theta).*chi.*h;

    [~, max_ind] = max(pb);
    xi_max = xi_in(max_ind);
    xi_out = xi_in(1:max_ind);
    xi_out(xi_out>1) = 1+(xi_out(xi_out>1)-1)/(xi_max-1);
    lambda_out = xi_crit+(lambda-xi_crit)*(xi_max-1);
    y_out = y_in(:,1:max_ind);
    
    out_vec = vertcat(xi_out,y_out);
    out_name = strcat(filename(1:end-4),'_chop.txt');
    save("Results/"+out_name,"out_vec","-ascii")
    write_record("Results/full_record.csv",out_name,{"chop_pres_h","water",Fr,theta,alpha,d,tau0,u_w,lambda,xi_crit})
end