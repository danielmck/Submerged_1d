function [xi_final,y_final] = vis_to_sinful(custom_init,params,master_xi,master_y)
% Does not work, an attempt to convert the no pe wave with no viscosity
% into the full wave
% Likely doesn't work as there is no way to find the flow height gradient
% at the critical point.
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    if ~exist("custom_init","var")
        custom_init = false;
    end
    % If there are custom initial parameters, the function may have to call
    % viscous_Iv_bvp_from_master to make the input
    if custom_init
        param_cell = num2cell(params);
        [Fr,theta,h_alt,alpha,d,tau0] = param_cell{:};
    else
        master_name = "non_vis_convert.txt";
        master_file = load("../Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('../Results/wave_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        theta = record.theta(in_table);
        lambda = record.lambda(in_table);
        Fr = record.Fr(in_table);
        tau0 = record.tau0(in_table);
        stat_len = 0;
        if reverse
            alpha = record.alpha(in_table);
            d = record.d(in_table);
        else
            alpha = 1e-5;
            d = 1e-4;
        end
    end
    
    
    rho = rho_p*phi_c+rho_f*(1-phi_c);
    P = (rho-rho_f)/rho;
    chi = (rho_f+3*rho)/(4*rho);
    
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    
    crit_pb = rho_f*g*cosd(theta)*h0;
    
    p_tot = rho*g*cosd(theta);
    
    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    tau0_dl = tau0/p_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 
    
    u_eq_dl = u_eq/v_scale;
    
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale;
    d_dl = d/z_scale;
    
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    
    
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    
%     stable = viscous_stability(theta,Fr,nu,lambda);

    h_pert = 5e-4;

    master_h = master_y(3,:);
    [~,max_ind] = max(master_h);
    
    master_y = master_y(:,1:max_ind);
    master_y = vertcat(master_y(1:3,:),master_y(5:end,:));
    master_xi = master_xi(1:max_ind);
    lambda = master_xi(end);
    master_xi = master_xi/master_xi(end);
    
    master_h = master_h(1:max_ind);
    master_Q1 = master_y(2,:);
    h_alt = master_h(1)-master_Q1(1)/master_y(1,1);
    crit_vec = master_h.^3-master_Q1.^2*Fr^2;
    root_ind = sum(crit_vec<0);
    root_xi = master_xi(root_ind) + (1-crit_vec(root_ind+1)/(crit_vec(root_ind+1)-crit_vec(root_ind)))*(master_xi(root_ind+1)-master_xi(root_ind));
    h_crit = master_h(root_ind) + (1-crit_vec(root_ind+1)/(crit_vec(root_ind+1)-crit_vec(root_ind)))*(master_h(root_ind+1)-master_h(root_ind));
    
    master_y_l = interp1(master_xi,master_y',linspace(0,root_xi,100))';
    master_y_r = interp1(master_xi,master_y',linspace(root_xi,1,100))';
    
    xi_in = linspace(0,1,100);
    y_in = vertcat(master_y_l(1,:),lambda*ones(1,100),root_xi*lambda*ones(1,100),master_y_l(2:end,:),master_y_r(2:end,:));
    
    tol = 1e-3;
    h_p_low = 1e-3;
    h_p_step = 5e-4;
%     while h_p_step > 5e-6
        h_p_new = max(h_p_low-h_p_step,h_p_low/2);
        y_in(5,end) = h_crit - h_p_new;
        y_in(10,1) = h_crit + h_p_new;
        
        solInit1=bvpinit(xi_in,@bvp_guess);
        opts = bvpset('RelTol',tol,'NMax',1000);

        solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1,opts);
        resid = solN1.stats.maxres;
        if resid<tol
            y_in = solN1.y;
            xi_in = solN1.x;
            h_p_low = h_p_new;
        else
            h_p_step = min(h_p_low/3, h_p_step/2);
        end
%     end
    h_p_low
    
    out_vec = vertcat(xi_in,y_in);
    filename = "no_vis_full_master.txt";
    save("Results/"+filename,"out_vec","-ascii")
    write_record("Results/wave_record.csv",filename,{"full","water",Fr,theta,alpha,d,tau0})
        
    function guess = bvp_guess(xi)
        % Initial guess function from the ode soln
        guess_ind = sum(xi_in<xi)+1;
        if guess_ind == max(size(xi_in))
            guess = y_in(:,end);
        else
            gap = xi_in(guess_ind+1)-xi_in(guess_ind);
            guess = y_in(:,guess_ind)*(xi-xi_in(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
        end  
    end

    function dydxi = viscous_syst(xi,y)
        % The system from the viroulet paper with a different friction
        % law. Need to solve for u_w and Q1 too.
        u_w = y(1);
        lambda = y(2);
        h_crit_pos = y(3);

        Q1_l = y(4);
        h_l = y(5);
        u_l = (-Q1_l + h_l.*u_w)./h_l;
        m_l = y(6);
        phi_l = y(7)./Q1_l;
        pb_l = y(8) + g_dl*cosd(theta)*rho_dl*chi.*h_l;

        zeta_l = 3/(2*alpha_dl*h_l) + P/4;
        p_p_l = p_tot_grad_dl*h_l-pb_l;
        Iv_l = 3*eta_f_dl*abs(u_l)/h_l/p_p_l;
        D_l = -2/beta_dl/h_l*(pb_l-h_l);

        denom_l = (h_l.^3/Fr^2-Q1_l.^2);
        fric_coeff_l = p_p_l/(p_tot_grad_dl*h_l)*mu_Iv_fn(Iv_l);
        fb_val_l = tand(theta)-fric_coeff_l-tau0_dl*rho_f/rho/h_l;
        dhdxi_l = 1/Fr^2.*h_l.^3.*(fb_val_l+P*D_l*h_l)./denom_l;

        R_w3_l = -phi_l*rho_f_dl/rho_dl*D_l;
        R_w4_l = (-P*chi+zeta_l)*D_l - 2*3/alpha_dl/h_l*u_l*(phi_l - phi_c./(1+sqrt(Iv_l)));

        dQdxi_l = -P*D_l;
        dmdxi_l = h_l/(lambda+stat_len)*u_l;
        dy6dxi_l = -R_w3_l;
        dy7dxi_l = R_w4_l/(u_l-u_w);
        dydxi_l = [dQdxi_l,dhdxi_l,dmdxi_l,dy6dxi_l,dy7dxi_l]*h_crit_pos;

        Q1_r = y(9);
        h_r = y(10);
        u_r = (-Q1_r + h_r.*u_w)./h_r;
        m_r = y(11);
        phi_r = y(12)./Q1_r;
        pb_r = y(13) + g_dl*cosd(theta)*rho_dl*chi.*h_r;

        zeta_r = 3/(2*alpha_dl*h_r) + P/4;
        p_p_r = p_tot_grad_dl*h_r-pb_r;
        Iv_r = 3*eta_f_dl*abs(u_r)/h_r/p_p_r;
        D_r = -2/beta_dl/h_r*(pb_r-h_r);

        denom_r = (h_r.^3/Fr^2-Q1_r.^2);
        fric_coeff_r = p_p_r/(p_tot_grad_dl*h_r)*mu_Iv_fn(Iv_r);
        fb_val_r = tand(theta)-fric_coeff_r-tau0_dl*rho_f/rho/h_r;
        dhdxi_r = 1/Fr^2.*h_r.^3.*(fb_val_r+P*D_r*h_r)./denom_r;

        R_w3_r = -phi_r*rho_f_dl/rho_dl*D_r;
        R_w4_r = (-P*chi+zeta_r)*D_r - 2*3/alpha_dl/h_r*u_r*(phi_r - phi_c./(1+sqrt(Iv_r)));

        dQdxi_r = -P*D_r;
        dmdxi_r = h_r/(lambda+stat_len)*u_r;
        dy6dxi_r = -R_w3_r;
        dy7dxi_r = R_w4_r/(u_r-u_w);
        dydxi_r = [dQdxi_r,dhdxi_r,dmdxi_r,dy6dxi_r,dy7dxi_r]*(lambda-h_crit_pos);

        dydxi = [0,0,0,dydxi_l,dydxi_r]';
    end
        
    function resid = bc_vals(ya,yb)
        h_crit = (yb(4)*Fr)^(2/3);
        h_start = ya(4)/ya(1)+h_alt;
        h_stop = (-h_start + sqrt(h_start.^2+8*ya(4)^2./(h_start/Fr^2)))/2;
        resid = [ya(4)-yb(9), ya(5)-h_start, yb(10)-h_stop, ya(6), yb(11)-1,...
            (ya(7)-yb(12)), (ya(8)-yb(13)), yb(4)-ya(9), yb(6)-ya(11), yb(7)-ya(12), yb(8)-ya(13), ...
            yb(5)-(h_crit-h_pert),ya(10)-(h_crit+h_pert)]; 
    end
end