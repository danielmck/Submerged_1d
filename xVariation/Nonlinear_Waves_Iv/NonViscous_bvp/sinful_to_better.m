function [xi_final,y_final] = sinful_to_better(custom_init,reverse,params,provide_init,master_xi,master_y)
% Does not work, an attempt to convert the no pe wave with no viscosity
% into the full wave
% Likely doesn't work as there is no way to find the flow height gradient
% at the critical point.
    [phi_c,rho_f,rho_p,rho,eta_f,g] = get_params_water();
    P = (rho-rho_f)/rho;
    if ~exist("custom_init","var")
        custom_init = false;
    end
    if ~exist("reverse","var")
        reverse = false;
    end
    if ~exist("provide_init","var")
        provide_init = false;
    end
    % If there are custom initial parameters, the function may have to call
    % viscous_Iv_bvp_from_master to make the input
    if custom_init
        param_cell = num2cell(params);
        [Fr,theta,h_alt,alpha,d,tau0] = param_cell{:};
    else
        master_name = "no_vis_full_master.txt";
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/wave_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        theta = record.theta(in_table);
%         lambda = record.lambda(in_table);
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
    h_crit = (master_y(5,end)+master_y(10,1))/2;
    master_y(5,end) = h_crit;
    master_y(10,1) = h_crit;
    g_l = (master_y(5,end-1)-master_y(5,end-2))/(master_xi(end-1)-master_xi(end-2))*master_y(3,1);
    g_r = (master_y(10,3)-master_y(10,2))/(master_xi(3)-master_xi(2))*(master_y(2,1)-master_y(3,1));
    cp_grad =  (g_l+g_r)/2;
    h_alt = master_y(5,1)-master_y(4,1)/master_y(1,1);
    
    xi_in = horzcat(master_xi,master_xi+1);
    y_in = horzcat(master_y(1:8,:),vertcat(master_y(1:3,:),master_y(9:end,:)));
    y_in = vertcat(y_in,cp_grad*ones(size(xi_in)));
    
    tol=1e-4;
    xi_eps = 1e-3;
    solInit1=bvpinit(xi_in,@bvp_guess);
    opts = bvpset('RelTol',tol,'NMax',2000);

    solN1 = bvp4c(@(x,y,r) viscous_syst(x,y,r),@bc_vals,solInit1,opts);
    resid = solN1.stats.maxres;
    
    out_vec = vertcat(xi_in,y_in);
    filename = "no_vis_better_master.txt";
    save("Results/"+filename,"out_vec","-ascii")
    write_record("Results/wave_record.csv",filename,{"full","water",Fr,theta,alpha,d,tau0})
        
    function guess = bvp_guess(xi,region)
        % Initial guess function from the ode soln
        guess_ind = sum(xi_in<xi)+1;
        if guess_ind == max(size(xi_in))
            guess = y_in(:,end);
        else
            gap = xi_in(guess_ind+1)-xi_in(guess_ind);
            if abs(gap)<1e-8
                guess = y_in(:,guess_ind);
            else
                guess = y_in(:,guess_ind)*(xi-xi_in(guess_ind))/gap + y_in(:,guess_ind+1)*(xi_in(guess_ind+1)-xi)/gap;
            end
        end  
    end

    function dydxi = viscous_syst(xi,y, region)
        % The system from the viroulet paper with a different friction
        % law. Need to solve for u_w and Q1 too.
        u_w = y(1);
        lambda = y(2);
        h_crit_pos = y(3);

        Q1 = y(4);
        h = y(5);
        u = (-Q1 + h.*u_w)./h;
        m = y(6);
        phi = y(7)./Q1;
        pb = y(8) + g_dl*cosd(theta)*rho_dl*chi.*h;
        crit_h_grad = y(9);

        zeta = 3/(2*alpha_dl*h) + P/4;
        p_p = p_tot_grad_dl*h-pb;
        Iv = 3*eta_f_dl*abs(u)/h/p_p;
        D = -2/beta_dl/h*(pb-h);

        denom = (h.^3/Fr^2-Q1.^2);
        fric_coeff = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv);
        fb_val = tand(theta)-fric_coeff-tau0_dl*rho_f/rho/h;
        dhdxi = 1/Fr^2.*h.^3.*(fb_val+P*D*h)./denom;
        
        R_w3 = -phi*rho_f_dl/rho_dl*D;
        R_w4 = (-P*chi+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

        dQdxi = -P*D;
        dmdxi = h/(lambda+stat_len)*u;
        dy6dxi = -R_w3;
        dy7dxi = R_w4/(u-u_w);
        dydxi = [0,0,0,dQdxi,dhdxi,dmdxi,dy6dxi,dy7dxi,0];
        
        dpbdxi = dy7dxi + g_dl*cosd(theta)*rho_dl*chi.*dhdxi;
        dppdxi = p_tot_grad_dl*dhdxi-dpbdxi;
        dudxi = -dQdxi./h+Q1/h^2*dhdxi;
        dIvdxi = 3*eta_f_dl*dudxi/h/p_p-3*eta_f_dl*abs(u)/h^2/p_p*dhdxi-3*eta_f_dl*abs(u)/h/p_p^2*dppdxi;
        ddenomdxi = (3*h.^2/Fr^2.*dhdxi-2*Q1.*dQdxi);
        dnumerdxi = (-dppdxi/(p_tot_grad_dl*h)+p_p/(p_tot_grad_dl*h^2)*dhdxi)*mu_Iv_fn(Iv)-p_p/(p_tot_grad_dl*h)*dmudIv_fn(Iv)*dIvdxi...
            +tau0_dl*rho_f/rho/h^2*dhdxi+P*2/beta_dl*(dpbdxi-dhdxi);
        d2hdxi2 = 3/Fr^2.*h.^2.*(fb_val+P*D*h)./denom.*dhdxi+(dnumerdxi*denom-(fb_val+P*D*h)*ddenomdxi)/denom^2;
        switch region
            case 1    % x in [0 1]
                dydxi = dydxi*h_crit_pos;
            case 2    % x in [1 2]
                dydxi = dydxi*(lambda-h_crit_pos);
        end
        if abs(xi-1)<xi_eps
            dydxi(5) = crit_h_grad;
        elseif abs(xi-1)<5*xi_eps && abs(d2hdxi2)>1
            dydxi(5) = crit_h_grad;
        end
        dydxi = dydxi';
    end
        
    function resid = bc_vals(ya,yb)
        h_start = ya(4,1)/ya(1,1)+h_alt;
        h_stop = (-h_start + sqrt(h_start.^2+8*ya(4,1)^2./(h_start/Fr^2)))/2;
        
        h_mid = ya(5,2);
        u_mid = (-ya(4,2) + h_mid.*ya(1,2))./h_mid;
        denom_mid = (h_mid.^3/Fr^2-ya(4,2).^2);
        pb_mid = ya(8,2) + g_dl*cosd(theta)*rho_dl*chi.*h_mid;
        p_p_mid = p_tot_grad_dl*h_mid-pb_mid;
        Iv_mid = 3*eta_f_dl*abs(u_mid)/h_mid/p_p_mid;
        fric_coeff_mid = p_p_mid/(p_tot_grad_dl*h_mid)*mu_Iv_fn(Iv_mid);
        fb_val_mid = tand(theta)-fric_coeff_mid-tau0_dl*rho_f/rho/h_mid;
        
        cont_resid = (ya(:,2) - yb(:,1))';  %[ya(1,2)-yb(1,1),ya(2,2)-yb(2,1),ya(3,2)-yb(3,1),ya(4,2)-yb(4,1), ...
             %ya(8,2)-yb(8,1),ya(8,2)-yb(8,1),ya(8,2)-yb(8,1),ya(8,2)-yb(8,1)]
        resid = [ya(4,1)-yb(4,2), ya(5,1)-h_start, yb(5,2)-h_stop, ya(6,1), yb(6,2)-1,...
            (ya(7,1)-yb(7,2)), (ya(8,1)-yb(8,2)), denom_mid, fb_val_mid, cont_resid]; 
    end
end