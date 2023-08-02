function [xi_final,y_final] = convert_to_rauter_nonvis(custom_init,params,provide_init,master_xi,master_y)
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
        master_name = "nonvis_convert.txt";
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/full_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        pres_h = 0; %(wave_type=="full_pres_h");
        theta = record.theta(in_table);
        lambda_vis = record.lambda(in_table);
        Fr = record.Fr(in_table);
        tau0 = record.tau0(in_table);
        master_u_w = record.u_w(in_table);
        stat_len = 0;
        a = record.a(in_table);
        d = record.d(in_table);
        phi_rlp = record.phi_rlp(in_table);
        phi_rcp = record.phi_rcp(in_table);
    end
    
    [h0, crit_Iv] = crit_Iv_rauter(theta, rho_p, rho_f, eta_f, Fr, a, phi_rlp, phi_rcp, tau0, false);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    pp0 = 3*eta_f*u_eq/h0/crit_Iv; 
    rho_eq = pp0/h0/g/cosd(theta)+rho_f;
    phi_eq = (rho_eq-rho_f)/(rho_p-rho_f);
    % rho_eq = phi0*rho_p+(1-phi0)*rho_f;

    p_tot = rho_eq*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;

    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    tau0_dl = tau0/p_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    a_dl = a/p_scale;
    g_dl = g*t_scale/v_scale; 

    u_eq_dl = u_eq/v_scale;
    
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_eq_dl = rho_eq*v_scale^2/p_scale;
    d_dl = d/z_scale;

    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    master_h = master_y(3,:);
    [~,max_ind] = max(master_h);
    
    master_y = master_y(:,1:max_ind);
    master_y = vertcat(master_y(1:2,:),master_y(4:end,:));
    master_xi = master_xi(1:max_ind);
    lambda = master_xi(end);
    denom_val = (master_y(2,:).^3/Fr^2-master_y(1,:).^2);
    crit_ind = sum(denom_val<0);
    exact = (master_xi(crit_ind)*((abs(denom_val(crit_ind))+abs(denom_val(crit_ind+1)))-abs(denom_val(crit_ind)))+master_xi(crit_ind+1)*(abs(denom_val(crit_ind))))/(abs(denom_val(crit_ind))+abs(denom_val(crit_ind+1)));
    crit_y = interp1(master_xi,master_y',exact)';
    y_in = horzcat(master_y(:,1:crit_ind),crit_y,crit_y,master_y(:,crit_ind+1:end));
    y_in(3,:) = y_in(3,:)/lambda*lambda_vis;
    xi_in = horzcat(master_xi(:,1:crit_ind)/exact,1,1,1+(master_xi(:,crit_ind+1:end)-exact)/(lambda-exact));
    p_in = [master_u_w,lambda,exact];
    rf_in = y_in(3,end);

    h_alt = y_in(2,1);
    tol=1e-4;
    denom_eps = 1e-4;
    solInit1=bvpinit(xi_in,@bvp_guess,p_in);
    opts = bvpset('RelTol',tol,'NMax',1000);

    solN1 = bvp4c(@(x,y,r,p) viscous_syst(x,y,r,p),@bc_vals,solInit1,opts);
    resid = solN1.stats.maxres;
    
    y_out = solN1.y;
    xi_out = solN1.x;
    p_out = solN1.parameters;
    out_vec = vertcat(xi_out,y_out);
    filename = "var_rho_master_pres_h.txt";
    save("Results/"+filename,"out_vec","-ascii")
    write_record("Results/full_record.csv",filename,{"var_rho","water",Fr,theta,alpha,d,tau0,p_out(1),p_out(2),p_out(3),0})
        
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

    function dydxi = viscous_syst(xi,y,region,p)
        % The system from the viroulet paper with a different friction
        % law. Need to solve for u_w and Q1 too.
        u_w = p(1);
        lambda = p(2);
        h_crit_pos = p(3);

        Q1 = y(1);
        h = y(2);
        u = u_w - Q1/h;
        m = y(3);
        phi = y(4)./Q1;

        rho_dl = phi*rho_p_dl+(1-phi)*rho_f_dl;
        p_tot_grad_dl = rho_dl*g_dl*cosd(theta);
        p_p_contact = max(a_dl*(phi-phi_rlp)/(phi_rcp-phi),0);
        Iv_phi = (phi_c/phi-1)^2;
        p_p_shear = 3*u*eta_f_dl/Iv_phi/h;
        p_p = (p_p_contact+p_p_shear);
        pb = p_tot_grad_dl*h-p_p;
        D = -2/beta_dl/h*(pb-h);
        Iv = abs(3*eta_f_dl*u/h/p_p);
        R_w3 = -phi*rho_f_dl/rho_dl*D;

        n_coeff = 1-Q1.^2.*Fr^2/h^3;
        mu_val = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv);
        fb_val = tand(theta)-sign(u).*mu_val;
        numer = fb_val+(u_w-u)*P*D;
        dQdxi = -P*D;
        dhdxi = (numer)./n_coeff;
        dmdxi = h*rho_dl/rho_eq_dl/lambda*u^(1-pres_h);
        dy4dxi = phi*rho_f_dl/rho_dl*D;
        dydxi = [dQdxi,dhdxi,dmdxi,dy4dxi];
        
        if abs(n_coeff)<denom_eps
            dppdh = 3*eta_f_dl/Iv_phi/h*(Q1/h^2)-3*u*eta_f_dl/Iv_phi/h^2;
            dppdQ = 3*eta_f_dl/Iv_phi/h*(-1/h);
            dppdphi = a_dl*(phi_rcp-phi_rlp)/(phi_rcp-phi)^2+6*u*eta_f_dl/h/(phi_c/phi-1)^3*phi_c/phi^2;
            dphidxi = (dy4dxi-phi*dQdxi)/Q1;
            dudh = Q1/h^2;
            dudQ = -1/h;
            
            fb_pp_coeff = -mu_Iv_fn(Iv)./(p_tot_grad_dl*h)+p_p./(p_tot_grad_dl*h).*dmudIv_fn(Iv).*Iv./p_p;
            full_pp_coeff = fb_pp_coeff+2*(u_w-u)/beta_dl*P/h;
%             fb_pb_coeff = -fb_pp_coeff; %-P.*2./beta_dl -P.*2./beta_dl
            full_u_coeff = -p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv)*Iv/u-P*D;
            fb_h_coeff = p_p./(p_tot_grad_dl.*h^2)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv/h;
            full_h_coeff = fb_h_coeff-(u_w-u)*P*(p_p/h^2)+full_pp_coeff*dppdh+full_u_coeff*dudh;
            fb_phi_coeff = p_p./(rho_dl^2.*g_dl*cosd(theta).*h^2)*mu_Iv_fn(Iv)*(rho_p_dl-rho_f_dl);
            full_phi_coeff = fb_phi_coeff+(u_w-u)*((rho_p_dl-rho_f_dl)*rho_f_dl/rho^2*D-P*(rho_p_dl-rho_f_dl)*g_dl*cosd(theta))+full_pp_coeff*dppdphi;
            full_Q_coeff = full_u_coeff*dudQ+full_pp_coeff*dppdQ;
            
%             full_h_coeff = (h.^2.*(full_fb_h_coeff+other_h_coeff)+2.*h.*(fb_val./Fr^2+P*D*(u_w-u)));
%                 full_h_coeff = (h.^3./Fr^2.*(g_dl*cosd(theta)*rho_dl.*chi.*A_pb_coeff+A_h_coeff-p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*Q1./h.^2)+3.*h.^2./Fr^2.*(fb_val+P.*D.*h));
%             numer_h_term = full_h_coeff.*dhdxi;
%             numer_other_term = h.^2.*((fb_pp_coeff*(dppdQ*dQdxi+dppdphi*dphidxi)+p_p./(p_tot_grad_dl).*dmudIv_fn(Iv).*Iv./u.*dQdxi./h)./Fr^2+P*D*dQdxi./h+u*P*((rho_dl-rho_f_dl)*g_dl*cosd(theta)*dphidxi-(dppdQ*dQdxi+dppdphi*dphidxi)));

            h_quad_roots = roots([3*Q1^2/Fr^2/h^4, -full_h_coeff-2*Q1/Fr^2/h^3*dQdxi, -full_phi_coeff*dphidxi-full_Q_coeff*dQdxi]);
            if imag(h_quad_roots(1))==0
                dydxi(2) = min(h_quad_roots(h_quad_roots>0));
            else
                dydxi(2) = 0;
            end
        end
        switch region
            case 1    % x in [0 1]
                dydxi = dydxi*h_crit_pos;
            case 2    % x in [1 2]
                dydxi = dydxi*(lambda-h_crit_pos);
        end
        dydxi = dydxi';
    end
        
    function resid = bc_vals(ya,yb,p)
        h_start = h_alt;
        h_stop = (-h_start + sqrt(h_start.^2+8*ya(1,1)^2./(h_start/Fr^2)))/2;
        
        h_mid = ya(2,2);
        u_mid = (-ya(1,2) + h_mid.*p(1))./h_mid;
        denom_mid = (h_mid.^3/Fr^2-ya(1,2).^2);
        phi_mid = ya(4,2)/ya(1,2);
        rho_mid = (rho_p_dl*phi_mid+rho_f_dl*(1-phi_mid));
        P_mid = (rho_mid-rho_f_dl)/rho_mid;
        p_tot_grad_mid = rho_mid*g_dl*cosd(theta);
        p_p_contact_mid = max(a_dl*(phi_mid-phi_rlp)/(phi_rcp-phi_mid),0);
        Iv_phi_mid = (phi_c/phi_mid-1)^2;
        p_p_shear_mid = 3*u_mid*eta_f_dl/Iv_phi_mid/h_mid;
        p_p_mid = (p_p_contact_mid+p_p_shear_mid);
        pb_mid = p_tot_grad_mid*h_mid-p_p_mid;
        D_mid = -2/beta_dl/h_mid*(pb_mid-h_mid);
        Iv_mid = 3*eta_f_dl*abs(u_mid)/h_mid/p_p_mid;
        fric_coeff_mid = p_p_mid/(p_tot_grad_mid*h_mid)*mu_Iv_fn(Iv_mid);
        fb_val_mid = tand(theta)-fric_coeff_mid-tau0_dl*rho_f_dl/rho_mid/h_mid+P_mid*D_mid*h_mid;
        
        cont_resid = (ya(:,2) - yb(:,1))';  %[ya(1,2)-yb(1,1),ya(2,2)-yb(2,1),ya(3,2)-yb(3,1),ya(4,2)-yb(4,1), ...
             %ya(8,2)-yb(8,1),ya(8,2)-yb(8,1),ya(8,2)-yb(8,1),ya(8,2)-yb(8,1)]
        resid = [ya(1,1)-yb(1,2), ya(2,1)-h_start, yb(2,2)-h_stop, ya(3,1), yb(3,2)-rf_in,...
            denom_mid, fb_val_mid, cont_resid]; 
    end
end