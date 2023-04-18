function [xi_final,y_final] = better_to_var_rho(custom_init,reverse,params,provide_init,master_xi,master_y)
    [phi_c,rho_f,rho_p,~,eta_f,g] = get_params_water();
%     P = (rho-rho_f)/rho;
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
        master_name = "no_vis_better_master.txt";
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/full_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        theta = record.theta(in_table);
        master_lambda = record.lambda(in_table);
        Fr = record.Fr(in_table);
        tau0 = record.tau0(in_table);
        master_u_w = record.u_w(in_table);
        master_crit_xi = record.crit_xi(in_table);
        
        stat_len = 0;
        if reverse
            alpha = record.alpha(in_table);
            d = record.d(in_table);
        else
            alpha = 1e-5;
            d = 1e-4;
        end
    end
    
    filename = "var_rho_master.txt";
    master_p = [master_u_w,master_lambda,master_crit_xi];

    
    
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    
    crit_pb = rho_f*g*cosd(theta)*h0;
    
%     p_tot = rho*g*cosd(theta);
    
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
    
%     p_tot_grad_dl = p_tot/p_scale*z_scale;
%     rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    
    
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    
    h_alt = master_y(2,1);
    y_in = master_y;
    xi_in=master_xi;
    
    tol=1e-4;
    denom_eps = 1e-4;
    
    nstep=100;
    [y_out,xi_out,p_out] = recurse_step(nstep,master_y,master_xi,master_p);
    out_vec = vertcat(xi_out,y_out);
    save("Results/"+filename,"out_vec","-ascii")
    write_record("Results/full_record.csv",filename,{"full","water",Fr,theta,alpha,d,tau0,p_out(1),p_out(2),p_out(3)})
    
    function [y_in,xi_in,p_in] = recurse_step(nstep,y_in,xi_in,p_in,curr,target,counter) %,
        if ~exist("curr","var")
            curr = 0;
        end
        if ~exist("target","var")
            target = 1;
        end
        if ~exist("counter","var")
            counter = 1;
        end
        if (counter>10)
            error("Max iteration depth reached")
        end
        if (counter == 1)
            phi_param_list = [0,logspace(-6,0,nstep-1)];
        else
            phi_param_list = linspace(curr,target,nstep);
        end
        n_min = 2-(counter==1);
        for i=n_min:nstep
            phi_param = phi_param_list(i);
            solInit1=bvpinit(xi_in,@bvp_guess,p_in);
            opts = bvpset('RelTol',tol,'NMax',2000);
            solN1 = bvp4c(@(x,y,r,p) viscous_syst(x,y,r,p),@bc_vals,solInit1,opts);
            resid = solN1.stats.maxres;
            if (resid>tol)
                [y_in,xi_in,p_in] = recurse_step(3,y_in,xi_in,p_in,phi_param_list(i-1),phi_param,counter+1);
            else
                y_in = solN1.y;
                xi_in = solN1.x;
                p_in = solN1.parameters;
            end 
        end
        
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

        function dydxi = viscous_syst(xi,y, region, p)
            % The system from the viroulet paper with a different friction
            % law. Need to solve for u_w and Q1 too.
            u_w = p(1);
            lambda = p(2);
            h_crit_pos = p(3);

            Q1 = y(1);
            h = y(2);
            u = (-Q1 + h.*u_w)./h;
            m = y(3);
            phi = y(4)./Q1;
            rho_dl = phi_param*(rho_p_dl*phi+rho_f_dl*(1-phi))+(1-phi_param)*(rho_p_dl*phi_c+rho_f_dl*(1-phi_c));
            P = (rho_dl-rho_f_dl)/rho_dl;
            chi = (rho_f_dl+3*rho_dl)/(4*rho_dl);
            pb = y(5) + g_dl*cosd(theta)*rho_dl*chi.*h;

            p_tot_grad_dl = rho_dl*g_dl*cosd(theta);
            zeta = 3/(2*alpha_dl*h) + P/4;
            p_p = p_tot_grad_dl*h-pb;
            Iv = 3*eta_f_dl*abs(u)/h/p_p;
            D = -2/beta_dl/h*(pb-h);

            denom = (h.^3/Fr^2-Q1.^2);
            fric_coeff = p_p/(p_tot_grad_dl*h)*mu_Iv_fn(Iv);
            fb_val = tand(theta)-fric_coeff-tau0_dl*rho_f_dl/rho_dl/h;
            dhdxi = 1/Fr^2.*h.^3.*(fb_val+P*D*h)./denom;

            R_w3 = -phi*rho_f_dl/rho_dl*D;
            R_w4 = ((-P/4)*(phi_param)-P*chi*(1-phi_param)+zeta)*D - 2*3/alpha_dl/h*u*(phi - phi_c./(1+sqrt(Iv)));

            dQdxi = -P*D;
            dmdxi = h/(lambda+stat_len)*u;
            dy6dxi = -R_w3;
            dy7dxi = R_w4/(u-u_w);
            dydxi = [dQdxi,dhdxi,dmdxi,dy6dxi,dy7dxi];

    %         dpbdxi = dy7dxi + g_dl*cosd(theta)*rho_dl*chi.*dhdxi;
    %         dppdxi = p_tot_grad_dl*dhdxi-dpbdxi;
    %         dudxi = -dQdxi./h+Q1/h^2*dhdxi;
    %         dIvdxi = 3*eta_f_dl*dudxi/h/p_p-3*eta_f_dl*abs(u)/h^2/p_p*dhdxi-3*eta_f_dl*abs(u)/h/p_p^2*dppdxi;
    %         ddenomdxi = (3*h.^2/Fr^2.*dhdxi-2*Q1.*dQdxi);
    %         dnumerdxi = (-dppdxi/(p_tot_grad_dl*h)+p_p/(p_tot_grad_dl*h^2)*dhdxi)*mu_Iv_fn(Iv)-p_p/(p_tot_grad_dl*h)*dmudIv_fn(Iv)*dIvdxi...
    %             +tau0_dl*rho_f/rho/h^2*dhdxi+P*2/beta_dl*(dpbdxi-dhdxi);
    %         d2hdxi2 = 3/Fr^2.*h.^2.*(fb_val+P*D*h)./denom.*dhdxi+(dnumerdxi*denom-(fb_val+P*D*h)*ddenomdxi)/denom^2;
            if abs(denom)<denom_eps
                A_pp_coeff = -mu_Iv_fn(Iv)./(p_tot_grad_dl.*h)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./p_p;
                A_pb_coeff = -P.*2./beta_dl-A_pp_coeff;
                A_h_coeff = p_p./(p_tot_grad_dl.*h.^2).*mu_Iv_fn(Iv)+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*3*eta_f_dl.*abs(u)./h.^2./p_p+tau0_dl.*rho_f_dl./rho_dl./h.^2-P.*2./beta_dl+A_pp_coeff.*rho_dl.*g_dl.*cosd(theta);

                full_h_coeff = (h.^3./Fr^2.*(g_dl*cosd(theta)*rho_dl.*chi.*A_pb_coeff+A_h_coeff-p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*Q1./h.^2)+3.*h.^2./Fr^2.*(fb_val+P.*D.*h));
                numer_h_term = full_h_coeff.*dhdxi;
                numer_other_term = h.^3./Fr^2.*(dy7dxi.*A_pb_coeff+p_p./(p_tot_grad_dl.*h).*dmudIv_fn(Iv).*Iv./u.*dQdxi./h.^2);

                h_quad_roots = roots([3*h.^2./Fr^2, -full_h_coeff-2.*Q1.*dQdxi, -numer_other_term]);
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
            phi_mid = ya(4,2);
            denom_mid = (h_mid.^3/Fr^2-ya(1,2).^2);
            rho_mid = phi_param*(rho_p_dl*phi_mid+rho_f_dl*(1-phi_mid))+(1-phi_param)*(rho_p_dl*phi_c+rho_f_dl*(1-phi_c));
            P_mid = (rho_mid-rho_f_dl)/rho_mid;
            chi_mid = (rho_f_dl+3*rho_mid)/(4*rho_mid);
            pb_mid = ya(5,2) + g_dl*cosd(theta)*rho_mid*chi_mid.*h_mid;
            p_tot_grad_dl_mid = rho_mid*g_dl*cosd(theta);
            p_p_mid = p_tot_grad_dl_mid*h_mid-pb_mid;
            D_mid = -2/beta_dl/h_mid*(pb_mid-h_mid);
            Iv_mid = 3*eta_f_dl*abs(u_mid)/h_mid/p_p_mid;
            fric_coeff_mid = p_p_mid/(p_tot_grad_dl_mid*h_mid)*mu_Iv_fn(Iv_mid);
            fb_val_mid = tand(theta)-fric_coeff_mid-tau0_dl*rho_f_dl/rho_mid/h_mid+P_mid*D_mid*h_mid;

            cont_resid = (ya(:,2) - yb(:,1))';
            resid = [ya(1,1)-yb(1,2), ya(2,1)-h_start, yb(2,2)-h_stop, ya(3,1), yb(3,2)-1,...
                (ya(5,1)-yb(5,2)), denom_mid, fb_val_mid, cont_resid]; 
        end
    end
end