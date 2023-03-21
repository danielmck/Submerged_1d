function [xi_final,y_final] = no_pe_full_depo_sinful(custom_init,reverse,params,provide_init,master_xi,master_y)
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
        if ~provide_init
            if reverse
                [master_xi,master_y] = bvp_full_from_master(true,params(1:4));
            else
                [master_xi,master_y] = viscous_Iv_bvp_from_master(true,params(1:4));
            end
        end
    else
        master_name = "no_vis_no_pe_break.txt";
        master_file = load("Results/"+master_name);
        master_xi = master_file(1,:);
        master_y = master_file(2:end,:);
        record = readtable('Results/wave_record.csv');

        in_table = strcmp(record.Name, master_name);
        wave_type = record.wave_type(in_table);
        theta = record.theta(in_table);
        h_alt = record.h_alt(in_table);
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
    
    [h0, crit_Iv] = crit_Iv_tau0(theta, rho_p, rho_f, eta_f, Fr, tau0);
    u_eq = Fr*sqrt(g*cosd(theta)*h0);
    phi_eq = phi_c/(1+sqrt(crit_Iv));
    p_tot = rho*g*cosd(theta);
    crit_pb = rho_f*g*cosd(theta)*h0;
    
    z_scale = h0;
    v_scale = u_eq;
    p_scale = crit_pb;
    t_scale = z_scale/v_scale;

    tau0_dl = tau0/p_scale;
    eta_f_dl = eta_f/(p_scale*t_scale);
    alpha_dl = alpha*p_scale;
    g_dl = g*t_scale/v_scale; 
    
    u_eq_dl = u_eq/v_scale;
    p_tot_grad_dl = p_tot/p_scale*z_scale;
    rho_f_dl = rho_f*v_scale^2/p_scale;
    rho_p_dl = rho_p*v_scale^2/p_scale; 
    rho_dl = rho_p_dl*phi_c+rho_f_dl*(1-phi_c);
    d_dl = d/z_scale;
    
    beta_dl = 150*phi_c.^2.*eta_f_dl./((1-phi_c).^3.*d_dl^2);
    chi = (rho_f+3*rho)/(4*rho);
    
    delta_vals = logspace(-6,0,100);
    delta_vals = horzcat(0,delta_vals);
    if reverse
        delta_vals = wrev(delta_vals(2:end));
    end
%     stable = viscous_stability(theta,Fr,nu,lambda);
    h_pert = 1e-4;
    if ~reverse
        Q1_in = master_y(4,:);
        h_in_l = master_y(5,:);
        h_in_r = master_y(7,:);
        y6_in_l= phi_eq*Q1_in;
        y6_in_r= phi_eq*Q1_in;
        
        pb_st = 0.85*h_in_l(1);
        pb_end = pb_st+rho_dl*g_dl*cosd(theta)*chi*(h_in_r(end)-h_in_l(1));
        pb_in_l = pb_st+(pb_end-pb_st)*(h_in_l-h_in_l(1))/(h_in_r(end)-h_in_l(1));
        pb_in_r = pb_st+(pb_end-pb_st)*(h_in_r-h_in_l(1))/(h_in_r(end)-h_in_l(1));
        
        y7_in_l = zeros(size(master_xi)); %pb_in_l - rho_dl*g_dl*cosd(theta)*chi.*h_in_l;
        y7_in_r = zeros(size(master_xi)); %pb_in_r - rho_dl*g_dl*cosd(theta)*chi.*h_in_r;
        
        y_in = vertcat(master_y(1:6,:),y6_in_l,y7_in_l,master_y(4,:),master_y(7:8,:),y6_in_r,y7_in_r);
        xi_in=master_xi;
%         master_y = vertcat(master_y(1,:),master_y(3,:),master_y(2,:),master_y(4:end,:));
%         y6_in = phi_eq*master_y(3,:);
%         pb_st = 0.85*master_y(4,1);
%         pb_end = pb_st+rho_dl*g_dl*cosd(theta)*chi*(master_y(4,end)-master_y(4,1));
%         pb_in = pb_st+(pb_end-pb_st)*(master_y(4,:)-master_y(4,1))/(master_y(4,end)-master_y(4,1));
%         y7_in = master_y(4,:) - rho_dl*g_dl*cosd(theta)*chi.*master_y(4,:);
%         y_full = vertcat(master_y,y6_in,y7_in);
%         
%         xi_crit = interp1(master_y(4,:),master_xi,(master_y(3,1)*Fr)^(2/3));
%         y_l_in = interp1(master_xi,y_full',linspace(0,xi_crit,100))';
%         y_l_in(4,end)= (master_y(3,1)*Fr)^(2/3)-h_pert;
%         y_r_in = interp1(master_xi,y_full',linspace(xi_crit,1,100))';
%         y_r_in(4,1) = (master_y(3,1)*Fr)^(2/3)+h_pert;
%         y_in = vertcat(y_l_in(1:2,:),xi_crit*y_l_in(2,1)*ones(1,100),y_l_in(3:end,:),y_r_in(3:end,:));
%         xi_in = linspace(0,1,100);
%         y_in=master_y;

        
    else
        y_in = master_y;
    end
    [xi_out,y_out] = run_bvp_iter(delta_vals, xi_in,y_in);
    if reverse
        no_pe_param = [Fr,theta,h_alt,tau0];
        [xi_out,y_out] = viscous_Iv_bvp_from_master(true,no_pe_param,true,xi_out,y_out(1:5,:),no_pe_param);
    end
    if ~custom_init && ~reverse && ~provide_init
        out_vec = vertcat(xi_out,y_out);
        filename = "no_vis_full_master.txt";
        save("Results/"+filename,"out_vec","-ascii")
        write_record("Results/wave_record.csv",filename,{"full_no_vis","water",Fr,theta,h_alt,alpha,d,tau0})
    end
    
    function [xi_in,y_in] = run_bvp_iter(delta_vals, xi_in, y_in, tol,counter)
        if ~exist('counter','var')
            counter = 1;
        end
        if ~exist('tol','var')
            tol = 1e-3;
        end
        if counter > 10
            error("Max iteration depth reached, non convergence")
        end
        n_step = size(delta_vals,2);
        for i = 2:n_step
            delta_in = delta_vals(i);
            solInit1=bvpinit(xi_in,@bvp_guess);
            opts = bvpset('RelTol',tol,'NMax',1000);
            try
                solN1 = bvp4c(@viscous_syst,@bc_vals,solInit1,opts);
                resid = solN1.stats.maxres;
            catch ME
                switch ME.identifier
                    case 'MATLAB:bvp4c:SingJac'
                        warning('Singular Jacobian encountered, reducing step size');
                        resid = tol+1;
                    otherwise
                        rethrow(ME)
                end    
            end
            if resid < tol
                y_in = solN1.y;
                xi_in = solN1.x;
            else
                [xi_in,y_in] = run_bvp_iter(linspace(delta_vals(i-1),delta_in,3)...
                    ,xi_in, y_in, tol,counter+1);
            end
        end
        
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
            pb_l = y(8) + g_dl*cosd(theta)*(rho_f_dl+(rho_dl*chi-rho_f_dl)*delta_in).*h_l;
            
            zeta_l = 3/(2*alpha_dl*h_l) + P/4;
            p_p_l = p_tot_grad_dl*h_l-pb_l;
            Iv_l = 3*eta_f_dl*abs(u_l)/h_l/p_p_l;
            D_l = -2/beta_dl/h_l*(pb_l-h_l);
            
            denom_l = (h_l.^3/Fr^2-Q1_l.^2);
            fric_coeff_l = delta_in*p_p_l/(p_tot_grad_dl*h_l)*mu_Iv_fn(Iv_l)+(1-delta_in)*P*mu_Iv_fn(crit_Iv*abs(u_l)/h_l^2);
            fb_val_l = tand(theta)-fric_coeff_l-tau0_dl*rho_f/rho/h_l;
            dhdxi_l = 1/Fr^2.*h_l.^3.*fb_val_l./denom_l;
            
            R_w3_l = -phi_l*rho_f_dl/rho_dl*D_l;
            R_w4_l = (-P*chi+zeta_l)*D_l - 2*3/alpha_dl/h_l*u_l*(phi_l - phi_c./(1+sqrt(Iv_l)));
            
            dQdxi_l = -delta_in*P*D_l;
            dmdxi_l = h_l/(lambda+stat_len)*u_l;
            dy6dxi_l = -delta_in*R_w3_l;
            dy7dxi_l = delta_in*R_w4_l/(u_l-u_w);
            dydxi_l = [dQdxi_l,dhdxi_l,dmdxi_l,dy6dxi_l,dy7dxi_l]*h_crit_pos;
            
            Q1_r = y(9);
            h_r = y(10);
            u_r = (-Q1_r + h_r.*u_w)./h_r;
            m_r = y(11);
            phi_r = y(12)./Q1_r;
            pb_r = y(13) + g_dl*cosd(theta)*(rho_f_dl+(rho_dl*chi-rho_f_dl)*delta_in).*h_r;
            
            zeta_r = 3/(2*alpha_dl*h_r) + P/4;
            p_p_r = p_tot_grad_dl*h_r-pb_r;
            Iv_r = 3*eta_f_dl*abs(u_r)/h_r/p_p_r;
            D_r = -2/beta_dl/h_r*(pb_r-h_r);
            
            denom_r = (h_r.^3/Fr^2-Q1_r.^2);
            fric_coeff_r = delta_in*p_p_r/(p_tot_grad_dl*h_r)*mu_Iv_fn(Iv_r)+(1-delta_in)*P*mu_Iv_fn(crit_Iv*abs(u_r)/h_r^2);
            fb_val_r = tand(theta)-fric_coeff_r-tau0_dl*rho_f/rho/h_r;
            dhdxi_r = 1/Fr^2.*h_r.^3.*fb_val_r./denom_r;
            
            R_w3_r = -phi_r*rho_f_dl/rho_dl*D_r;
            R_w4_r = (-P*chi+zeta_r)*D_r - 2*3/alpha_dl/h_r*u_r*(phi_r - phi_c./(1+sqrt(Iv_r)));
            
            dQdxi_r = -delta_in*P*D_r;
            dmdxi_r = h_r/(h_crit_pos+stat_len)*u_r;
            dy6dxi_r = -delta_in*R_w3_r;
            dy7dxi_r = delta_in*R_w4_r/(u_r-u_w);
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
end